
TGraph *getWave(char *fileName);

void plotWaveAndFFT()
{
  gSystem->Load("/sw/lib/libfftw3.so");
  gSystem->Load("libRootFftwWrapper.so");
  plotWaveAndFFT("/home/rjn/saltStuff/salt_trip3_data/above_ground_data/TEK00000.CSV");
}


void plotWaveAndFFT(char *fileName)
{
  TGraph *gr=getWave(fileName);
  TCanvas *can = new TCanvas("can","can");
  can->Divide(1,3);
  can->cd(1);
  gr->Draw("al");
  can->cd(2);
  TGraph *grPower=FFTtools::makePowerSpectrumVoltsSeconds(gr);
  grPower->Draw("al");
  can->cd(3);  
  TGraph *grPeriod=FFTtools::makePowerSpectrumPeriodogram(gr);
  grPeriod->Draw("al");
  
}

TGraph *getWave(char *fileName)
{
  ifstream CsvFile(fileName);
  
  int numLines=0;
  int numLines2=0;
  TString line;
  Double_t xVals[100000];
  Double_t yVals[100000];
  Double_t xVals2[100000];
  Double_t yVals2[100000];
  while(line.ReadLine(CsvFile)) {

    //    cout << line.Data() << endl;
    TObjArray *vals = (TObjArray*) line.Tokenize(",");
    //    vals.Dump();
    TObjString *x=vals->At(0);
    TObjString *y=vals->At(1);
    //    cout << vals[0].Data() << "\t" << vals[1].Data() << endl;
    //    cout << x->GetString()->Data() << "\t" << y->GetString()->Data() << endl;
    xVals[numLines]=x->GetString().Atof();
    yVals[numLines]=y->GetString().Atof();
    //    yVals[numLines]=TMath::Sin(TMath::TwoPi()*500e6*xVals[numLines]);
    if(numLines%2==0) {
      xVals2[numLines2]=xVals[numLines];
      yVals2[numLines2]=yVals[numLines];
      numLines2++;
    }
    numLines++;
    
  }
  TGraph *grWave = new TGraph(numLines,xVals,yVals);
  return grWave;
}
  
