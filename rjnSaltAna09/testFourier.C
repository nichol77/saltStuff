
TGraph *getWave(char *fileName);
Int_t integratePower(TGraph *gr,Int_t firstBin=-1,Int_t lastBin=-1);
Int_t integratePowerProper(TGraph *gr,Int_t firstBin=-1,Int_t lastBin=-1);
Double_t sumVoltageSquared(TGraph *gr,Int_t firstBin=-1,Int_t lastBin=-1);
Double_t integrateVoltageSquared(TGraph *gr,Int_t firstBin=-1,Int_t lastBin=-1);

void testFourier()
{
  gSystem->Load("/sw/lib/libfftw3.so");
  gSystem->Load("libRootFftwWrapper.so");
  doTest();
}

void doTest()
{

  TGraph *gr1=getWave("/home/rjn/saltStuff/salt_trip3_data/disk_4/TEK00024.CSV");
  TGraph *gr2=getWave("/home/rjn/saltStuff/salt_trip3_data/disk_4/TEK00025.CSV");


  Int_t N1=gr1->GetN();
  Double_t *x1=gr1->GetX();  
  Double_t dt1=x1[1]-x1[0];
  Int_t peakBin1=FFTtools::getPeakBin(gr1);
  //  cout << dt1  << "\t" << N1*dt1 << "\n";
 
  Int_t N2=gr2->GetN();
  Double_t *x2=gr2->GetX();  
  Double_t *y2=gr2->GetY();  
  Double_t dt2=x2[1]-x2[0];
  Int_t peakBin2=FFTtools::getPeakBin(gr2);
  //  cout << dt2  << "\t" << N2*dt2 << "\n";
  
  Double_t firstPoint1=x1[0]-x1[peakBin1];
  Double_t lastPoint1=x1[N1-1]-x1[peakBin1];
  
  Int_t newFirstBin2=Int_t(peakBin2+(firstPoint1/dt2));
  Int_t newN2=Int_t(TMath::Ceil(peakBin2+lastPoint1/dt2))-newFirstBin2;
  TGraph *gr2Zoom = new TGraph(newN2,&(x2[newFirstBin2]),&(y2[newFirstBin2]));
  

  //  cout <<  << "\t" << TMath::Ceil(peakBin2+lastPoint1/dt2) << "\n";

  TCanvas *canWave = new TCanvas("canWave","canWave");
  canWave->Divide(1,2);
  canWave->cd(1);
  gr1->SetTitle("Example Pulse -- 251 Samples");
  gr1->Draw("al");  
  gr1->SetLineColor(getNiceColour(0));
  canWave->cd(2);
  gr2Zoom->SetTitle("Example Pulse -- 126 Samples");
  gr2Zoom->Draw("al");
  gr2Zoom->SetLineColor(getNiceColour(1));

  TGraph *grPower1=FFTtools::makePowerSpectrum(gr1);
  TGraph *grPower2=FFTtools::makePowerSpectrum(gr2Zoom);
  
  TCanvas *canPower = new TCanvas("canPower","canPower");
  canPower->Divide(1,2);
  canPower->cd(1);
  grPower1->SetTitle("Unnormalised -- 251 Samples");
  grPower1->Draw("al");  
  grPower1->SetLineColor(getNiceColour(0));
  canPower->cd(2);
  grPower2->SetTitle("Unnormalised -- 126 Samples");
  grPower2->Draw("al");		
  grPower2->SetLineColor(getNiceColour(1));
  
  cout << sumVoltageSquared(gr1) << "\t" << sumVoltageSquared(gr2Zoom) << "\n";
  cout << "Unnormalised\t" << integratePower(grPower1) << "\t" << integratePower(grPower2) << "\n";


  TGraph *grPeriod1=FFTtools::makePowerSpectrumPeriodogram(gr1);
  TGraph *grPeriod2=FFTtools::makePowerSpectrumPeriodogram(gr2Zoom);
  
  TCanvas *canPeriod = new TCanvas("canPeriod","canPeriod");
  canPeriod->Divide(1,2);
  canPeriod->cd(1);
  grPeriod1->Draw("al"); 
  grPeriod1->SetTitle("Peiodogram -- 251 Samples"); 
  grPeriod1->SetLineColor(getNiceColour(0));
  canPeriod->cd(2);
  grPeriod2->Draw("al");		
  grPeriod2->SetTitle("Peiodogram -- 126 Samples"); 
  grPeriod2->SetLineColor(getNiceColour(1));
  
  cout << grPeriod1->GetN() << "\t" << grPeriod2->GetN() << "\n";
  cout << sumVoltageSquared(gr1)/gr1->GetN() << "\t" << sumVoltageSquared(gr2Zoom)/gr2Zoom->GetN() << "\n";
  cout << integrateVoltageSquared(gr1) << "\t" << integrateVoltageSquared(gr2Zoom) << "\n";
  cout << "Periodogram\t" << integratePower(grPeriod1) << "\t" << integratePower(grPeriod2) << "\n";


  TGraph *grdB1=FFTtools::makePowerSpectrumVoltsSeconds(gr1);
  TGraph *grdB2=FFTtools::makePowerSpectrumVoltsSeconds(gr2Zoom);
  
  TCanvas *candB = new TCanvas("candB","candB");
  candB->Divide(1,2);
  candB->cd(1);
  grdB1->Draw("al");  
  grdB1->SetTitle("Normalised  -- 251 Samples"); 
  grdB1->SetLineColor(getNiceColour(0));
  candB->cd(2);
  grdB2->SetTitle("Normalised  -- 126 Samples"); 
  grdB2->Draw("al");		
  grdB2->SetLineColor(getNiceColour(1));
  //  cout << integratePower

  TGraph *grdBPadded1=FFTtools::makePowerSpectrumVoltsSecondsPadded(gr1);
  TGraph *grdBPadded2=FFTtools::makePowerSpectrumVoltsSecondsPadded(gr2Zoom);
  
  TCanvas *candBPadded = new TCanvas("candBPadded","candBPadded");
  candBPadded->Divide(1,2);
  candBPadded->cd(1);
  grdBPadded1->Draw("al"); 
  grdBPadded1->SetTitle("Normalised -- 251 Samples (x4 Zero Padding)");  
  grdBPadded1->SetLineColor(getNiceColour(0));
  candBPadded->cd(2);
  grdBPadded2->SetTitle("Normalised -- 126 Samples (x4 Zero Padding)"); 
  grdBPadded2->Draw("al");		
  grdBPadded2->SetLineColor(getNiceColour(1));

  cout << "Sum\n";
  cout << integratePower(grdB1) << "\t" << integratePower(grdB2) << "\n";
  cout << integratePower(grdBPadded1) << "\t"  << integratePower(grdBPadded2) << "\n";
  cout << "Integral\n";
  cout << integratePowerProper(grdB1) << "\t" << integratePowerProper(grdB2) << "\n";
  cout << integratePowerProper(grdBPadded1) << "\t"  << integratePowerProper(grdBPadded2) << "\n";

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
    yVals[numLines]=y->GetString().Atof()+100;
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
  
