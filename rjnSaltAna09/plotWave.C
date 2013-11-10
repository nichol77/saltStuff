
void plotWave()
{
  TGraph *grT12=getWave(3,2);
  TGraph *grT12Crop = cropWave(grT12,-10e-9,24.70253e-9);
  TGraph *grT13=getWave(4,0);
  TGraph *grT13Crop = cropWave(grT13,3.7067e-9,40e-9);

  TCanvas *canCrop = new TCanvas("canCrop","canCrop",600,600);
  canCrop->Divide(1,2);
  canCrop->cd(1);
  grT12->SetLineColor(getNiceColour(1));
  grT12->Draw("al");
  grT12Crop->SetLineColor(getNiceColour(0));
  grT12Crop->Draw("l");
  canCrop->cd(2);
  grT13->SetLineColor(getNiceColour(1));
  grT13->Draw("al");
  grT13Crop->SetLineColor(getNiceColour(0));
  grT13Crop->Draw("l");
}

TGraph *cropWave(TGraph *grIn, Double_t minTime, Double_t maxTime)
{
  Int_t numPoints=grIn->GetN();
  Double_t *times = grIn->GetX();
  Double_t *oldVolts= grIn->GetY();
  Double_t *newVolts = new Double_t[numPoints];

  for(int i=0;i<numPoints;i++) {
    if(times[i]>=minTime && times[i]<=maxTime)
      newVolts[i]=oldVolts[i];
    else
      newVolts[i]=0;
  }
  TGraph *grCrop = new TGraph(numPoints,times,newVolts);
  delete [] newVolts;
  return grCrop;
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
  


TGraph *getWave(int diskNum, int fileNum) 
{
  char filename[180];
  sprintf(filename,"/home/rjn/saltStuff/salt_trip3_data/disk_%d/TEK%05d.CSV",diskNum,fileNum);
  cout << filename << endl;
  getWave(filename);
}
