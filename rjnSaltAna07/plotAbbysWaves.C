
Int_t depths[10];
Int_t disks[10][2];
Int_t files[10][2];
Double_t minTime[10][2];
Double_t maxTime[10][2];
Double_t peakTime[10][2];


float dh1=163*.3048;//in meters
float dh2=553*.3048;//in meters

Double_t timeChange[10][1000];
Double_t attenLengthArray[10][1000];


void plotAbbysWaves()
{
  gSystem->Load("libMathMore.so");
  gSystem->Load("/sw/lib/libfftw3.so");
  gSystem->Load("libRootFftwWrapper.so");
  readAbbyTimeFile();
  char canName[180];
  char graphTitle[180];
  //  TCanvas *canAtten = new TCanvas("canAtten","canAtten",600,400);
  TGraph *grAtten[10];
  for(int i=0;i<9;i++) {
    sprintf(canName,"canCrop%d",depths[i]);
    //    cout << disks[i][0] << "\t" << files[i][0] << endl;
    //    cout << depths[i] << "\t" << peakTime[i][0] << "\t" << peakTime[i][1] << endl;

    //    cout << minTime[i][1] << "\t" << maxTime[i][1] << endl;
    

    TGraph *grT12=getWave(disks[i][0],files[i][0]);
    TGraph *grT13=getWave(disks[i][1],files[i][1]);
    
    int count=0;
    //    for(Double_t dt=-10e-9;dt<10e-9;dt+=0.5e-9) {
    TGraph *grT12Crop = cropWave(grT12,1e-9*peakTime[i][0]-10e-9,1e-9*peakTime[i][0]+27e-9);
    TGraph *grT13Crop = cropWave(grT13,1e-9*peakTime[i][1]-10e-9,1e-9*peakTime[i][1]+27e-9);

    //    TGraph *grT12Crop = cropWave(grT12,1e-9*minTime[i][0],1e-9*maxTime[i][0]);
    //    TGraph *grT13Crop = cropWave(grT13,1e-9*minTime[i][1],1e-9*maxTime[i][1]);
    
    TGraph *grPower12Raw = FFTtools::makePowerSpectrumVoltsSeconds(grT12Crop);
    TGraph *grPower13Raw = FFTtools::makePowerSpectrumVoltsSeconds(grT13Crop);
    TGraph *grPower12=convertToDb(grPower12Raw);
    TGraph *grPower13=convertToDb(grPower13Raw);
    sprintf(graphTitle,"Hole 1 - Hole 2 -- Depth %d (Christian)",depths[i]);
    grT12->SetTitle(graphTitle);
    grPower12->SetTitle(graphTitle);
    
    sprintf(graphTitle,"Hole 1 - Hole 3 -- Depth %d (Christian)",depths[i]);
    grT13->SetTitle(graphTitle);
    grPower13->SetTitle(graphTitle);
    
    Double_t power12=FFTtools::integrateVoltageSquared(grT12Crop);
    Double_t power13=FFTtools::integrateVoltageSquared(grT13Crop);

    TGraph *grPowInt12 = integratePowerSpectrum(grPower12,100,7);
    TGraph *grPowInt13 = integratePowerSpectrum(grPower13,100,7);
    sprintf(graphTitle,"grAtten%d",depths[i]);
    grAtten[i]=getAttenLengthGraph(grPowInt12,grPowInt13);
    grAtten[i]->SetName(graphTitle);
    //      cout << depths[i] << "\t" << power12 << "\t" << power13 << endl;
    //    Double_t *t12=grT12->GetX();
    //    Double_t *t13=grT13->GetX();
    //    cout << depths[i] << "\t" << t12[1]-t12[0] << "\t" << t13[1]-t13[0]
    //	 << endl;
    
    Double_t ratio=dh1*dh1*power12/(dh2*dh2*power13);
    Double_t attenLength=2*(dh2-dh1)/TMath::Log(ratio);
    cout << depths[i] << "\t" << attenLength << endl;
    //      timeChange[i][count]=dt;
    //      attenLengthArray[i][count]=attenLength;
      //      count++;
      //    }
    //    cout << i  << endl;
    //    grAtten[i]= new TGraph(count,timeChange[i],attenLengthArray[i]);
    //   }
    //   TH1F *framey = gPad->DrawFrame(-10e-9,0,10e-9,200);
    //   For(int i=0;i<9;i++) {
    //     grAtten[i]->SetLineColor(getNiceColour(i));
    //     grAtten[i]->Draw("l");
    //   }
    //   {
    //     return;
    //     exit(0);
    
    TCanvas *canCrop = new TCanvas(canName,canName,800,600);
    canCrop->Divide(3,2);
    canCrop->cd(1);
    grT12->SetLineColor(getNiceColour(1));
    grT12->Draw("al");
    grT12->GetXaxis()->SetRangeUser(1e-9*peakTime[i][0]-20e-9,1e-9*peakTime[i][0]+60e-9);
    grT12Crop->SetLineColor(getNiceColour(0));
    grT12Crop->Draw("l");
    canCrop->cd(2);
    grPower12->SetLineColor(getNiceColour(0));
    grPower12->Draw("al");
    canCrop->cd(3);
    grPowInt12->SetLineColor(getNiceColour(0));
    grPowInt12->Draw("al");

    canCrop->cd(4);
    grT13->SetLineColor(getNiceColour(1));
    grT13->Draw("al");
    grT13->GetXaxis()->SetRangeUser(1e-9*peakTime[i][1]-20e-9,1e-9*peakTime[i][1]+60e-9);
    grT13Crop->SetLineColor(getNiceColour(0));
    grT13Crop->Draw("l");
    canCrop->cd(5);
    grPower13->SetLineColor(getNiceColour(0));
    grPower13->Draw("al");
    canCrop->cd(6);
    grPowInt13->SetLineColor(getNiceColour(0));
    grPowInt13->Draw("al");
  }
  
  TCanvas *canAtten = new TCanvas("canAtten","canAtten");
  gPad->SetLogy(1);
  gPad->SetLogx(1);
  TH1F *framey = gPad->DrawFrame(10,10,1000,1000);
  for(int i=1;i<9;i++) {
    grAtten[i]->SetLineColor(getNiceColour(i));
    grAtten[i]->SetMarkerColor(getNiceColour(i));
    grAtten[i]->Draw("lp");
  }
			  

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
  //  cout << fileName << endl;
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
  //  cout << grWave << endl;
  return grWave;
}
  


TGraph *getWave(int diskNum, int fileNum) 
{
  char filename[180];
  sprintf(filename,"/home/rjn/saltStuff/salt_trip3_data/disk_%d/TEK%05d.CSV",diskNum,fileNum);
  return getWave(filename);
}


void readAbbyTimeFile()
{
  ifstream CsvFile("/home/rjn/saltStuff/abbyTimeWindows.csv");
  
  int numLines=0;
  TString line;
  
  line.ReadLine(CsvFile); // Get rid of inital headings
  
  while(line.ReadLine(CsvFile)) {

    //    cout << line.Data() << endl;
    TObjArray *vals = (TObjArray*) line.Tokenize(",");
    //    vals->Dump();
    TObjString *firstField=vals->At(0);
    //    exit(0);

    TObjString *disk=vals->At(1);
    TObjString *file=vals->At(2);
    TObjString *peakTimeStr=vals->At(6);
    TObjString *minTimeStr=vals->At(7);
    TObjString *maxTimeStr=vals->At(8);
    //    cout << vals[0].Data() << "\t" << vals[1].Data() << endl;
    //    cout << x->GetString()->Data() << "\t" << y->GetString()->Data() << endl;
    Int_t index=(firstField->GetString().Atoi()/10)-1;   
    depths[index]=firstField->GetString().Atoi();
    //    cout << numLines << "\t" << index <<  "\t" << depths[index] << endl;
    if(numLines<9) {
      disks[index][0]=disk->GetString().Atoi()+1;
      files[index][0]=file->GetString().Atoi();
      peakTime[index][0]=peakTimeStr->GetString().Atof();
      minTime[index][0]=minTimeStr->GetString().Atof();
      maxTime[index][0]=maxTimeStr->GetString().Atof();
    }
    else {
      disks[index][1]=disk->GetString().Atoi()+1;
      files[index][1]=file->GetString().Atoi();
      peakTime[index][1]=peakTimeStr->GetString().Atof();
      minTime[index][1]=minTimeStr->GetString().Atof();
      maxTime[index][1]=maxTimeStr->GetString().Atof();
    }

    numLines++;
    
  }

}
  
TGraph *integratePowerSpectrum(TGraph *grIn, Double_t freqBinWidth, Int_t numBins) {  
  //  cout << numBins << endl;
  Double_t *minFreq= new Double_t [numBins];
  Double_t *maxFreq= new Double_t [numBins];
  Double_t *centreFreqs = new Double_t [numBins];
  Double_t *powerInt = new Double_t [numBins];
  for(int i=0;i<numBins;i++) {
    minFreq[i]=i*freqBinWidth;
    centreFreqs[i]=(i+0.5)*freqBinWidth;
    maxFreq[i]=(i+1)*freqBinWidth;  
    powerInt[i]=0;
  }
  
  Double_t *freqs = grIn->GetX();
  Double_t *power = grIn->GetY();
  Int_t numPoints=grIn->GetN();
  Double_t deltaF=freqs[1]-freqs[0];
  //  cout << numPoints << "\t" << deltaF << endl;
  for(int i=0;i<numPoints;i++) {
    //    cout << i << "\t" << freqs[i]<< endl;
    for(int bin=0;bin<numBins;bin++) {      
      if(freqs[i]>=minFreq[bin] && freqs[i]<maxFreq[bin]) {
  	powerInt[bin]+=deltaF*power[i];
      }
    }
  }
  TGraph *grInt = new TGraph(numBins,centreFreqs,powerInt);
  return grInt;
} 

TGraph *getAttenLengthGraph(TGraph *gr12,TGraph *gr13)
{
  Int_t numPoints=gr12->GetN();
  Double_t *freq=gr12->GetX();
  Double_t *pow12=gr12->GetY();
  Double_t *pow13=gr13->GetY();

  Double_t *attenLength= new Double_t [numPoints];
  for(int i=0;i<numPoints;i++) {
    //    cout << freq[i] << "\t" << pow13[i] << "\t" << pow12[i] << endl;
    Double_t rat=(dh1*dh1*pow12[i])/(dh2*dh2*pow13[i]);
    attenLength[i]=2*(dh2-dh1)/TMath::Log(rat);
    //    cout << freq[i] << "\t" << attenLength[i] << endl;
  }
  
  TGraph *grAtten = new TGraph(numPoints,freq,attenLength);
  return grAtten;
}

TGraph *convertToDb(TGraph *grIn) {
  Double_t *freqs=grIn->GetX();
  Double_t *linPow =grIn->GetY();
  Int_t numPoints=grIn->GetN();
  Double_t *dbPow = new Double_t [numPoints];

  for(int i=0;i<numPoints;i++) {
    dbPow[i]=10*TMath::Log10(linPow[i]);
  }

  TGraph *grDb = new TGraph(numPoints,freqs,dbPow);
  delete [] dbPow;
  return grDb;
}
