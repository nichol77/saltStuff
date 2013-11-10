
Int_t depths[3][10];
Int_t disks[3][10][2];
Int_t files[3][10][2];
Double_t minTime[3][10][2];
Double_t maxTime[3][10][2];
Double_t peakTime[3][10][2];
Int_t numDepths[3]={9,3,3};

float dh1=163*.3048;//in meters
float dh2=553*.3048;//in meters

Double_t timeChange[10][1000];
Double_t attenLengthArray[10][1000];


Double_t antFreqCentres[3]={150,125,200};
Double_t antFreqWidth[3]={50,50,50};
Int_t antFreqNumBins[3]={9,2,17};

Double_t antDtLeft[3]={-10e-9,-15e-9,-5e-9};
Double_t antDtRight[3]={27e-9,23e-9,7e-9};

Int_t antT13Amp[3]={0,0,26};

Double_t funcPowerLaw(Double_t *x, Double_t *par) {
  return par[0]*TMath::Power(x[0]*1e6,par[1]);
}

void plotAbbyAtten()
{
  gSystem->Load("libMathMore.so");
  gSystem->Load("/sw/lib/libfftw3.so");
  gSystem->Load("libRootFftwWrapper.so");
  readAbbyTimeFile();
  char canName[180];
  char graphTitle[180];
  //  TCanvas *canAtten = new TCanvas("canAtten","canAtten",600,400);
  TGraph *grAtten[3][10];
  int antIndex=2;

  for(int i=0;i<numDepths[antIndex];i++) {
    sprintf(canName,"canCrop%d",depths[antIndex][i]);
    cout << depths[antIndex][i] << "\t" << disks[antIndex][i][0] << "\t" << files[antIndex][i][0] << endl;
    cout << depths[antIndex][i] << "\t" << disks[antIndex][i][1] << "\t" << files[antIndex][i][1] << endl;
    //    cout << depths[antIndex][i] << "\t" << peakTime[antIndex][i][0] << "\t" << peakTime[antIndex][i][1] << endl;

    //    cout << minTime[antIndex][i][1] << "\t" << maxTime[antIndex][i][1] << endl;
    

    TGraph *grT12=getWave(disks[antIndex][i][0],files[antIndex][i][0]);
    TGraph *grT13=getWave(disks[antIndex][i][1],files[antIndex][i][1]);
    
    int count=0;
    TGraph *grT12Crop = cropWave(grT12,1e-9*peakTime[antIndex][i][0]+antDtLeft[antIndex],1e-9*peakTime[antIndex][i][0]+antDtRight[antIndex]);
    TGraph *grT13Crop = cropWave(grT13,1e-9*peakTime[antIndex][i][1]+antDtLeft[antIndex],1e-9*peakTime[antIndex][i][1]+antDtRight[antIndex]);

    
    TGraph *grPower12 = FFTtools::makePowerSpectrumVoltsSeconds(grT12Crop);
    TGraph *grPower13 = FFTtools::makePowerSpectrumVoltsSeconds(grT13Crop);
    TGraph *grPower12dB = convertToDb(grPower12);
    TGraph *grPower13dB = convertToDb(grPower13);
    sprintf(graphTitle,"Hole 1 - Hole 2 -- Depth %d (Christian)",depths[antIndex][i]);
    grT12->SetTitle(graphTitle);
    grPower12->SetTitle(graphTitle);
    grPower12dB->SetTitle(graphTitle);
    
    sprintf(graphTitle,"Hole 1 - Hole 3 -- Depth %d (Christian)",depths[antIndex][i]);
    grT13->SetTitle(graphTitle);
    grPower13->SetTitle(graphTitle);
    grPower13dB->SetTitle(graphTitle);
    
    Double_t power12=FFTtools::integrateVoltageSquared(grT12Crop);
    Double_t power13=FFTtools::integrateVoltageSquared(grT13Crop);

    
    TGraphErrors *grPowInt12 = integratePowerSpectrum(grPower12,antFreqCentres[antIndex],antFreqWidth[antIndex],antFreqNumBins[antIndex]);
    TGraphErrors *grPowInt13 = integratePowerSpectrum(grPower13,antFreqCentres[antIndex],antFreqWidth[antIndex],antFreqNumBins[antIndex]);
    if(antT13Amp[antIndex]!=0) {
      //Need to take account of amplification
      Int_t numPows=grPowInt13->GetN();
      Double_t *powVals=grPowInt13->GetY();
      Double_t *freqVals=grPowInt13->GetX();
      Double_t *powErrs=grPowInt13->GetEY();
      Double_t *freqErrs=grPowInt13->GetEX();
      Double_t gainFactor=TMath::Power(10,-1*Double_t(antT13Amp[antIndex]/10));
      for(int pi=0;pi<numPows;pi++) {
	powVals[pi]*=gainFactor;
      }
      grPowInt13 = new TGraphErrors(numPows,freqVals,powVals,freqErrs,powErrs);
    }


    sprintf(graphTitle,"grAtten%d",depths[antIndex][i]);
    grAtten[antIndex][i]=getAttenLengthGraph(grPowInt12,grPowInt13);
    grAtten[antIndex][i]->SetName(graphTitle);
    //      cout << depths[antIndex][i] << "\t" << power12 << "\t" << power13 << endl;
    //    Double_t *t12=grT12->GetX();
    //    Double_t *t13=grT13->GetX();
    //    cout << depths[antIndex][i] << "\t" << t12[1]-t12[0] << "\t" << t13[1]-t13[0]
    //	 << endl;
    
    Double_t ratio=dh1*dh1*power12/(dh2*dh2*power13);
    Double_t attenLength=2*(dh2-dh1)/TMath::Log(ratio);
    cout << depths[antIndex][i] << "\t" << attenLength << endl;
    //      timeChange[antIndex][i][count]=dt;
    //      attenLengthArray[antIndex][i][count]=attenLength;
      //      count++;
      //    }
    //    cout << i  << endl;
    //    grAtten[antIndex][i]= new TGraph(count,timeChange[antIndex][i],attenLengthArray[antIndex][i]);
    //   }
    //   TH1F *framey = gPad->DrawFrame(-10e-9,0,10e-9,200);
    //   For(int i=0;i<9;i++) {
    //     grAtten[antIndex][i]->SetLineColor(getNiceColour(i));
    //     grAtten[antIndex][i]->Draw("l");
    //   }
    //   {
    //     return;
    //     exit(0);
    
    TCanvas *canCrop = new TCanvas(canName,canName,800,600);
    canCrop->Divide(3,2);
    canCrop->cd(1);
    grT12->SetLineColor(getNiceColour(1));
    grT12->Draw("al");
    grT12->GetXaxis()->SetRangeUser(1e-9*peakTime[antIndex][i][0]-20e-9,1e-9*peakTime[antIndex][i][0]+60e-9);
    grT12Crop->SetLineColor(getNiceColour(0));
    grT12Crop->Draw("l");
    canCrop->cd(2);
    grPower12dB->SetLineColor(getNiceColour(0));
    grPower12dB->Draw("al");
    canCrop->cd(3);
    grPowInt12->SetLineColor(getNiceColour(0));
    grPowInt12->Draw("alp");

    canCrop->cd(4);
    grT13->SetLineColor(getNiceColour(1));
    grT13->Draw("al");
    grT13->GetXaxis()->SetRangeUser(1e-9*peakTime[antIndex][i][1]-20e-9,1e-9*peakTime[antIndex][i][1]+60e-9);
    grT13Crop->SetLineColor(getNiceColour(0));
    grT13Crop->Draw("l");
    canCrop->cd(5);
    grPower13dB->SetLineColor(getNiceColour(0));
    grPower13dB->Draw("al");
    canCrop->cd(6);
    grPowInt13->SetLineColor(getNiceColour(0));
    grPowInt13->Draw("al");
  }
  TF1 *fpl = new TF1("fpl",funcPowerLaw,10,1000,2);
  fpl->SetParameters(3.281e6,-0.5484);
  TCanvas *canAtten = new TCanvas("canAtten","canAtten");
  gPad->SetLogy(1);
  gPad->SetLogx(1);
  TH1F *framey = gPad->DrawFrame(10,10,1000,1000);
  for(int i=1;i<numDepths[antIndex];i++) {
    grAtten[antIndex][i]->SetLineColor(getNiceColour(i));
    grAtten[antIndex][i]->SetMarkerColor(getNiceColour(i));
    grAtten[antIndex][i]->Draw("lp");
  }
  fpl->Draw("same");
			  

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
  Int_t depthIndex[3]={-1,-1,-1};
  Int_t lastDepth[3]={0,0,0};

  while(line.ReadLine(CsvFile)) {

    //    cout << line.Data() << endl;
    TObjArray *vals = (TObjArray*) line.Tokenize(",");
    //    vals->Dump();
    TObjString *firstField=vals->At(0);
    //    exit(0);
    TObjString *antenna=vals->At(1);
    TObjString *disk=vals->At(2);
    TObjString *file=vals->At(3);
    TObjString *hole=vals->At(6);
    TObjString *peakTimeStr=vals->At(7);
    TObjString *minTimeStr=vals->At(8);
    TObjString *maxTimeStr=vals->At(9);
    //    cout << vals[0].Data() << "\t" << vals[1].Data() << endl;
    //    cout << x->GetString()->Data() << "\t" << y->GetString()->Data() << endl;
    Int_t antIndex=antenna->GetString().Atoi()-1;
    Int_t holeIndex=hole->GetString().Atoi()-2;
    Int_t depth=firstField->GetString().Atoi();
    if(lastDepth[antIndex]!=depth) {
      depthIndex[antIndex]++;
      depths[antIndex][depthIndex[antIndex]]=depth;
      lastDepth[antIndex]=depth;
    }
    //    cout << antIndex << "\t" << depthIndex[antIndex] << "\t" << holeIndex << "\t" << depth << endl;
    //    cout << numLines << "\t" << index <<  "\t" << depths[depthIndex[antIndex]] << endl;
    disks[antIndex][depthIndex[antIndex]][holeIndex]=disk->GetString().Atoi()+1;
    files[antIndex][depthIndex[antIndex]][holeIndex]=file->GetString().Atoi();
    peakTime[antIndex][depthIndex[antIndex]][holeIndex]=peakTimeStr->GetString().Atof();
    minTime[antIndex][depthIndex[antIndex]][holeIndex]=minTimeStr->GetString().Atof();
    maxTime[antIndex][depthIndex[antIndex]][holeIndex]=maxTimeStr->GetString().Atof();

    numLines++;
    
  }

}
  
TGraphErrors *integratePowerSpectrum(TGraph *grIn, Double_t firstBinCentre, Double_t freqBinWidth, Int_t numBins) {  
  //  cout << numBins << endl;
  Double_t *minFreq= new Double_t [numBins];
  Double_t *maxFreq= new Double_t [numBins];
  Double_t *centreFreqs = new Double_t [numBins];
  Double_t *powerInt = new Double_t [numBins];
  Double_t *powerIntSq = new Double_t [numBins];
  Double_t *error = new Double_t [numBins];
  for(int i=0;i<numBins;i++) {
    minFreq[i]=firstBinCentre+(i-0.5)*freqBinWidth;
    centreFreqs[i]=firstBinCentre+i*freqBinWidth;
    maxFreq[i]=firstBinCentre+(i+0.5)*freqBinWidth;  
    powerInt[i]=0;
    powerIntSq[i]=0;
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
	//	cout << power[i] << endl;
  	powerInt[bin]+=deltaF*power[i]*1e20;
	powerIntSq[bin]+=(deltaF*power[i]*1e20)*(deltaF*power[i]*1e20);
      }
    }
  }
  for(int i=0;i<numBins;i++) {
    error[i]=0;
    //    error[i]=TMath::Sqrt(powerIntSq[i]-(powerInt[i]*powerInt[i]));
    //    cout << powerInt[i]*powerInt[i] << "\t" << powerIntSq[i]-(powerInt[i]*powerInt[i]) << endl;
  }
  TGraphErrors *grInt = new TGraphErrors(numBins,centreFreqs,powerInt,0,error);
  delete [] minFreq;
  delete [] maxFreq;
  delete [] centreFreqs;
  delete [] powerInt;
  delete [] powerIntSq;
  delete [] error;
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
