
Int_t depths[3][10];
Int_t disks[3][10][2];
Int_t files[3][10][2];
Double_t minTime[3][10][2];
Double_t maxTime[3][10][2];
Double_t peakTime[3][10][2];
Int_t numDepths[3]={9,3,3};
Int_t canDivide[3][2]={{3,3},{1,3},{1,3}};


float dh1=163*.3048;//in meters
float dh2=553*.3048;//in meters

Double_t timeChange[10][1000];
Double_t attenLengthArray[10][1000];

Double_t antFreqCentres[3]={150,125,200};
Double_t antFreqWidth[3]={50,50,100};
Int_t antFreqNumBins[3]={9,2,9};


Double_t antDtLeft[3]={-10e-9,-15e-9,-5e-9};
Double_t antDtRight[3]={27e-9,23e-9,7e-9};

char *antName[3]={"christian","long","baby"};

Double_t funcPowerLaw(Double_t *x, Double_t *par) {
  return par[0]*TMath::Power(x[0]*1e6,par[1]);
}

void plotAbbyTimeThings()
{
  gSystem->Load("libMathMore.so");
  gSystem->Load("/sw/lib/libfftw3.so");
  gSystem->Load("libRootFftwWrapper.so");
  readAbbyTimeFile();
  char canName[180];
  char graphTitle[180];

  TF1 *fpl = new TF1("fpl",funcPowerLaw,10,1000,2);
  fpl->SetParameters(3.281e6,-0.5484);
  TGraph *grAtten[1000];

  int antIndex=2; 

  char outName[180];
  sprintf(outName,"%sErrorTimeSlide.txt",antName[antIndex]);
  ofstream AbbyNumbers(outName);
  sprintf(outName,"%sErrorTimeSlideArrays.C",antName[antIndex]);
  ofstream AbbyArrays(outName);

  TCanvas *canAtten = new TCanvas("canAtten","canAtten",1000,800);
  TCanvas *canPower12 = new TCanvas("canPower12","canPower12",1000,800);
  TCanvas *canPower13 = new TCanvas("canPower13","canPower13",1000,800);
  canAtten->Divide(canDivide[antIndex][0],canDivide[antIndex][1]);
  canPower12->cd();
  gPad->SetLogy();
  TH1F *frame12 = canPower12->DrawFrame(10,1e4,1200,1e11);
  canPower13->cd();
  gPad->SetLogy();
  TH1F *frame13 = canPower13->DrawFrame(10,100,1200,1e8);
  for(int i=0;i<numDepths[antIndex];i++) {
    canAtten->cd(i+1);
    sprintf(canName,"Depth %d",depths[antIndex][i]);
    gPad->SetLogy(1);
    gPad->SetLogx(1);
    TH1F *framey = gPad->DrawFrame(30,10,1000,1000,canName);

    Double_t meanPower12[100]={0};
    Double_t meanSqPower12[100]={0};
    Double_t meanPower13[100]={0};
    Double_t meanSqPower13[100]={0};
    Double_t freqVal[100]={0};
    Int_t numPowPoints=0;


    TGraph *grT12=getWave(disks[antIndex][i][0],files[antIndex][i][0]);
    TGraph *grT13=getWave(disks[antIndex][i][1],files[antIndex][i][1]);

    int count=0;

    for(Double_t dt=-3e-9;dt<3e-9;dt+=0.25e-9) {
      TGraph *grT12Crop = cropWave(grT12,(1e-9*peakTime[antIndex][i][0]+antDtLeft[antIndex])+dt,1e-9*peakTime[antIndex][i][0]+antDtRight[antIndex]+dt);
      TGraph *grT13Crop = cropWave(grT13,(1e-9*peakTime[antIndex][i][1]+antDtLeft[antIndex])+dt,1e-9*peakTime[antIndex][i][1]+antDtRight[antIndex]+dt);

      //      TGraph *grT12Crop = cropWave(grT12,1e-9*minTime[antIndex][i][0]+dt,1e-9*maxTime[antIndex][i][0]+dt);
      //      TGraph *grT13Crop = cropWave(grT13,1e-9*minTime[antIndex][i][1]+dt,1e-9*maxTime[antIndex][i][1]+dt);
      
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
      

      sprintf(graphTitle,"grAtten%d_%d",depths[antIndex][i],count);
      grAtten[count]=getAttenLengthGraph(grPowInt12,grPowInt13);
      grAtten[count]->SetName(graphTitle);
      //      grAtten[count]->SetLineColor(getNiceColour(count));
      //      grAtten[count]->SetMarkerColor(getNiceColour(count));
      grAtten[count]->SetLineColor(count);
      grAtten[count]->SetMarkerColor(count);
      grAtten[count]->Draw("lp");
      fpl->Draw("same");

      numPowPoints=grPowInt13->GetN();
      Double_t *pow12Vals=grPowInt12->GetY();
      Double_t *pow13Vals=grPowInt13->GetY();
      Double_t *freq13Vals=grPowInt13->GetX();
      for(int powi=0;powi<numPowPoints;powi++) {
	meanPower12[powi]+=pow12Vals[powi];
	meanSqPower12[powi]+=pow12Vals[powi]*pow12Vals[powi];
	meanPower13[powi]+=pow13Vals[powi];
	meanSqPower13[powi]+=pow13Vals[powi]*pow13Vals[powi];
	freqVal[powi]=freq13Vals[powi];
      }

      count++;
    }
    //Now just slide the 12 window
 for(Double_t dt=-3e-9;dt<3e-9;dt+=0.25e-9) {
      TGraph *grT12Crop = cropWave(grT12,(1e-9*peakTime[antIndex][i][0]+antDtLeft[antIndex])+dt,1e-9*peakTime[antIndex][i][0]+antDtRight[antIndex]+dt);
      TGraph *grT13Crop = cropWave(grT13,(1e-9*peakTime[antIndex][i][1]+antDtLeft[antIndex]),1e-9*peakTime[antIndex][i][1]+antDtRight[antIndex]);

      //      TGraph *grT12Crop = cropWave(grT12,1e-9*minTime[antIndex][i][0]+dt,1e-9*maxTime[antIndex][i][0]+dt);
      //      TGraph *grT13Crop = cropWave(grT13,1e-9*minTime[antIndex][i][1]+dt,1e-9*maxTime[antIndex][i][1]+dt);
      
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
      

      sprintf(graphTitle,"grAtten%d_%d",depths[antIndex][i],count);
      grAtten[count]=getAttenLengthGraph(grPowInt12,grPowInt13);
      grAtten[count]->SetName(graphTitle);
      //      grAtten[count]->SetLineColor(getNiceColour(count));
      //      grAtten[count]->SetMarkerColor(getNiceColour(count));
      grAtten[count]->SetLineColor(count);
      grAtten[count]->SetMarkerColor(count);
      grAtten[count]->Draw("lp");
      fpl->Draw("same");

      numPowPoints=grPowInt13->GetN();
      Double_t *pow12Vals=grPowInt12->GetY();
      Double_t *pow13Vals=grPowInt13->GetY();
      Double_t *freq13Vals=grPowInt13->GetX();
      for(int powi=0;powi<numPowPoints;powi++) {
	meanPower12[powi]+=pow12Vals[powi];
	meanSqPower12[powi]+=pow12Vals[powi]*pow12Vals[powi];
	meanPower13[powi]+=pow13Vals[powi];
	meanSqPower13[powi]+=pow13Vals[powi]*pow13Vals[powi];
	freqVal[powi]=freq13Vals[powi];
      }

      count++;
    }
 cout << count << endl;
 //Now just slide the 13 window
 for(Double_t dt=-3e-9;dt<3e-9;dt+=0.25e-9) {
      TGraph *grT12Crop = cropWave(grT12,(1e-9*peakTime[antIndex][i][0]+antDtLeft[antIndex]),1e-9*peakTime[antIndex][i][0]+antDtRight[antIndex]);
      TGraph *grT13Crop = cropWave(grT13,(1e-9*peakTime[antIndex][i][1]+antDtLeft[antIndex])+dt,1e-9*peakTime[antIndex][i][1]+antDtRight[antIndex]+dt);

      //      TGraph *grT12Crop = cropWave(grT12,1e-9*minTime[antIndex][i][0]+dt,1e-9*maxTime[antIndex][i][0]+dt);
      //      TGraph *grT13Crop = cropWave(grT13,1e-9*minTime[antIndex][i][1]+dt,1e-9*maxTime[antIndex][i][1]+dt);
      
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
      

      sprintf(graphTitle,"grAtten%d_%d",depths[antIndex][i],count);
      grAtten[count]=getAttenLengthGraph(grPowInt12,grPowInt13);
      grAtten[count]->SetName(graphTitle);
      //      grAtten[count]->SetLineColor(getNiceColour(count));
      //      grAtten[count]->SetMarkerColor(getNiceColour(count));
      grAtten[count]->SetLineColor(count);
      grAtten[count]->SetMarkerColor(count);
      grAtten[count]->Draw("lp");
      fpl->Draw("same");

      numPowPoints=grPowInt13->GetN();
      Double_t *pow12Vals=grPowInt12->GetY();
      Double_t *pow13Vals=grPowInt13->GetY();
      Double_t *freq13Vals=grPowInt13->GetX();
      for(int powi=0;powi<numPowPoints;powi++) {
	meanPower12[powi]+=pow12Vals[powi];
	meanSqPower12[powi]+=pow12Vals[powi]*pow12Vals[powi];
	meanPower13[powi]+=pow13Vals[powi];
	meanSqPower13[powi]+=pow13Vals[powi]*pow13Vals[powi];
	freqVal[powi]=freq13Vals[powi];
      }

      count++;
    }

    Double_t rms12Val[100];
    Double_t rms13Val[100];
    Double_t err12Val[100];
    Double_t err13Val[100];
    for(int powi=0;powi<numPowPoints;powi++) {
      meanPower12[powi]/=count;
      meanSqPower12[powi]/=count;
      meanPower13[powi]/=count;
      meanSqPower13[powi]/=count;
      rms12Val[powi]=TMath::Sqrt(meanSqPower12[powi]-meanPower12[powi]*meanPower12[powi]);
      rms13Val[powi]=TMath::Sqrt(meanSqPower13[powi]-meanPower13[powi]*meanPower13[powi]);
      err12Val[powi]=rms12Val[powi]/TMath::Sqrt(count);
      err13Val[powi]=rms13Val[powi]/TMath::Sqrt(count);

    }
    AbbyNumbers << antIndex << "\t" << depths[antIndex][i] << "\t" << 2 << "\t";
    AbbyArrays << "meanPowerArray_" << antIndex << "_" << depths[antIndex][i] << "_" << 2 << "[" << numPowPoints << "]={";
    for(int powi=0;powi<numPowPoints;powi++) {
      AbbyNumbers << meanPower12[powi] << "\t" << rms12Val[i] << "\t";
      AbbyArrays << meanPower12[powi];
      if(powi<numPowPoints-1) 
	AbbyArrays << ",";
      else 
	AbbyArrays << "};\n";
    }
    AbbyArrays << "rmsPowerArray_" << antIndex << "_" << depths[antIndex][i] << "_" << 2 << "[" << numPowPoints << "]={";
     for(int powi=0;powi<numPowPoints;powi++) {
      AbbyArrays << rms12Val[powi];
      if(powi<numPowPoints-1) 
	AbbyArrays << ",";
      else 
	AbbyArrays << "};\n";
    }
    AbbyArrays << "meanPowerArray_" << antIndex << "_" << depths[antIndex][i] << "_" << 3 << "[" << numPowPoints << "]={";
    for(int powi=0;powi<numPowPoints;powi++) {
      AbbyArrays << meanPower13[powi];
      if(powi<numPowPoints-1) 
	AbbyArrays << ",";
      else 
	AbbyArrays << "};\n";
    }
    AbbyArrays << "rmsPowerArray_" << antIndex << "_" << depths[antIndex][i] << "_" << 3 << "[" << numPowPoints << "]={";
     for(int powi=0;powi<numPowPoints;powi++) {
      AbbyArrays << rms13Val[powi];
      if(powi<numPowPoints-1) 
	AbbyArrays << ",";
      else 
	AbbyArrays << "};\n";
    }
        
    AbbyNumbers << "\n";
    AbbyNumbers << antIndex << "\t" << depths[antIndex][i] << "\t" << 3 << "\t";
    for(int powi=0;powi<numPowPoints;powi++) {
      AbbyNumbers << meanPower13[powi] << "\t" << rms13Val[i] << "\t";
    }
    AbbyNumbers << "\n";
    
    
    TGraphErrors *grPowErr12 = new TGraphErrors(numPowPoints,freqVal,meanPower12,0,rms12Val);
    TGraphErrors *grPowErr13 = new TGraphErrors(numPowPoints,freqVal,meanPower13,0,rms13Val);
    canPower12->cd();
    grPowErr12->SetLineColor(getNiceColour(i));
    grPowErr12->SetMarkerColor(getNiceColour(i));    
    grPowErr12->Draw("lp");
    canPower13->cd();
    grPowErr13->SetLineColor(getNiceColour(i));
    grPowErr13->SetMarkerColor(getNiceColour(i));
    grPowErr13->Draw("lp");
    
    
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
