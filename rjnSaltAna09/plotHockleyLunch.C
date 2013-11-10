
const Double_t lmr200Velocity=TMath::C()*0.83;
const Double_t lmr200NsPerFt=1.22;
const Double_t lmr600Velocity=TMath::C()*0.87;
const Double_t lmr600NsPerFt=1.17;
const Double_t lmr240Velocity=TMath::C()*0.84;
const Double_t lmr240NsPerFt=1.21;
const Double_t rg58Velocity=TMath::C()*0.66;
const Double_t rg58NsPerFt=0.97;
const Double_t syncLag=31.3;
const Double_t cps3InsertionDelayNs=20; //Should look up
const Double_t distBToHole1=12; //In m
const Double_t saltRefractiveIndex=2.413; //Take value from hole data

Double_t holeOneTimes[4]={-15e-9,-15e-9,-3e-9,-3e-9};

Char_t *antNames[4]={"Skinny Christian","Short Raven","Long Pole","Baby"};

#define WINDOW_BEFORE_PEAK 20e-9  //ns
#define WINDOW_AFTER_PEAK 30e-9  //ns

typedef enum EAntType {
  kSkinnyChristian = 0,
  kShortRaven = 1,
  kLongPole = 2,
  kBaby = 3 
} AntType_t;


#define NUM_FILES 51
Int_t diskNums[NUM_FILES];
Int_t fileNums[NUM_FILES];
Int_t rxLocation[NUM_FILES];
Int_t txLocation[NUM_FILES];
Int_t antNum[NUM_FILES];
Int_t ampInDb[NUM_FILES];


TGraph *getWave(char *fileName);
void plotWaveAndFFT(char *fileName, char *graphTitle=0);
void readInputFile();
Double_t getTriggerTime();
Double_t getRFCableDelay();
void plotAllLunchStuff();
void plotThisOne(Int_t diskNum, Int_t fileNum);
void plotThisOne(Int_t index);
char *getAntName(EAntType antType);
Int_t getWindowedPeakBin(TGraph *gr, Double_t low, Double_t high);
Double_t integratePowerSpectrum(TGraph *gr, Double_t low, Double_t high);
Double_t attenLengthLogDist(Double_t *x, Double_t *par);
Double_t refractiveIndex(Double_t *x, Double_t *par);
Double_t getLMR600Atten(Double_t freqInMHz, Double_t cableLength);
Double_t getLMR240Atten(Double_t freqInMHz, Double_t cableLengthInFt);

void plotThisOne(Int_t index) 
{

  char graphTitle[180];
  char fileName[180];
  sprintf(fileName,"/home/rjn/saltStuff/hockley2009/disk_%d/TEK%05d.txt",diskNums[index],fileNums[index]);
  sprintf(graphTitle,"Ant %d, Hole %d",antNum[index],rxLocation[index]);
  plotWaveAndFFT(fileName,graphTitle);
}

void plotThisOne(Int_t diskNum, Int_t fileNum) 
{
  char graphTitle[180];
  char fileName[180];
  sprintf(fileName,"/home/rjn/saltStuff/hockley2009/disk_%d/TEK%05d.txt",diskNum,fileNum);
  //  cout << fileName  << endl;
  sprintf(graphTitle,"Disk %d, File %d",diskNum,fileNum);
  plotWaveAndFFT(fileName,graphTitle);

}


void plotHockleyLunch()
{
  gSystem->Load("/sw/lib/libfftw3.so");
  gSystem->Load("libRootFftwWrapper.so");
  readInputFile();
  plotAllLunchStuff();
}


void plotWaveAndFFT(char *fileName, char *graphTitle)
{
  TGraph *gr=getWave(fileName);
  TCanvas *can = new TCanvas("can","can");
  can->Divide(1,3);
  can->cd(1);
  if(graphTitle) gr->SetTitle(graphTitle);
  gr->SetLineColor(getNiceColour(0));
  gr->Draw("al");
  can->cd(2);
  TGraph *grHilbert = FFTtools::getHilbertEnvelope(gr);
  if(graphTitle) grHilbert->SetTitle(graphTitle);
  grHilbert->SetLineColor(getNiceColour(0));
  grHilbert->Draw("al");
  Int_t peakBin=FFTtools::getPeakBin(grHilbert);
  Double_t *timeVals=grHilbert->GetX();
    cout << timeVals[peakBin] << "\t" << getPeakToPeak(gr)*1000 << "\n";
  can->cd(3);
  TGraph *grPower=FFTtools::makePowerSpectrumVoltsSecondsdB(gr);
  if(graphTitle) grPower->SetTitle(graphTitle);
  grPower->SetLineColor(getNiceColour(0));
  grPower->Draw("al");

}

TGraph *getWave(char *fileName)
{
  ifstream TxtFile(fileName);
  
  int numLines=0;
  int numLines2=0;
  TString line;
  Double_t xVals[1000000];
  Double_t yVals[1000000];
  Double_t xVals2[1000000];
  Double_t yVals2[100000];
  while(line.ReadLine(TxtFile)) {

    //    cout << line.Data() << endl;
    TObjArray *vals = (TObjArray*) line.Tokenize(" ");
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
  
void readInputFile()
{

  ifstream TxtFile("hockleyLunchAreaFileLog.csv");  
  int numLines=0;
  int numLines2=0;
  TString line;
  line.ReadLine(TxtFile);
  while(line.ReadLine(TxtFile)) {
    TObjArray *vals = (TObjArray*) line.Tokenize(",");
    //    vals.Dump();
    TObjString *tempDisk=vals->At(0);
    TObjString *tempFile=vals->At(1);
    TObjString *tempAnt=vals->At(2);
    TObjString *tempTx=vals->At(3);
    TObjString *tempRx=vals->At(4);

    diskNums[numLines]=tempDisk->GetString().Atoi();
    fileNums[numLines]=tempFile->GetString().Atoi();
    rxLocation[numLines]=tempRx->GetString().Atoi();
    TString tempString=tempAnt->GetString();
    if(tempString.Contains("SC")) {
      antNum[numLines]=kSkinnyChristian;
    }
    if(tempString.Contains("SR")) {
      antNum[numLines]=kShortRaven;
    }
    if(tempString.Contains("LP")) {
      antNum[numLines]=kLongPole;
    }
    if(tempString.Contains("B")) {
      antNum[numLines]=kBaby;
    }
    
    TString tempString2=tempTx->GetString();
    if(tempString2.Contains("A")) {
      txLocation[numLines]=1;
    }
    if(tempString2.Contains("B")) {
      txLocation[numLines]=2;
    }
    //    std::cout << rxLocation[numLines] << "\t" << numLines << "\n";
    numLines++;

  }
}


Double_t getTriggerTime()
{
  return syncLag+200*lmr200NsPerFt;


}

Double_t getRFCableDelay()
{
  Double_t rg58Delay=3*rg58NsPerFt;
  Double_t lmr240Delay=10*lmr240NsPerFt;
  Double_t lmr600Delay=100*lmr600NsPerFt;
  return rg58Delay+cps3InsertionDelayNs+lmr240Delay+lmr600Delay;
}

void plotAllLunchStuff()
{
  Double_t pulseDists[NUM_FILES];
  Double_t rxNumber[NUM_FILES];
  Double_t rawPeakTime[NUM_FILES];
  Double_t correctedPeakTime[NUM_FILES];
  Double_t peakToPeak[NUM_FILES];
  Double_t distTimesPeakToPeakOverTwo[NUM_FILES];
  Double_t distTimesVoltageRatio[NUM_FILES];

  Int_t antCount[4]={0};
  Int_t antFreqCount[4][8]={0};
  Double_t antPeakToPeak[4][NUM_FILES];
  Double_t antLogHalfPeakToPeakTimesDist[4][NUM_FILES];
  Double_t antLogSqrtPowerTimesDist[4][NUM_FILES];
  Double_t antLogSqrtIntPowerFreqTimesDist[4][NUM_FILES];
  Double_t antLogSqrtIntPowerFreqTimesDistFreqBins[4][8][NUM_FILES];
  Double_t antHoleDist[4][NUM_FILES];

  Int_t openDisk=5;
  Int_t openFile=0;
  char openName[180];
  sprintf(openName,"/home/rjn/saltStuff/hockley2009/disk_%d/TEK%05d.txt",openDisk,openFile);
  TGraph *grOpen =getWave(openName);
  TGraph *grOpenHilbert=FFTtools::getHilbertEnvelope(grOpen);
  Int_t peakBinOpen=FFTtools::getPeakBin(grOpenHilbert);
  Double_t *timeValsOpen=grOpenHilbert->GetX();
  Double_t openPeakTime=timeValsOpen[peakBinOpen];
  TGraph *grOpenCrop = FFTtools::cropWave(grOpen,openPeakTime-WINDOW_BEFORE_PEAK,openPeakTime+WINDOW_AFTER_PEAK);
  //  new TCanvas();
  //  grOpenCrop->Draw("al");
  new TCanvas();
  grOpenCrop->Draw("al");
  grOpenCrop->GetXaxis()->SetTitle("Time (s)");
  grOpenCrop->GetYaxis()->SetTitle("Voltgae (v)");
  Double_t openPeakToPeak=getPeakToPeak(grOpenCrop);
  Double_t totalAttendB=50+getLMR240Atten(150,15)+getLMR600Atten(150,50);
  //  std::cout << openPeakToPeak*0.5 << "\t" << totalAttendB << "\n";
  Double_t transmittedVoltage=0.5*openPeakToPeak*TMath::Power(10,totalAttendB/20.);
  std::cout << openPeakToPeak*0.5 << "\t" << totalAttendB << "\t" << transmittedVoltage << "\n";
  for(int index=0;index<NUM_FILES;index++) {    
    char fileName[180];
    sprintf(fileName,"/home/rjn/saltStuff/hockley2009/disk_%d/TEK%05d.txt",diskNums[index],fileNums[index]);
    //    std::cout << index << "\t" << fileName << "\n";
    TGraph *gr=getWave(fileName);
    TGraph *grHilbert = FFTtools::getHilbertEnvelope(gr);
    //    Int_t peakBin=getWindowedPeakBin(grHilbert,windowLow[rxLocation[index]][antNum[index]],windowHigh[rxLocation[index]][antNum[index]]);
    Int_t peakBin=FFTtools::getPeakBin(grHilbert);
    Double_t *timeVals=grHilbert->GetX();
    rawPeakTime[index]=timeVals[peakBin];
    correctedPeakTime[index]=(timeVals[peakBin]*1e9)+getTriggerTime()-getRFCableDelay()
;
    
    TGraph *grCrop = FFTtools::cropWave(gr,rawPeakTime[index]-WINDOW_BEFORE_PEAK,rawPeakTime[index]+WINDOW_AFTER_PEAK);

    pulseDists[index]=distBToHole1 + (rawPeakTime[index]-holeOneTimes[antNum[index]])*TMath::C()/saltRefractiveIndex;
    rxNumber[index]=rxLocation[index];
    peakToPeak[index]=getPeakToPeak(grCrop);
    //    std::cout << peakToPeak[index]*0.5 << "\n";
    distTimesPeakToPeakOverTwo[index]=pulseDists[index]*peakToPeak[index]*0.5;
    distTimesVoltageRatio[index]=distTimesPeakToPeakOverTwo[index]/transmittedVoltage;
  }

  {
    TCanvas *canDist = new TCanvas("canDist","canDist",600,600);
    canDist->Divide(1,2);
    canDist->cd(1);
    TMultiGraph *mgTime = new TMultiGraph();
    TGraph *grSRB = new TGraph(12,&(rxNumber[0]),&(rawPeakTime[0]));
    grSRB->SetMarkerStyle(getMarker(0));
    grSRB->SetMarkerColor(getNiceColour(1));
    mgTime->Add(grSRB,"p");
    TGraph *grSRA = new TGraph(13,&(rxNumber[12]),&(rawPeakTime[12]));
    grSRA->SetMarkerStyle(getMarker(0));
    grSRA->SetMarkerColor(getNiceColour(0));
    mgTime->Add(grSRA,"p");
    TGraph *grLPA = new TGraph(13,&(rxNumber[25]),&(rawPeakTime[25]));
    grLPA->SetMarkerStyle(getMarker(1));
    grLPA->SetMarkerColor(getNiceColour(0));
    mgTime->Add(grLPA,"p");
    TGraph *grLPB = new TGraph(13,&(rxNumber[38]),&(rawPeakTime[38]));
    grLPB->SetMarkerStyle(getMarker(1));
    grLPB->SetMarkerColor(getNiceColour(1));
    mgTime->Add(grLPB,"p");
    mgTime->Draw("ap");
    mgTime->SetTitle("Time of Peak"); 
    mgTime->GetXaxis()->SetTitle("Receiver Position Number");
    mgTime->GetYaxis()->SetTitle("Peak Time (s)");
    TLegend *leggy = new TLegend(0.7,0.2,0.9,0.5);
    leggy->SetFillColor(0);
    leggy->SetFillStyle(0);
    leggy->SetBorderSize(0);
    leggy->AddEntry(grSRA,"SR -- Tx A","p");
    leggy->AddEntry(grSRB,"SR -- Tx B","p");
    leggy->AddEntry(grLPA,"LP -- Tx A","p");
    leggy->AddEntry(grLPB,"LP -- Tx B","p");
    leggy->Draw();


    canDist->cd(2);
    TMultiGraph *mg = new TMultiGraph();
    TGraph *grSRB = new TGraph(12,&(rxNumber[0]),&(pulseDists[0]));
    grSRB->SetMarkerStyle(getMarker(0));
    grSRB->SetMarkerColor(getNiceColour(1));
    mg->Add(grSRB,"p");
    TGraph *grSRA = new TGraph(13,&(rxNumber[12]),&(pulseDists[12]));
    grSRA->SetMarkerStyle(getMarker(0));
    grSRA->SetMarkerColor(getNiceColour(0));
    mg->Add(grSRA,"p");
    TGraph *grLPA = new TGraph(13,&(rxNumber[25]),&(pulseDists[25]));
    grLPA->SetMarkerStyle(getMarker(1));
    grLPA->SetMarkerColor(getNiceColour(0));
    mg->Add(grLPA,"p");
    TGraph *grLPB = new TGraph(13,&(rxNumber[38]),&(pulseDists[38]));
    grLPB->SetMarkerStyle(getMarker(1));
    grLPB->SetMarkerColor(getNiceColour(1));
    mg->Add(grLPB,"p");
    mg->Draw("ap");
    mg->SetTitle("Estimated Path Length");
    mg->GetXaxis()->SetTitle("Receiver Position Number");
    mg->GetYaxis()->SetTitle("Distance (m)");
    canDist->Update();
    canDist->Modified();
  }
  
  {
    TCanvas *canDistTimesV = new TCanvas("canDistTimesV","canDistTimesV",600,400);
    canDistTimesV->SetLogy(1);
    TMultiGraph *mg = new TMultiGraph();
    TGraph *grSRB = new TGraph(12,&(pulseDists[0]),&(distTimesPeakToPeakOverTwo[0]));
    grSRB->SetMarkerStyle(getMarker(0));
    grSRB->SetMarkerColor(getNiceColour(1));
    mg->Add(grSRB,"p");
    TGraph *grSRA = new TGraph(13,&(pulseDists[12]),&(distTimesPeakToPeakOverTwo[12]));
    grSRA->SetMarkerStyle(getMarker(0));
    grSRA->SetMarkerColor(getNiceColour(0));
    mg->Add(grSRA,"p");
    TGraph *grLPA = new TGraph(13,&(pulseDists[25]),&(distTimesPeakToPeakOverTwo[25]));
    grLPA->SetMarkerStyle(getMarker(1));
    grLPA->SetMarkerColor(getNiceColour(0));
    mg->Add(grLPA,"p");
    TGraph *grLPB = new TGraph(13,&(pulseDists[38]),&(distTimesPeakToPeakOverTwo[38]));
    grLPB->SetMarkerStyle(getMarker(1));
    grLPB->SetMarkerColor(getNiceColour(1));
    mg->Add(grLPB,"p");
    mg->Draw("ap");
    mg->SetTitle();
    mg->GetYaxis()->SetTitle("d #times V_{rx} (V m)");
    mg->GetXaxis()->SetTitle("Distance (m)");
    mg->GetYaxis()->SetRangeUser(0.01,10);
    canDistTimesV->Modified();
    canDistTimesV->Update();
    TLegend *leggy = new TLegend(0.1,0.2,0.3,0.45);
    leggy->SetFillColor(0);
    leggy->SetFillStyle(0);
    leggy->SetBorderSize(0);
    leggy->AddEntry(grSRA,"SR -- Tx A","p");
    leggy->AddEntry(grSRB,"SR -- Tx B","p");
    leggy->AddEntry(grLPA,"LP -- Tx A","p");
    leggy->AddEntry(grLPB,"LP -- Tx B","p");
    leggy->Draw();
  }


{
    TCanvas *canFig5 = new TCanvas("canFig5","canFig5",600,400);
    canFig5->SetLogy(1);
    TMultiGraph *mg = new TMultiGraph();
    TGraph *grSRB = new TGraph(12,&(pulseDists[0]),&(distTimesVoltageRatio[0]));
    grSRB->SetMarkerStyle(getMarker(0));
    grSRB->SetMarkerColor(getNiceColour(1));
    mg->Add(grSRB,"p");
    TGraph *grSRA = new TGraph(13,&(pulseDists[12]),&(distTimesVoltageRatio[12]));
    grSRA->SetMarkerStyle(getMarker(0));
    grSRA->SetMarkerColor(getNiceColour(0));
    mg->Add(grSRA,"p");
    TGraph *grLPA = new TGraph(13,&(pulseDists[25]),&(distTimesVoltageRatio[25]));
    grLPA->SetMarkerStyle(getMarker(1));
    grLPA->SetMarkerColor(getNiceColour(0));
    mg->Add(grLPA,"p");
    TGraph *grLPB = new TGraph(13,&(pulseDists[38]),&(distTimesVoltageRatio[38]));
    grLPB->SetMarkerStyle(getMarker(1));
    grLPB->SetMarkerColor(getNiceColour(1));
    mg->Add(grLPB,"p");
    mg->Draw("ap");
    mg->SetTitle();
    mg->GetYaxis()->SetTitle("d #times V_{rx}/V_{tx} (m)");
    mg->GetXaxis()->SetTitle("Distance (m)");
    mg->GetYaxis()->SetRangeUser(1e-4,1);
    TF1 *polly = new TF1("pol0","pol0",0,50);
    polly->SetParameter(0,0.13*TMath::C()/(150e6*saltRefractiveIndex));
    polly->Draw("same");
    canFig5->Modified();
    canFig5->Update();
    TLegend *leggy = new TLegend(0.1,0.2,0.3,0.45);
    leggy->SetFillColor(0);
    leggy->SetFillStyle(0);
    leggy->SetBorderSize(0);
    leggy->AddEntry(grSRA,"SR -- Tx A","p");
    leggy->AddEntry(grSRB,"SR -- Tx B","p");
    leggy->AddEntry(grLPA,"LP -- Tx A","p");
    leggy->AddEntry(grLPB,"LP -- Tx B","p");
    leggy->Draw();
  }


  {
    TCanvas *canDvRxNum = new TCanvas("canDvRxNum","canDvRxNum",600,400);
    canDvRxNum->SetLogy(1);
    TMultiGraph *mg = new TMultiGraph();
    TGraph *grSRB = new TGraph(12,&(rxNumber[0]),&(distTimesPeakToPeakOverTwo[0]));
    grSRB->SetMarkerStyle(getMarker(0));
    grSRB->SetMarkerColor(getNiceColour(1));
    mg->Add(grSRB,"p");
    TGraph *grSRA = new TGraph(13,&(rxNumber[12]),&(distTimesPeakToPeakOverTwo[12]));
    grSRA->SetMarkerStyle(getMarker(0));
    grSRA->SetMarkerColor(getNiceColour(0));
    mg->Add(grSRA,"p");
    TGraph *grLPA = new TGraph(13,&(rxNumber[25]),&(distTimesPeakToPeakOverTwo[25]));
    grLPA->SetMarkerStyle(getMarker(1));
    grLPA->SetMarkerColor(getNiceColour(0));
    mg->Add(grLPA,"p");
    TGraph *grLPB = new TGraph(13,&(rxNumber[38]),&(distTimesPeakToPeakOverTwo[38]));
    grLPB->SetMarkerStyle(getMarker(1));
    grLPB->SetMarkerColor(getNiceColour(1));
    mg->Add(grLPB,"p");
    mg->Draw("ap");
    mg->SetTitle();
    mg->GetYaxis()->SetTitle("d #times V_{rx} (V m)");
    mg->GetXaxis()->SetTitle("Receiver Position Number");
    mg->GetYaxis()->SetRangeUser(0.01,10);
    canDvRxNum->Modified();
    canDvRxNum->Update();
    TLegend *leggy = new TLegend(0.7,0.2,0.9,0.45);
    leggy->SetFillColor(0);
    leggy->SetFillStyle(0);
    leggy->SetBorderSize(0);
    leggy->AddEntry(grSRA,"SR -- Tx A","p");
    leggy->AddEntry(grSRB,"SR -- Tx B","p");
    leggy->AddEntry(grLPA,"LP -- Tx A","p");
    leggy->AddEntry(grLPB,"LP -- Tx B","p");
    leggy->Draw();
  }


}


Int_t getWindowedPeakBin(TGraph *gr, Double_t low, Double_t high) 
{
  Double_t x,y;
  gr->GetPoint(0,x,y);
  Double_t peakVal=y;
  Int_t peakBin=0;
  for(int i=1;i<gr->GetN();i++) {
    gr->GetPoint(i,x,y);
    if(x>=low && x<=high) {
      if(peakVal<y) {
	peakVal=y;
	peakBin=i;
      }      
    }
  }
  return peakBin;
}

Double_t getPeakToPeak(TGraph *gr)
{
  Double_t *xVals=gr->GetX();
  Double_t *yVals=gr->GetY();
  Int_t numPoints=gr->GetN();
  if(numPoints<5) return 0;

  Int_t *posPeakBin = new Int_t[numPoints];
  Int_t *negPeakBin = new Int_t[numPoints];

  Int_t numPos=0;
  Int_t numNeg=0;

  for(int i=1; i<numPoints-1; i++ ) {
    if(yVals[i]>yVals[i-1] && yVals[i]>yVals[i+1]) {
      posPeakBin[numPos]=i;
      numPos++;
    }
    if(yVals[i]<yVals[i-1] && yVals[i]<yVals[i+1]) {
      negPeakBin[numNeg]=i;
      numNeg++;
    }    
  }
  Double_t maxPos=0;
  for(int i=0;i<numPos;i++) {
    if(maxPos<yVals[posPeakBin[i]])
      maxPos=yVals[posPeakBin[i]];
  }
  Double_t maxNeg=0;
  for(int i=0;i<numNeg;i++) {
    if(maxNeg>yVals[negPeakBin[i]])
      maxNeg=yVals[negPeakBin[i]];
  }
  return maxPos-maxNeg;
    


}


Double_t integratePowerSpectrum(TGraph *gr, Double_t low, Double_t high)
{
  Double_t *freqVals=gr->GetX();
  Double_t *powerVals=gr->GetY();
  Int_t numPoints=gr->GetN();
  Double_t retVal=0;
  Double_t df=freqVals[1]-freqVals[0];
  for(int i=0;i<numPoints;i++) {
    //    std::cout << freqVals[i] << "\t" << freqLow[antType] << "\n";
    if(freqVals[i]>=low && freqVals[i]<=high) 
      retVal+=powerVals[i]*df*1e6;
    //    std::cout << freqVals[i] << "\t" << retVal << "\t" << freqLow[antType] << "\n";
  }
  //  std::cout << retVal << "\n";
  return retVal;
}

Double_t attenLengthLogDist(Double_t *x, Double_t *par)
{
  Double_t p0=par[0];
  Double_t L=par[1];
  return p0-x[0]/L;
}

Double_t refractiveIndex(Double_t *x, Double_t *par)
{
  Double_t n=par[0];
  return x[0]*1e9*n/TMath::C();
}

Double_t getLMR600Atten(Double_t freqInMHz, Double_t cableLengthInFt)
{
  Double_t attenPer100ft=0.075550*sqrt(freqInMHz) + 0.000260*freqInMHz;
  Double_t atten=cableLengthInFt*attenPer100ft/100;
  return atten;
}


Double_t getLMR240Atten(Double_t freqInMHz, Double_t cableLengthInFt)
{
  Double_t attenPer100ft=0.242080*sqrt(freqInMHz) + 0.000330*freqInMHz;
  Double_t atten=cableLengthInFt*attenPer100ft/100;
  return atten;
}
