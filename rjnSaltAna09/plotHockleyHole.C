
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
Double_t holeDistances[5]={0,0,44.32,206.30,320.30};
Double_t windowLow[5][4]={{0,0,0,0},{0,0,0,0},{-40.000001e-9,-40.000001e-9,-50.000001e-9,-40.000001e-9},{780.000001e-9,780.000001e-9,790.000001e-9,780.000001e-9},{1120.000001e-9,1120.000001e-9,1150.000001e-9,1120.000001e-9}};
Double_t windowHigh[5][4]={{0,0,0,0},{0,0,0,0},{0.000001e-9,0.000001e-9,50.000001e-9,0.000001e-9},{930.000001e-9,930.000001e-9,890.000001e-9,930.000001e-9},{1270.000001e-9,1270.000001e-9,1250.000001e-9,1270.000001e-9}};

Double_t freqLow[4]={150,150,100,600};
Double_t freqHigh[4]={350,450,150,800};

Double_t noiseOffset[5]={0,0,-100e-9,-1000e-9,-1450e-9};

Char_t *antNames[4]={"Skinny Christian","Short Raven","Long Pole","Baby"};

typedef enum EAntType {
  kSkinnyChristian = 0,
  kShortRaven = 1,
  kLongPole = 2,
  kBaby = 3 
} AntType_t;


#define NUM_FILES 89
Int_t diskNums[NUM_FILES];
Int_t fileNums[NUM_FILES];
Int_t rxHole[NUM_FILES];
Int_t antNum[NUM_FILES];
Int_t ampInDb[NUM_FILES];


TGraph *getWave(char *fileName);
void plotWaveAndFFT(char *fileName, char *graphTitle=0);
void readInputFile();
Double_t getTriggerTime(Int_t holeInd);
Double_t getRFCableDelay();
void plotAllPeakTimes();
void plotThisOne(Int_t diskNum, Int_t fileNum);
void plotThisOne(Int_t index);
char *getAntName(EAntType antType);
Int_t getWindowedPeakBin(TGraph *gr, Double_t low, Double_t high);
Double_t integratePowerSpectrum(TGraph *gr, Double_t low, Double_t high);
Double_t attenLengthLogDist(Double_t *x, Double_t *par);
Double_t refractiveIndex(Double_t *x, Double_t *par);

void plotThisOne(Int_t index) 
{

  char graphTitle[180];
  char fileName[180];
  sprintf(fileName,"/home/rjn/saltStuff/hockley2009/disk_%d/TEK%05d.txt",diskNums[index],fileNums[index]);
  sprintf(graphTitle,"Ant %d, Hole %d",antNum[index],rxHole[index]);
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


void plotHockleyHole()
{
  gSystem->Load("/sw/lib/libfftw3.so");
  gSystem->Load("libRootFftwWrapper.so");
  readInputFile();
  plotAllPeakTimes();
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
  //  cout << timeVals[peakBin] << "\n";
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

  ifstream TxtFile("hockleyFileLog.csv");  
  int numLines=0;
  int numLines2=0;
  TString line;
  line.ReadLine(TxtFile);
  while(line.ReadLine(TxtFile)) {
    TObjArray *vals = (TObjArray*) line.Tokenize(",");
    //    vals.Dump();
    TObjString *tempDisk=vals->At(0);
    TObjString *tempFile=vals->At(1);
    TObjString *tempHole=vals->At(2);
    TObjString *tempAnt=vals->At(3);
    TObjString *tempAmp=vals->At(4);

    diskNums[numLines]=tempDisk->GetString().Atoi();
    fileNums[numLines]=tempFile->GetString().Atoi();
    rxHole[numLines]=tempHole->GetString().Atoi();
    ampInDb[numLines]=tempAmp->GetString().Atoi();
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
//     std::cout << diskNums[numLines] << "\t" 
// 	      << fileNums[numLines] << "\t"
// 	      << rxHole[numLines] << "\t" 
// 	      << antNum[numLines] << "\n";
//    std::cout << diskNums[numLines] << "\t" 
//	      << fileNums[numLines] << "\t"
//	      << rxHole[numLines] << "\t" 
//	      << antNum[numLines] << "\t" 
//	      << ampInDb[numLines] << "\n";
    numLines++;
  }
  //  std::cout << numLines << endl;
}


Double_t getTriggerTime(Int_t holeInd)
{
  Double_t lmr200Length[5]={0,0,400,800,1200};
  
  return syncLag+lmr200Length[holeInd]*lmr200NsPerFt;


}

Double_t getRFCableDelay()
{
  Double_t rg58Delay=3*rg58NsPerFt;
  Double_t lmr240Delay=10*lmr240NsPerFt;
  Double_t lmr600Delay=100*lmr600NsPerFt;
  return rg58Delay+cps3InsertionDelayNs+lmr240Delay+lmr600Delay;
}

void plotAllPeakTimes()
{
  Double_t holeDists[NUM_FILES];
  Double_t rawPeak[NUM_FILES];
  Double_t correctedPeak[NUM_FILES];
  Double_t integratedPower[NUM_FILES];
  Double_t sqrtPower[NUM_FILES];
  Double_t sqrtPowerTimesDist[NUM_FILES];
  Double_t logSqrtPowerTimesDist[NUM_FILES];
  Double_t intPowerFreq[NUM_FILES];
  Double_t logSqrtIntPowerFreqTimesDist[NUM_FILES];

  Int_t antCount[4]={0};
  Int_t antFreqCount[4][8]={0};
  Double_t antPeakToPeak[4][NUM_FILES];
  Double_t antLogHalfPeakToPeakTimesDist[4][NUM_FILES];
  Double_t antLogSqrtPowerTimesDist[4][NUM_FILES];
  Double_t antLogSqrtIntPowerFreqTimesDist[4][NUM_FILES];
  Double_t antLogSqrtIntPowerFreqTimesDistFreqBins[4][8][NUM_FILES];
  Double_t antHoleDist[4][NUM_FILES];

  for(int index=0;index<NUM_FILES;index++) {    
    char fileName[180];
    sprintf(fileName,"/home/rjn/saltStuff/hockley2009/disk_%d/TEK%05d.txt",diskNums[index],fileNums[index]);
    TGraph *gr=getWave(fileName);
    TGraph *grHilbert = FFTtools::getHilbertEnvelope(gr);
    Int_t peakBin=getWindowedPeakBin(grHilbert,windowLow[rxHole[index]][antNum[index]],windowHigh[rxHole[index]][antNum[index]]);
    Double_t *timeVals=grHilbert->GetX();
    rawPeak[index]=timeVals[peakBin];
    correctedPeak[index]=(timeVals[peakBin]*1e9)+getTriggerTime(rxHole[index])-getRFCableDelay();
    holeDists[index]=holeDistances[rxHole[index]];


    TGraph *grCrop = FFTtools::cropWave(gr,windowLow[rxHole[index]][antNum[index]],windowHigh[rxHole[index]][antNum[index]]);
    TGraph *grNoiseCrop = FFTtools::cropWave(gr,windowLow[rxHole[index]][antNum[index]]+noiseOffset[rxHole[index]],windowHigh[rxHole[index]][antNum[index]]+noiseOffset[rxHole[index]]);

    integratedPower[index]=FFTtools::integrateVoltageSquared(grCrop);
    integratedPower[index]/=TMath::Power(10,ampInDb[index]/10.);
    sqrtPower[index]=TMath::Sqrt(integratedPower[index]);
    sqrtPowerTimesDist[index]=sqrtPower[index]*holeDists[index];
    logSqrtPowerTimesDist[index]=TMath::Log(sqrtPowerTimesDist[index]);
//     if(antNum[index]==3) {
//       cout << index << "\t" << diskNums[index] << "\t" << fileNums[index] << "\t" << rxHole[index] << "\t" << correctedPeak[index] << "\t" << ampInDb[index] << "\t" << timeVals[peakBin] << "\n";
//       TCanvas *can1 = new TCanvas();
//       gr->Draw("al");
//       TCanvas *can = new TCanvas();
//       can->Divide(2,2);
//       can->cd(1);
//       grCrop->Draw("al");
//       can->cd(2);
//       TGraph *grCropPowerDb = FFTtools::makePowerSpectrumVoltsSecondsdB(grCrop);
//       grCropPowerDb->Draw("al");
//       can->cd(3);
//       grNoiseCrop->Draw("al");
//       can->cd(4);
//       TGraph *grNoiseCropPowerDb = FFTtools::makePowerSpectrumVoltsSecondsdB(grNoiseCrop);
//       grNoiseCropPowerDb->Draw("al");
//     }


    //Now integrate the power in the frequency domain to use in the attenuation length calculation.
    TGraph *grCropPower = FFTtools::makePowerSpectrumVoltsSeconds(grCrop);
    TGraph *grNoiseCropPower = FFTtools::makePowerSpectrumVoltsSeconds(grNoiseCrop);
    if(grCrop->GetN()!=grNoiseCrop->GetN()) {
      std::cout << "Size mismatch\t" << grCrop->GetN() << "\t" << grNoiseCrop->GetN() << "\n";
      
    }
    intPowerFreq[index]=integratePowerSpectrum(grCropPower,freqLow[antNum[index]],freqHigh[antNum[index]]);
    intPowerFreq[index]/=TMath::Power(10,ampInDb[index]/10.);
    
    intPowerFreq[index]-=(integratePowerSpectrum(grNoiseCropPower,freqLow[antNum[index]],freqHigh[antNum[index]]))/TMath::Power(10,ampInDb[index]/10.);
    if(intPowerFreq[index]<0) {
      //No signal above noise.
      continue;
    }
    //    std::cout << intPowerFreq[index] << "\t";
    logSqrtIntPowerFreqTimesDist[index]=TMath::Log(holeDists[index]*TMath::Sqrt(intPowerFreq[index]));
    //    std::cout << intPowerFreq[index] << "\t" << logSqrtIntPowerFreqTimesDist[index] << "\n";

    //    if(rxHole[index]==4)
    //      continue;

        if(antNum[index]==2 && rxHole[index]==3)  {
              continue;
	  //Saturated long pole
        }
    antHoleDist[antNum[index]][antCount[antNum[index]]]=holeDists[index];
    antLogSqrtPowerTimesDist[antNum[index]][antCount[antNum[index]]]=logSqrtPowerTimesDist[index];
    antLogSqrtIntPowerFreqTimesDist[antNum[index]][antCount[antNum[index]]]=logSqrtIntPowerFreqTimesDist[index];
    antPeakToPeak[antNum[index]][antCount[antNum[index]]]=getPeakToPeak(grCrop);
    antPeakToPeak[antNum[index]][antCount[antNum[index]]]/=TMath::Power(20,ampInDb[index]/10.);
    antLogHalfPeakToPeakTimesDist[antNum[index]][antCount[antNum[index]]]=TMath::Log(0.5*holeDists[index]*antPeakToPeak[antNum[index]][antCount[antNum[index]]]);


    antCount[antNum[index]]++;    
    //Now play with frequency bins
    for(int freqBin=0;freqBin<8;freqBin++) {
      Double_t minFreq=100*freqBin+50;
      Double_t maxFreq=100*freqBin+150;
      Double_t signalPower=integratePowerSpectrum(grCropPower,minFreq,maxFreq);
      Double_t noisePower=integratePowerSpectrum(grNoiseCropPower,minFreq,maxFreq);
      signalPower-=noisePower;
      if(signalPower<0) continue;
      signalPower/=TMath::Power(10,ampInDb[index]/10.);
      antLogSqrtIntPowerFreqTimesDistFreqBins[antNum[index]][freqBin][antFreqCount[antNum[index]][freqBin]]=TMath::Log(TMath::Sqrt(signalPower));
      antFreqCount[antNum[index]][freqBin]++;
    }								     
  }
  TCanvas *canDist = new TCanvas("canDist","canDist",600,400);
  //  canDist->Divide(1,2);
  //  canDist->cd(1);
  //  TGraph *grRaw = new TGraph(NUM_FILES,holeDists,rawPeak);
  //  grRaw->Draw("ap");
  //  canDist->cd(2);
  TF1 *fitRef = new TF1("fitRef",refractiveIndex,0,400,1);
  fitRef->SetParName(0,"n_{salt}");
  TGraph *grCor = new TGraph(NUM_FILES,holeDists,correctedPeak);
  grCor->SetTitle("Time Of Flight -- Refractive Index");
  grCor->Draw("ap");
  grCor->GetXaxis()->SetTitle("Hole Distance (m)");
  grCor->GetYaxis()->SetTitle("Peak Time (ns)");
  grCor->Fit("fitRef");
  canDist->Update();
  canDist->Modified();
  sortOutTitle();

  TF1 *fitty = new TF1("fitty",attenLengthLogDist,0,400,2);
  fitty->SetParameters(-5,40);
  fitty->SetParNames("P_0","L_{att}");

  gStyle->SetStatH(0.3);
  gStyle->SetStatW(0.2);
  gStyle->SetTitleSize(0.08);
  gStyle->SetTitleOffset(0.6,"y");
  gStyle->SetTitleOffset(0.8,"x");
  TCanvas *can = new TCanvas("canAnt","canAnt",800,800);
  //  TMultiGraph *mg = new TMultiGraph();
  can->Divide(1,4);
  for(int ant=0;ant<4;ant++) {						
    can->cd(ant+1);
    //    TGraph *grSlope = new TGraph(antCount[ant],antHoleDist[ant],antLogSqrtPowerTimesDist[ant]);
    //    TGraph *grSlope = new TGraph(antCount[ant],antHoleDist[ant],antLogSqrtIntPowerFreqTimesDist[ant]);
    TGraph *grSlope = new TGraph(antCount[ant],antHoleDist[ant],antLogHalfPeakToPeakTimesDist[ant]);
    grSlope->SetMarkerStyle(getMarker(ant));
    grSlope->SetMarkerColor(getNiceColour(ant));
    grSlope->SetTitle(antNames[ant]);
    //    mg->Add(grSlope,"p");
    grSlope->Draw("ap");
    grSlope->Fit("fitty");
    grSlope->GetXaxis()->SetTitle("Hole Distance (m)");
    grSlope->GetYaxis()->SetTitle("Log(d #sqrt{P})");
    can->Update();
    can->Modified();
    sortOutTitle(0.09);
  }


//   TCanvas *canRaven = new TCanvas("canRaven","canRaven");
//   //  TMultiGraph *mg = new TMultiGraph();
//   canRaven->Divide(1,8);
//   for(int ant=1;ant<2;ant++) {						
//     for(int freqBin=0;freqBin<8;freqBin++) {
//       canRaven->cd(freqBin+1);
//       TGraph *grSlope = new TGraph(antFreqCount[ant][freqBin],antHoleDist[ant],antLogSqrtIntPowerFreqTimesDistFreqBins[ant][freqBin]);
//       grSlope->SetMarkerStyle(getMarker(ant));
//       grSlope->SetMarkerColor(getNiceColour(freqBin));
//       grSlope->Draw("ap");
//       grSlope->Fit("fitty");
//     }
//   }
  
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
