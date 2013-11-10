#include <vector>
#include "Math/InterpolationTypes.h"


TGraph *getWave(char *fileName);

void plotGpr()
{
  gSystem->Load("/sw/lib/libfftw3.so");
  gSystem->Load("libRootFftwWrapper.so");
  gSystem->Load("libMathMore.so");
  //  gStyle->SetHistLineWidth(2);
  //  gStyle->SetLineWidth(2);
  //  TGraph *grTemp=createTemplate(5,14,50000,0,21e-9);
  plotCorrelation(6,23);
  //plotFile(5,21);
}



void plotCorrelation(int disk, int fileNum) {

   gStyle->SetPadTopMargin(0.1);
   gStyle->SetPadBottomMargin(0.15);
   gStyle->SetPadLeftMargin(0.1);
   gStyle->SetPadRightMargin(0.1);
   char filename[180];
   TGraph *grGpr; 
   
   sprintf(filename,"/home/rjn/saltStuff/salt_trip3_data/disk_%d/TEK%05d.CSV",disk,fileNum);
   grGpr=getWave(filename); 
   cout << grGpr->GetN() << endl;
   Double_t *xVals=grGpr->GetX();
   Double_t deltaT=xVals[1]-xVals[0];
   cout << "Delta T Waveform:\t" << deltaT << endl;
   TGraph *grTemp=createTemplate(5,14,grGpr->GetN(),0,21e-9,deltaT);
   TGraph *grCorr=FFTtools::getCorrelationGraph(grGpr,grTemp);
   TGraph *grEnv=FFTtools::getSimplePowerEnvelopeGraph(grGpr);
   TGraph *grCorrEnv=FFTtools::getSimplePowerEnvelopeGraph(grCorr);

   TCanvas *can = new TCanvas("canGpr","canGpr",800,800);
   can->Divide(1,4);
   can->cd(1);
   grGpr->SetLineColor(getNiceColour(1));
   grGpr->SetTitle("Waveform");
   grGpr->Draw("al");
   can->cd(2);
   grEnv->SetLineColor(getNiceColour(1));
   grEnv->SetTitle("Power Envelope");
   grEnv->Draw("al");
   can->cd(3);
   grCorr->SetLineColor(getNiceColour(2));
   grCorr->SetTitle("Correlation With S12 Template");
   grCorr->Draw("al");
   grCorr->GetXaxis()->SetRangeUser(0,1e-5);
   can->cd(4);
   grCorrEnv->SetLineColor(getNiceColour(2));
   grCorrEnv->SetTitle("Correlation Power Envelope");
   grCorrEnv->Draw("al");
   grCorrEnv->GetXaxis()->SetRangeUser(0,1e-5);
   
}

TGraph *createTemplate(int disk, int fileNum, int numPoints, Double_t startTime, Double_t endTime, Double_t deltaT) {
  //for now just assume deltaT is the same
   char filename[180];
   TGraph *grGpr; 
   
   sprintf(filename,"/home/rjn/saltStuff/salt_trip3_data/disk_%d/TEK%05d.CSV",disk,fileNum);
   grGpr=getWave(filename); 
   //   cout << grGpr->GetN() << "\n";
   Double_t *xVals=grGpr->GetX();
   Double_t *yVals=grGpr->GetY();
   Double_t deltaTTemp=xVals[1]-xVals[0];
   cout << "Delta T Template:\t" << deltaTTemp << endl;
   Int_t first=-1,last=-1;
   for(int i=0;i<grGpr->GetN();i++) {
     if(first==-1)
       if(xVals[i]>startTime)
	 first=i;
     
     if(xVals[i]>endTime) {
       last=i;
       break;
     }
   }

  
     //Will use the ROOT::Math::Interpolator function to do this.
   std::vector<double> tVec;
   std::vector<double> vVec;
   for (int samp=0;samp<grGpr->GetN();samp++) {
     tVec.push_back(xVals[samp]);
     vVec.push_back(yVals[samp]);
   }
   ROOT::Math::Interpolator chanInterp(tVec,vVec,ROOT::Math::Interpolation::kAKIM
A);


   //   cout << first << "\t" << last<< endl;
   Double_t *newX = new Double_t [numPoints];
   Double_t *newY = new Double_t [numPoints];

   for(int i=0;i<numPoints;i++) {
     newX[i]=i*deltaT;
     newY[i]=0;
     if(newX[i]<(endTime-startTime)) {
       newY[i]=chanInterp.Eval(newX[i]);
     }
   }
     
   TGraph *grTemp = new TGraph(numPoints,newX,newY);
   delete [] newX;
   delete [] newY;
   delete grGpr;
   //   delete [] xVals;
   //   delete [] yVals;

   //   TCanvas *can = new TCanvas("canGpr","canGpr",600,400);
   //   grTemp->SetLineColor(getNiceColour(1));
   //   grTemp->Draw("al");
   return grTemp;
}

void plotFile(int disk, int fileNum) {
   char filename[180];
   TGraph *grGpr; 
   
   sprintf(filename,"/home/rjn/saltStuff/salt_trip3_data/disk_%d/TEK%05d.CSV",disk,fileNum);
   grGpr=getWave(filename);   
   


   TCanvas *can = new TCanvas("canGpr","canGpr",600,400);
   grGpr->SetLineColor(getNiceColour(1));
   grGpr->Draw("al");
}

TGraph *getWave(char *fileName)
{
  ifstream CsvFile(fileName);
  
  int numLines=0;
  int numLines2=0;
  TString line;
  Double_t xVals[1000000];
  Double_t yVals[1000000];
  Double_t xVals2[1000000];
  Double_t yVals2[1000000];
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
  
TGraph *getScopePeriodiogram(char *fileName)
{
  ifstream CsvFile(fileName);
  
  int numLines=0;
  int numLines2=0;
  TString line;
  Double_t xVals[1000000];
  Double_t yVals[1000000];
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
    numLines++;
    
  }
  Double_t N=numLines*2;
  Double_t deltaF=(xVals[1]-xVals[0]);
  Double_t deltaT=1./(N*deltaF);
  //  cout << deltaT << endl;
  for(int i=0;i<numLines;i++) {
     xVals[i]=i*deltaF;
  }
  //  cout << xVals[0] << "\t" << xVals[numLines-1] << "\t" << xVals[1]-xVals[0] << "\n";

  TGraph *grWave = new TGraph(numLines,xVals,yVals);
  return grWave;
}
  
  
TGraph *getScopeFFTNormToArea(char *fileName)
{
  ifstream CsvFile(fileName);
  
  int numLines=0;
  int numLines2=0;
  TString line;
  Double_t xVals[1000000];
  Double_t yVals[1000000];
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
    numLines++;
    
  }
  Double_t N=numLines*2;
  Double_t deltaF=(xVals[1]-xVals[0])/1e6;
  Double_t deltaT=1e-6./(N*deltaF);
  //  cout << deltaT << endl;
  for(int i=0;i<numLines;i++) {
     xVals[i]=i*deltaF;
     yVals[i]=TMath::Power(10,yVals[i]/10.)*N*deltaT/deltaF;
  }
  //  cout << xVals[0] << "\t" << xVals[numLines-1] << "\t" << xVals[1]-xVals[0] << "\n";

  TGraph *grWave = new TGraph(numLines,xVals,yVals);
  return grWave;
}
