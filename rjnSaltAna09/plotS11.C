#define NUM_OPENS 7
int openFiles[NUM_OPENS][2]={{1,5},{1,12},{1,21},{1,48},{1,49},{2,4},{5,0}};
int cropOpen[NUM_OPENS]={0,0,1,1,1,1,1};
Double_t minTimeOpen[NUM_OPENS]={0,0,0.34e-6,0.34e-6,0.34e-6,0.12e-6,0.12e-6};
Double_t maxTimeOpen[NUM_OPENS]={0,0,0.39e-6,0.39e-6,0.39e-6,0.19e-6,0.19e-6};


TGraph *getWaveFromCsv(char *fileName);
TGraph *getWaveFromTxt(char *fileName);
TGraph *getOpenWave(int whichOpen);
TGraph *cropLikeOpen(TGraph *grIn, int whichOpen);

void plotS11()
{
  ///  gROOT->SetStyle("BABAR");
  gSystem->Load("libgsl.so");
  gSystem->Load("libfftw3.so");
  gSystem->Load("libMathMore.so");
  gSystem->Load("libRootFftwWrapper.so");
  //  Gstyle->SetHistLineWidth(2);
  //  gStyle->SetLineWidth(2);
  //    plotS11Christian();
  //  plotS11Raven();
  //  plotS11SkinnyChristian();
  //  plotS11LongPole();
  //  plotS11Baby();
  plotOpens();
}

TGraph *getOpenWave(int whichOpen) {
   char filename[180];
  TGraph *grOpen;
  sprintf(filename,"/Users/rjn/saltStuff/hockley2009/disk_%d/TEK%05d.txt",openFiles[whichOpen][0],openFiles[whichOpen][1]);
  grOpen=getWaveFromTxt(filename); 
  if(cropOpen[whichOpen]) {
    TGraph *grTemp = grOpen;
    grOpen=FFTtools::cropWave(grTemp,minTimeOpen[whichOpen],maxTimeOpen[whichOpen]);
    delete grTemp;
  }
  return grOpen;  
}

TGraph *cropLikeOpen(TGraph *grIn, int whichOpen)
{
  if(cropOpen[whichOpen]) {
    TGraph *grTemp = FFTtools::cropWave(grIn,minTimeOpen[whichOpen],maxTimeOpen[whichOpen]);
    return grTemp;
  }
  return grIn;
}


void plotOpens() {
   TGraph *grOpen[NUM_OPENS]; 
   TGraph *grOpenPower[NUM_OPENS];
   TGraph *grOpenPowerBart[NUM_OPENS];
   TGraph *grOpenPowerBartSeg[NUM_OPENS];
   for(int i=0;i<NUM_OPENS;i++) {
     grOpen[i]=getOpenWave(i);
      cout << grOpen[i]->GetN() << endl;
      grOpenPower[i]=FFTtools::makePowerSpectrumVoltsSeconds(grOpen[i]);
      grOpenPowerBart[i]=FFTtools::makePowerSpectrumVoltsSecondsBartlett(grOpen[i]);
      grOpenPowerBartSeg[i]=FFTtools::makePSVSBartlettPaddedOverlap(grOpen[i],32,1024);
   }
   TCanvas *canOpen  = new TCanvas("canOpen","canOpen",600,600);
   canOpen->Divide(1,2);
   canOpen->cd(1);
   grOpen[0]->SetTitle("Open -- Inside Office (x2)");
   //   grOpen[1]->SetTitle("Open -- Below Ground");
   for(int i=0;i<2;i++) {
      
      grOpen[i]->SetLineColor(getNiceColour(i));
      grOpen[i]->SetLineWidth(2);
      //      grOpen[i]->SetTitle("Open");
      if(i==0) {
	 grOpen[i]->Draw("al");
	 grOpen[i]->GetXaxis()->SetTitle("Time (s)");
	 grOpen[i]->GetYaxis()->SetTitle("Voltage (V)");
      }
      else {
	 grOpen[i]->Draw("l");
      }
   }
   grOpenPower[0]->SetTitle("Open -- Inside Office (x2)");
   //   grOpenPower[0]->SetTitle("Open -- Above Ground");
   //   grOpenPower[1]->SetTitle("Open -- Below Ground");
   canOpen->cd(2);
   for(int i=0;i<2;i++) {
      gPad->SetLogy();
      grOpenPower[i]->SetLineColor(getNiceColour(i));
      grOpenPower[i]->SetLineWidth(2);
      //      grOpenPower[i]->SetTitle("Open");
      if(i==0) {
	 grOpenPower[i]->Draw("al");
	 grOpenPower[i]->SetMaximum(1e-12);
	 grOpenPower[i]->SetMinimum(1e-16);
	 grOpenPower[i]->GetXaxis()->SetRangeUser(0,1000);
	 grOpenPower[i]->GetXaxis()->SetTitle("Frequency (MHz)");
	 grOpenPower[i]->GetYaxis()->SetTitle("Power (normalised to area)");
      }
      else
	 grOpenPower[i]->Draw("l");
   }




//    canOpen->cd(3);
//    for(int i=0;i<2;i++) {
//       gPad->SetLogy();
//       grOpenPowerBart[i]->SetLineColor(getNiceColour(i));
//       grOpenPowerBart[i]->SetLineWidth(2);
//       //      grOpenPowerBart[i]->SetTitle("Open");
//       if(i==0) {
// 	 grOpenPowerBart[i]->Draw("al");
// 	 grOpenPowerBart[i]->SetMaximum(1e-12);
// 	 grOpenPowerBart[i]->SetMinimum(1e-16);
// 	 grOpenPowerBart[i]->GetXaxis()->SetRangeUser(0,1000);
// 	 grOpenPowerBart[i]->GetXaxis()->SetTitle("Frequency (MHz)");
// 	 grOpenPowerBart[i]->GetYaxis()->SetTitle("Power (normalised to area)");
//       }
//       else
// 	 grOpenPowerBart[i]->Draw("l");
//    }

//    canOpen->cd(4);
//    for(int i=0;i<2;i++) {
//       gPad->SetLogy();
//       grOpenPowerBartSeg[i]->SetLineColor(getNiceColour(i));
//       grOpenPowerBartSeg[i]->SetLineWidth(2);
//       //      grOpenPowerBartSeg[i]->SetTitle("Open");
//       if(i==0) {
// 	 grOpenPowerBartSeg[i]->Draw("al");
// 	 grOpenPowerBartSeg[i]->SetMaximum(1e-12);
// 	 grOpenPowerBartSeg[i]->SetMinimum(1e-16);
// 	 grOpenPowerBartSeg[i]->GetXaxis()->SetRangeUser(0,1000);
// 	 grOpenPowerBartSeg[i]->GetXaxis()->SetTitle("Frequency (MHz)");
// 	 grOpenPowerBartSeg[i]->GetYaxis()->SetTitle("Power (normalised to area)");
//       }
//       else
// 	 grOpenPowerBartSeg[i]->Draw("l");
//    }


   TCanvas *canOpen2  = new TCanvas("canOpen2","canOpen2",600,600);
   canOpen2->Divide(1,2);
   canOpen2->cd(1);
   grOpen[2]->SetTitle("Open -- Outside (x2)");
   //   grOpen[1]->SetTitle("Open -- Below Ground");
   for(int i=0;i<2;i++) {
     
      grOpen[i+2]->SetLineColor(getNiceColour(i));
      grOpen[i+2]->SetLineWidth(2);
      //      grOpen[i]->SetTitle("Open");
      if(i==0) {
	 grOpen[i+2]->Draw("al");
	 grOpen[i+2]->GetXaxis()->SetTitle("Time (s)");
	 grOpen[i+2]->GetYaxis()->SetTitle("Voltage (V)");
      }
      else {
	 grOpen[i+2]->Draw("l");
      }
   }
   grOpenPower[2]->SetTitle("Open -- Outside (x2)");
   //   grOpenPower[0]->SetTitle("Open -- Above Ground");
   //   grOpenPower[1]->SetTitle("Open -- Below Ground");
   canOpen2->cd(2);
   for(int i=0;i<2;i++) {
      gPad->SetLogy();
      grOpenPower[i+2]->SetLineColor(getNiceColour(i));
      grOpenPower[i+2]->SetLineWidth(2);
      //      grOpenPower[i+2]->SetTitle("Open");
      if(i==0) {
	 grOpenPower[i+2]->Draw("al");
	 grOpenPower[i+2]->SetMaximum(1e-12);
	 grOpenPower[i+2]->SetMinimum(1e-16);
	 grOpenPower[i+2]->GetXaxis()->SetRangeUser(0,1000);
	 grOpenPower[i+2]->GetXaxis()->SetTitle("Frequency (MHz)");
	 grOpenPower[i+2]->GetYaxis()->SetTitle("Power (normalised to area)");
      }
      else
	 grOpenPower[i+2]->Draw("l");
   }

//    canOpen2->cd(3);
//    for(int i=0;i<2;i++) {
//       gPad->SetLogy();
//       grOpenPowerBart[i]->SetLineColor(getNiceColour(i));
//       grOpenPowerBart[i]->SetLineWidth(2);
//       //      grOpenPowerBart[i]->SetTitle("Open");
//       if(i==0) {
// 	 grOpenPowerBart[i]->Draw("al");
// 	 grOpenPowerBart[i]->SetMaximum(1e-12);
// 	 grOpenPowerBart[i]->SetMinimum(1e-16);
// 	 grOpenPowerBart[i]->GetXaxis()->SetRangeUser(0,1000);
// 	 grOpenPowerBart[i]->GetXaxis()->SetTitle("Frequency (MHz)");
// 	 grOpenPowerBart[i]->GetYaxis()->SetTitle("Power (normalised to area)");
//       }
//       else
// 	 grOpenPowerBart[i]->Draw("l");
//    }

//    canOpen2->cd(4);
//    for(int i=0;i<2;i++) {
//       gPad->SetLogy();
//       grOpenPowerBartSeg[i]->SetLineColor(getNiceColour(i));
//       grOpenPowerBartSeg[i]->SetLineWidth(2);
//       //      grOpenPowerBartSeg[i]->SetTitle("Open");
//       if(i==0) {
// 	 grOpenPowerBartSeg[i]->Draw("al");
// 	 grOpenPowerBartSeg[i]->SetMaximum(1e-12);
// 	 grOpenPowerBartSeg[i]->SetMinimum(1e-16);
// 	 grOpenPowerBartSeg[i]->GetXaxis()->SetRangeUser(0,1000);
// 	 grOpenPowerBartSeg[i]->GetXaxis()->SetTitle("Frequency (MHz)");
// 	 grOpenPowerBartSeg[i]->GetYaxis()->SetTitle("Power (normalised to area)");
//       }
//       else
// 	 grOpenPowerBartSeg[i]->Draw("l");
//    }


   TCanvas *canOpen3  = new TCanvas("canOpen3","canOpen3",600,600);
   canOpen3->Divide(1,2);
   canOpen3->cd(1);
   grOpen[5]->SetTitle("Open -- Salt Mine (Day's I & III)");
   //   grOpen[1]->SetTitle("Open -- Below Ground");
   for(int i=0;i<2;i++) {
     
      grOpen[i+5]->SetLineColor(getNiceColour(i));
      grOpen[i+5]->SetLineWidth(2);
      //      grOpen[i]->SetTitle("Open");
      if(i==0) {
	 grOpen[i+5]->Draw("al");
	 grOpen[i+5]->GetXaxis()->SetTitle("Time (s)");
	 grOpen[i+5]->GetYaxis()->SetTitle("Voltage (V)");
      }
      else {
	 grOpen[i+5]->Draw("l");
      }
   }
   grOpenPower[5]->SetTitle("Open -- Salt Mine (Day's I & III)");
   //   grOpenPower[0]->SetTitle("Open -- Above Ground");
   //   grOpenPower[1]->SetTitle("Open -- Below Ground");
   canOpen3->cd(2);
   for(int i=0;i<2;i++) {
      gPad->SetLogy();
      grOpenPower[i+5]->SetLineColor(getNiceColour(i));
      grOpenPower[i+5]->SetLineWidth(2);
      //      grOpenPower[i+5]->SetTitle("Open");
      if(i==0) {
	 grOpenPower[i+5]->Draw("al");
	 grOpenPower[i+5]->SetMaximum(1e-11);
	 grOpenPower[i+5]->SetMinimum(1e-16);
	 grOpenPower[i+5]->GetXaxis()->SetRangeUser(0,1000);
	 grOpenPower[i+5]->GetXaxis()->SetTitle("Frequency (MHz)");
	 grOpenPower[i+5]->GetYaxis()->SetTitle("Power (normalised to area)");
      }
      else
	 grOpenPower[i+5]->Draw("l");
   }

}


void plotS11Raven() {

   const int numFiles=9;
   char filename[180];
   int antFiles[numFiles][2]={{1,13},{1,14},{1,26},{1,25},{1,31},{2,6},{2,9},{5,4},{5,5}};
   int whichOpen[numFiles]={1,1,2,2,3,5,5,6,6};
   int plotThis[numFiles]={1,1,1,1,1,2,2,2,2};
   char *graphDesc[numFiles]={"Office (1)","Office (2)","Outside (1)","Outside (2)","Outside (2)","Salt Day I (I)","Salt Day I (2)","Salt Day III (1)","Salt Day III (2)"};
   TGraph *grOpen[NUM_OPENS]; 
   TGraph *grOpenPower[NUM_OPENS];
   for(int i=0;i<NUM_OPENS;i++) {
      grOpen[i]=getOpenWave(i); 
      //      grOpenPower[i]=FFTtools::makePowerSpectrumVoltsSecondsPadded(grOpen[i],8);
      grOpenPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grOpen[i],32,1024);
      Double_t *freq=grOpenPower[i]->GetX();
      //      cout << i << "\t" << grOpenPower[i]->GetN() << "\t" <<  freq[1]-freq[0] << "\n";
   }
   
   TGraph *grAnt[numFiles];
   TGraph *grAntPower[numFiles];
   TGraph *grAntRatio[numFiles];
   TGraph *grAntTrans[numFiles];
   for(int i=0;i<numFiles;i++) {
       sprintf(filename,"/home/rjn/saltStuff/hockley2009/disk_%d/TEK%05d.txt",antFiles[i][0],antFiles[i][1]);
       grAnt[i]=getWaveFromTxt(filename);
       grAnt[i]=cropLikeOpen(grAnt[i],whichOpen[i]);
       if(i<3) {
	 //	  grAntPower[i]=FFTtools::makePowerSpectrumVoltsSecondsPadded(grAnt[i],8);
	  grAntPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grAnt[i],32,1024);
       }
       else {
	  grAntPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grAnt[i],32,1024);
	  //	  TGraph *grTemp=FFTtools::makePowerSpectrumVoltsSecondsPadded(grAnt[i],8);
	  //	  grAntPower[i]=FFTtools::smoothFFT(grTemp,2);
       }
       grAntRatio[i]=FFTtools::dbGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
       grAntTrans[i]=FFTtools::ratioSubtractOneGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
       //       Double_t *freq=grAntPower[i]->GetX();
       //       cout << i << "\t" << freq[1]-freq[0] << "\n";
       //       cout << i << "\t" << grAntPower[i]->GetN() << "\t" <<  freq[1]-freq[0] << "\n";
       //       cout << grAntRatio[i] << "\t" << grAntTrans[i] << endl;
   }


  TCanvas *canWave = new TCanvas("canRavenWave","canRavenWave",800,800);
  canWave->Divide(1,numFiles);
  for(int i=0;i<numFiles;i++) {
     canWave->cd(i+1);
     grAnt[i]->SetLineColor(getNiceColour(i));     
     grAnt[i]->Draw("al");
     grAnt[i]->SetTitle(graphDesc[i]);
     grAnt[i]->GetXaxis()->SetTitle("Time (s)");
     grAnt[i]->GetYaxis()->SetTitle("Voltage (V)");

  }

  TCanvas *can = new TCanvas("canRaven","canRaven",900,600);
  can->Divide(1,2);
  for(int subPad=1;subPad<=2;subPad++) {
    can->cd(subPad);
    TLegend *leggy = new TLegend(0.6,0.2,0.8,0.5);
    leggy->SetBorderSize(0);
    leggy->SetFillColor(0);
    leggy->SetFillStyle(0);  
    TMultiGraph *mg = new TMultiGraph();
    for(int i=0;i<numFiles;i++) {
      if(plotThis[i]!=subPad) continue;
      grAntRatio[i]->SetLineColor(getNiceColour(i));
      grAntRatio[i]->SetLineWidth(2);
      mg->Add(grAntRatio[i],"l");
      leggy->AddEntry(grAntRatio[i],graphDesc[i],"l");
    }
    mg->Draw("al");
    mg->SetTitle("S11 Measurements Short Raven");
    mg->SetMaximum(0);
    mg->SetMinimum(-30);
    mg->GetXaxis()->SetRangeUser(0,1000);
    mg->GetXaxis()->SetTitle("Frequency (MHz)");
    mg->GetYaxis()->SetTitle("Reflected Power (dB)");
    leggy->Draw();
  }


  TCanvas *canTrans = new TCanvas("canTransRaven","canTransRaven",900,600);
  canTrans->Divide(1,2);
  for(int subPad=1;subPad<=2;subPad++) {
    canTrans->cd(subPad);
    TLegend *leggy2 = new TLegend(0.3,0.2,0.5,0.5);
    leggy2->SetBorderSize(0);
    leggy2->SetFillColor(0);
    leggy2->SetFillStyle(0);  
    TMultiGraph *mg =new TMultiGraph();
    for(int i=0;i<numFiles;i++) {
      if(plotThis[i]!=subPad) continue;
      grAntTrans[i]->SetLineColor(getNiceColour(i));
      grAntTrans[i]->SetLineWidth(2);
      mg->Add(grAntTrans[i],"l");
      leggy2->AddEntry(grAntTrans[i],graphDesc[i],"l");
    }
    mg->Draw("al");
    mg->SetTitle("S11 Measurements Short Raven");
    mg->SetMaximum(1);
    mg->SetMinimum(0);
    mg->GetXaxis()->SetRangeUser(0,1000);
    mg->GetXaxis()->SetTitle("Frequency (MHz)");
    mg->GetYaxis()->SetTitle("Transmission (Ratio)");
    leggy2->Draw();
  }

}

void plotS11FatChristian() {

   const int numFiles=7;
   char filename[180];
   int antFiles[numFiles][2]={{1,3},{1,7},{1,10},{1,9},{1,19},{1,20},{1,22}};
   int whichOpen[numFiles]={0,0,1,1,1,2,2};
   int plotThis[numFiles]={1,1,1,1,1,1,1};
   char *graphDesc[numFiles]={"Office (1)","Office (2)","Office (1)","Office (2)","Outside (1)","Outside (1)","Outside (2)"};
   TGraph *grOpen[NUM_OPENS]; 
   TGraph *grOpenPower[NUM_OPENS];
   for(int i=0;i<NUM_OPENS;i++) {
      grOpen[i]=getOpenWave(i); 
      //      grOpenPower[i]=FFTtools::makePowerSpectrumVoltsSecondsPadded(grOpen[i],8);
      grOpenPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grOpen[i],32,1024);
      Double_t *freq=grOpenPower[i]->GetX();
      //      cout << i << "\t" << grOpenPower[i]->GetN() << "\t" <<  freq[1]-freq[0] << "\n";
   }
   
   TGraph *grAnt[numFiles];
   TGraph *grAntPower[numFiles];
   TGraph *grAntRatio[numFiles];
   TGraph *grAntTrans[numFiles];
   for(int i=0;i<numFiles;i++) {
       sprintf(filename,"/home/rjn/saltStuff/hockley2009/disk_%d/TEK%05d.txt",antFiles[i][0],antFiles[i][1]);
       grAnt[i]=getWaveFromTxt(filename);
       grAnt[i]=cropLikeOpen(grAnt[i],whichOpen[i]);
       if(i<3) {
	 //	  grAntPower[i]=FFTtools::makePowerSpectrumVoltsSecondsPadded(grAnt[i],8);
	  grAntPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grAnt[i],32,1024);
       }
       else {
	  grAntPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grAnt[i],32,1024);
	  //	  TGraph *grTemp=FFTtools::makePowerSpectrumVoltsSecondsPadded(grAnt[i],8);
	  //	  grAntPower[i]=FFTtools::smoothFFT(grTemp,2);
       }
       grAntRatio[i]=FFTtools::dbGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
       grAntTrans[i]=FFTtools::ratioSubtractOneGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
       //       Double_t *freq=grAntPower[i]->GetX();
       //       cout << i << "\t" << freq[1]-freq[0] << "\n";
       //       cout << i << "\t" << grAntPower[i]->GetN() << "\t" <<  freq[1]-freq[0] << "\n";
       //       cout << grAntRatio[i] << "\t" << grAntTrans[i] << endl;
   }


  TCanvas *canWave = new TCanvas("canFatChristianWave","canFatChristianWave",800,800);
  canWave->Divide(1,numFiles);
  for(int i=0;i<numFiles;i++) {
     canWave->cd(i+1);
     grAnt[i]->SetLineColor(getNiceColour(i));     
     grAnt[i]->Draw("al");
     grAnt[i]->SetTitle(graphDesc[i]);
     grAnt[i]->GetXaxis()->SetTitle("Time (s)");
     grAnt[i]->GetYaxis()->SetTitle("Voltage (V)");

  }

  TCanvas *can = new TCanvas("canFatChristian","canFatChristian",900,600);
  TLegend *leggy = new TLegend(0.6,0.2,0.8,0.5);
  leggy->SetBorderSize(0);
  leggy->SetFillColor(0);
  leggy->SetFillStyle(0);  
  for(int i=0;i<numFiles;i++) {
     if(!plotThis[i]) continue;
     grAntRatio[i]->SetLineColor(getNiceColour(i));
     grAntRatio[i]->SetLineWidth(2);
     if(i==0) {
	grAntRatio[i]->Draw("al");
	grAntRatio[i]->SetTitle("S11 Measurements Fat Christian");
	grAntRatio[i]->SetMaximum(0);
	grAntRatio[i]->SetMinimum(-20);
	grAntRatio[i]->GetXaxis()->SetRangeUser(0,1000);
	grAntRatio[i]->GetXaxis()->SetTitle("Frequency (MHz)");
	grAntRatio[i]->GetYaxis()->SetTitle("Reflected Power (dB)");
     }
     grAntRatio[i]->Draw("l");
     leggy->AddEntry(grAntRatio[i],graphDesc[i],"l");
  }
  leggy->Draw();


  TCanvas *canTrans = new TCanvas("canTransFatChristian","canTransFatChristian",900,600);
  TLegend *leggy2 = new TLegend(0.3,0.2,0.5,0.5);
  leggy2->SetBorderSize(0);
  leggy2->SetFillColor(0);
  leggy2->SetFillStyle(0);  
  for(int i=0;i<numFiles;i++) {
     if(!plotThis[i]) continue;
     grAntTrans[i]->SetLineColor(getNiceColour(i));
     grAntTrans[i]->SetLineWidth(2);
     if(i==0) {
	grAntTrans[i]->Draw("al");
	grAntTrans[i]->SetTitle("S11 Measurements Fat Christian");
	grAntTrans[i]->SetMaximum(1);
	grAntTrans[i]->SetMinimum(0);
	grAntTrans[i]->GetXaxis()->SetRangeUser(0,1000);
	grAntTrans[i]->GetXaxis()->SetTitle("Frequency (MHz)");
	grAntTrans[i]->GetYaxis()->SetTitle("Transmission (Ratio)");
     }
     grAntTrans[i]->Draw("l");
     leggy2->AddEntry(grAntTrans[i],graphDesc[i],"l");


  }
  leggy2->Draw();
}

void plotS11FatChristian() {

   const int numFiles=4;
   char filename[180];
   int antFiles[numFiles][2]={{1,8},{1,11},{1,23},{1,24}};
   int whichOpen[numFiles]={1,1,2,2};
   int plotThis[numFiles]={1,1,1,1};
   char *graphDesc[numFiles]={"Office (1)","Office (2)","Outside (1)","Outside (2)"};
   TGraph *grOpen[NUM_OPENS]; 
   TGraph *grOpenPower[NUM_OPENS];
   for(int i=0;i<NUM_OPENS;i++) {
      grOpen[i]=getOpenWave(i); 
      //      grOpenPower[i]=FFTtools::makePowerSpectrumVoltsSecondsPadded(grOpen[i],8);
      grOpenPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grOpen[i],32,1024);
      Double_t *freq=grOpenPower[i]->GetX();
      //      cout << i << "\t" << grOpenPower[i]->GetN() << "\t" <<  freq[1]-freq[0] << "\n";
   }
   
   TGraph *grAnt[numFiles];
   TGraph *grAntPower[numFiles];
   TGraph *grAntRatio[numFiles];
   TGraph *grAntTrans[numFiles];
   for(int i=0;i<numFiles;i++) {
       sprintf(filename,"/home/rjn/saltStuff/hockley2009/disk_%d/TEK%05d.txt",antFiles[i][0],antFiles[i][1]);
       grAnt[i]=getWaveFromTxt(filename);
       grAnt[i]=cropLikeOpen(grAnt[i],whichOpen[i]);
       if(i<3) {
	 //	  grAntPower[i]=FFTtools::makePowerSpectrumVoltsSecondsPadded(grAnt[i],8);
	  grAntPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grAnt[i],32,1024);
       }
       else {
	  grAntPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grAnt[i],32,1024);
	  //	  TGraph *grTemp=FFTtools::makePowerSpectrumVoltsSecondsPadded(grAnt[i],8);
	  //	  grAntPower[i]=FFTtools::smoothFFT(grTemp,2);
       }
       grAntRatio[i]=FFTtools::dbGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
       grAntTrans[i]=FFTtools::ratioSubtractOneGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
       //       Double_t *freq=grAntPower[i]->GetX();
       //       cout << i << "\t" << freq[1]-freq[0] << "\n";
       //       cout << i << "\t" << grAntPower[i]->GetN() << "\t" <<  freq[1]-freq[0] << "\n";
       //       cout << grAntRatio[i] << "\t" << grAntTrans[i] << endl;
   }


  TCanvas *canWave = new TCanvas("canFatChristianWave","canFatChristianWave",800,800);
  canWave->Divide(1,numFiles);
  for(int i=0;i<numFiles;i++) {
     canWave->cd(i+1);
     grAnt[i]->SetLineColor(getNiceColour(i));     
     grAnt[i]->Draw("al");
     grAnt[i]->SetTitle(graphDesc[i]);
     grAnt[i]->GetXaxis()->SetTitle("Time (s)");
     grAnt[i]->GetYaxis()->SetTitle("Voltage (V)");

  }

  TCanvas *can = new TCanvas("canFatChristian","canFatChristian",900,600);
  TLegend *leggy = new TLegend(0.6,0.2,0.8,0.5);
  leggy->SetBorderSize(0);
  leggy->SetFillColor(0);
  leggy->SetFillStyle(0);  
  for(int i=0;i<numFiles;i++) {
     if(!plotThis[i]) continue;
     grAntRatio[i]->SetLineColor(getNiceColour(i));
     grAntRatio[i]->SetLineWidth(2);
     if(i==0) {
	grAntRatio[i]->Draw("al");
	grAntRatio[i]->SetTitle("S11 Measurements Fat Christian");
	grAntRatio[i]->SetMaximum(0);
	grAntRatio[i]->SetMinimum(-20);
	grAntRatio[i]->GetXaxis()->SetRangeUser(0,1000);
	grAntRatio[i]->GetXaxis()->SetTitle("Frequency (MHz)");
	grAntRatio[i]->GetYaxis()->SetTitle("Reflected Power (dB)");
     }
     grAntRatio[i]->Draw("l");
     leggy->AddEntry(grAntRatio[i],graphDesc[i],"l");
  }
  leggy->Draw();


  TCanvas *canTrans = new TCanvas("canTransFatChristian","canTransFatChristian",900,600);
  TLegend *leggy2 = new TLegend(0.3,0.2,0.5,0.5);
  leggy2->SetBorderSize(0);
  leggy2->SetFillColor(0);
  leggy2->SetFillStyle(0);  
  for(int i=0;i<numFiles;i++) {
     if(!plotThis[i]) continue;
     grAntTrans[i]->SetLineColor(getNiceColour(i));
     grAntTrans[i]->SetLineWidth(2);
     if(i==0) {
	grAntTrans[i]->Draw("al");
	grAntTrans[i]->SetTitle("S11 Measurements Fat Christian");
	grAntTrans[i]->SetMaximum(1);
	grAntTrans[i]->SetMinimum(0);
	grAntTrans[i]->GetXaxis()->SetRangeUser(0,1000);
	grAntTrans[i]->GetXaxis()->SetTitle("Frequency (MHz)");
	grAntTrans[i]->GetYaxis()->SetTitle("Transmission (Ratio)");
     }
     grAntTrans[i]->Draw("l");
     leggy2->AddEntry(grAntTrans[i],graphDesc[i],"l");


  }
  leggy2->Draw();
}

void plotS11SkinnyChristian() {

   const int numFiles=9;
   char filename[180];
   int antFiles[numFiles][2]={{1,8},{1,11},{1,23},{1,24},{2,3},{2,5},{5,1},{5,3},{5,2}};
   int whichOpen[numFiles]={1,1,2,2,5,5,6,6,6};
   int plotThis[numFiles]={1,1,1,1,2,2,2,2,2};
   char *graphDesc[numFiles]={"Office (1)","Office (2)","Outside (1)","Outside (2)","Salt Day I (1)","Salt Day I (2)","Salt Day III (1)","Salt Day III (2)","Salt Top of Hole (1)"};
   TGraph *grOpen[NUM_OPENS]; 
   TGraph *grOpenPower[NUM_OPENS];
   for(int i=0;i<NUM_OPENS;i++) {
      grOpen[i]=getOpenWave(i); 
      //      grOpenPower[i]=FFTtools::makePowerSpectrumVoltsSecondsPadded(grOpen[i],8);
      grOpenPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grOpen[i],32,1024);
      Double_t *freq=grOpenPower[i]->GetX();
      //      cout << i << "\t" << grOpenPower[i]->GetN() << "\t" <<  freq[1]-freq[0] << "\n";
   }
   
   TGraph *grAnt[numFiles];
   TGraph *grAntPower[numFiles];
   TGraph *grAntRatio[numFiles];
   TGraph *grAntTrans[numFiles];
   for(int i=0;i<numFiles;i++) {
       sprintf(filename,"/home/rjn/saltStuff/hockley2009/disk_%d/TEK%05d.txt",antFiles[i][0],antFiles[i][1]);
       grAnt[i]=getWaveFromTxt(filename);
       grAnt[i]=cropLikeOpen(grAnt[i],whichOpen[i]);
       if(i<3) {
	 //	  grAntPower[i]=FFTtools::makePowerSpectrumVoltsSecondsPadded(grAnt[i],8);
	  grAntPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grAnt[i],32,1024);
       }
       else {
	  grAntPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grAnt[i],32,1024);
	  //	  TGraph *grTemp=FFTtools::makePowerSpectrumVoltsSecondsPadded(grAnt[i],8);
	  //	  grAntPower[i]=FFTtools::smoothFFT(grTemp,2);
       }
       grAntRatio[i]=FFTtools::dbGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
       grAntTrans[i]=FFTtools::ratioSubtractOneGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
       //       Double_t *freq=grAntPower[i]->GetX();
       //       cout << i << "\t" << freq[1]-freq[0] << "\n";
       //       cout << i << "\t" << grAntPower[i]->GetN() << "\t" <<  freq[1]-freq[0] << "\n";
       //       cout << grAntRatio[i] << "\t" << grAntTrans[i] << endl;
   }


  TCanvas *canWave = new TCanvas("canSkinnyChristianWave","canSkinnyChristianWave",800,800);
  canWave->Divide(1,numFiles);
  for(int i=0;i<numFiles;i++) {
     canWave->cd(i+1);
     grAnt[i]->SetLineColor(getNiceColour(i));     
     grAnt[i]->Draw("al");
     grAnt[i]->SetTitle(graphDesc[i]);
     grAnt[i]->GetXaxis()->SetTitle("Time (s)");
     grAnt[i]->GetYaxis()->SetTitle("Voltage (V)");

  }

  TCanvas *can = new TCanvas("canSkinnyChristian","canSkinnyChristian",800,800);
  can->Divide(1,2); 
  {
    for(int subPad=1;subPad<=2;subPad++) {
      can->cd(subPad);
      TMultiGraph *mg = new TMultiGraph();
      TLegend *leggy = new TLegend(0.6,0.2,0.8,0.5);
      leggy->SetBorderSize(0);
      leggy->SetFillColor(0);
      leggy->SetFillStyle(0);  
      for(int i=0;i<numFiles;i++) {
	  if(plotThis[i]!=subPad) continue;
	  grAntRatio[i]->SetLineColor(getNiceColour(i));
	  grAntRatio[i]->SetLineWidth(2);
	  mg->Add(grAntRatio[i],"l");
	  leggy->AddEntry(grAntRatio[i],graphDesc[i],"l");     
      }
      mg->Draw("al");
      mg->SetTitle("S11 Measurements Skinny Christian");
      mg->SetMaximum(0);
      mg->SetMinimum(-20);
      mg->GetXaxis()->SetRangeUser(0,1000);
      mg->GetXaxis()->SetTitle("Frequency (MHz)");
      mg->GetYaxis()->SetTitle("Reflected Power (dB)");
      leggy->Draw();
    }
  }


  TCanvas *canTrans = new TCanvas("canTransSkinnyChristian","canTransSkinnyChristian",800,800);
  canTrans->Divide(1,2); 
  for(int subPad=1;subPad<=2;subPad++) {
    canTrans->cd(subPad);
    
    TLegend *leggy2 = new TLegend(0.3,0.2,0.5,0.5);
    leggy2->SetBorderSize(0);
    leggy2->SetFillColor(0);
    leggy2->SetFillStyle(0);  
    TMultiGraph *mg = new TMultiGraph();
    for(int i=0;i<numFiles;i++) {	
      if(plotThis[i]!=subPad) continue;
      grAntTrans[i]->SetLineColor(getNiceColour(i));
      grAntTrans[i]->SetLineWidth(2);
      mg->Add(grAntTrans[i],"l");
      leggy2->AddEntry(grAntTrans[i],graphDesc[i],"l");
    }
    mg->Draw("al");
    mg->SetTitle("S11 Measurements Skinny Christian");
    mg->SetMaximum(1);
    mg->SetMinimum(0);
    mg->GetXaxis()->SetRangeUser(0,1000);
    mg->GetXaxis()->SetTitle("Frequency (MHz)");
    mg->GetYaxis()->SetTitle("Transmission (Ratio)");
    leggy2->Draw(); 
  }
}





void plotS11LongPole() {

   const int numFiles=8;
   char filename[180];
   int antFiles[numFiles][2]={{1,15},{1,16},{1,29},{1,30},{2,10},{2,11},{5,6},{5,7}};
   int whichOpen[numFiles]={1,1,2,2,5,5,6,6};
   int plotThis[numFiles]={1,1,1,1,2,2,2,2};
   char *graphDesc[numFiles]={"Office (1)","Office (2)","Outside (1)","Outside (2)","Salt Day I (1)","Salt Day I (2)","Salt Day III (1)","Salt Day III (2)"};
   TGraph *grOpen[NUM_OPENS]; 
   TGraph *grOpenPower[NUM_OPENS];
   for(int i=0;i<NUM_OPENS;i++) {
      grOpen[i]=getOpenWave(i); 
      //      grOpenPower[i]=FFTtools::makePowerSpectrumVoltsSecondsPadded(grOpen[i],8);
      grOpenPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grOpen[i],32,1024);
      Double_t *freq=grOpenPower[i]->GetX();
      //      cout << i << "\t" << grOpenPower[i]->GetN() << "\t" <<  freq[1]-freq[0] << "\n";
   }
   
   TGraph *grAnt[numFiles];
   TGraph *grAntPower[numFiles];
   TGraph *grAntRatio[numFiles];
   TGraph *grAntTrans[numFiles];
   for(int i=0;i<numFiles;i++) {
       sprintf(filename,"/home/rjn/saltStuff/hockley2009/disk_%d/TEK%05d.txt",antFiles[i][0],antFiles[i][1]);
       grAnt[i]=getWaveFromTxt(filename);
       grAnt[i]=cropLikeOpen(grAnt[i],whichOpen[i]);
       if(i<3) {
	 //	  grAntPower[i]=FFTtools::makePowerSpectrumVoltsSecondsPadded(grAnt[i],8);
	  grAntPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grAnt[i],32,1024);
       }
       else {
	  grAntPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grAnt[i],32,1024);
	  //	  TGraph *grTemp=FFTtools::makePowerSpectrumVoltsSecondsPadded(grAnt[i],8);
	  //	  grAntPower[i]=FFTtools::smoothFFT(grTemp,2);
       }
       grAntRatio[i]=FFTtools::dbGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
       grAntTrans[i]=FFTtools::ratioSubtractOneGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
       //       Double_t *freq=grAntPower[i]->GetX();
       //       cout << i << "\t" << freq[1]-freq[0] << "\n";
       //       cout << i << "\t" << grAntPower[i]->GetN() << "\t" <<  freq[1]-freq[0] << "\n";
       //       cout << grAntRatio[i] << "\t" << grAntTrans[i] << endl;
   }


  TCanvas *canWave = new TCanvas("canLongPoleWave","canLongPoleWave",800,800);
  canWave->Divide(1,numFiles);
  for(int i=0;i<numFiles;i++) {
     canWave->cd(i+1);
     grAnt[i]->SetLineColor(getNiceColour(i));     
     grAnt[i]->Draw("al");
     grAnt[i]->SetTitle(graphDesc[i]);
     grAnt[i]->GetXaxis()->SetTitle("Time (s)");
     grAnt[i]->GetYaxis()->SetTitle("Voltage (V)");

  }

  TCanvas *can = new TCanvas("canLongPole","canLongPole",900,600);
  can->Divide(1,2);
  for(int subPad=1;subPad<=2;subPad++) {
    can->cd(subPad);
    TMultiGraph *mg = new TMultiGraph();
    
    TLegend *leggy = new TLegend(0.6,0.2,0.8,0.5);
    leggy->SetBorderSize(0);
    leggy->SetFillColor(0);
    leggy->SetFillStyle(0);  
    for(int i=0;i<numFiles;i++) {
      if(plotThis[i]!=subPad) continue;
      grAntRatio[i]->SetLineColor(getNiceColour(i));
      grAntRatio[i]->SetLineWidth(2);
      mg->Add(grAntRatio[i],"l");
      leggy->AddEntry(grAntRatio[i],graphDesc[i],"l");
    }
    mg->Draw("al");
    mg->SetTitle("S11 Measurements Long Pole");
    mg->SetMaximum(0);
    mg->SetMinimum(-20);
    mg->GetXaxis()->SetRangeUser(0,1000);
    mg->GetXaxis()->SetTitle("Frequency (MHz)");
    mg->GetYaxis()->SetTitle("Reflected Power (dB)");
    leggy->Draw();
  }


  TCanvas *canTrans = new TCanvas("canTransLongPole","canTransLongPole",900,600);
  canTrans->Divide(1,2);
  for(int subPad=1;subPad<=2;subPad++) {
    canTrans->cd(subPad);
    TMultiGraph *mg = new TMultiGraph();
    TLegend *leggy2 = new TLegend(0.3,0.2,0.5,0.5);
    leggy2->SetBorderSize(0);
    leggy2->SetFillColor(0);
    leggy2->SetFillStyle(0);  
    for(int i=0;i<numFiles;i++) {
      if(plotThis[i]!=subPad) continue;
      grAntTrans[i]->SetLineColor(getNiceColour(i));
      grAntTrans[i]->SetLineWidth(2);
      mg->Add(grAntTrans[i],"l");
      leggy2->AddEntry(grAntTrans[i],graphDesc[i],"l");
    }
    mg->Draw("al");
    mg->SetTitle("S11 Measurements Long Pole");
    mg->SetMaximum(1);
    mg->SetMinimum(0);
    mg->GetXaxis()->SetRangeUser(0,1000);
    mg->GetXaxis()->SetTitle("Frequency (MHz)");
    mg->GetYaxis()->SetTitle("Transmission (Ratio)");    
    leggy2->Draw();
  }

}





void plotS11Baby() {

   const int numFiles=8;
   char filename[180];
   int antFiles[numFiles][2]={{1,17},{1,18},{1,28},{1,27},{2,7},{2,8},{5,8},{5,9}};
   int whichOpen[numFiles]={1,1,2,2,5,5,6,6};
   int plotThis[numFiles]={1,1,1,1,2,2,2,2};
   char *graphDesc[numFiles]={"Office (1)","Office (2)","Outside (1)","Outside (2)","Salt Day I (1)","Salt Day I (2)","Salt Day III (1)","Salt Day III (2)"};
   TGraph *grOpen[NUM_OPENS]; 
   TGraph *grOpenPower[NUM_OPENS];
   for(int i=0;i<NUM_OPENS;i++) {
      grOpen[i]=getOpenWave(i); 
      //      grOpenPower[i]=FFTtools::makePowerSpectrumVoltsSecondsPadded(grOpen[i],8);
      grOpenPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grOpen[i],32,1024);
      Double_t *freq=grOpenPower[i]->GetX();
      //      cout << i << "\t" << grOpenPower[i]->GetN() << "\t" <<  freq[1]-freq[0] << "\n";
   }
   
   TGraph *grAnt[numFiles];
   TGraph *grAntPower[numFiles];
   TGraph *grAntRatio[numFiles];
   TGraph *grAntTrans[numFiles];
   for(int i=0;i<numFiles;i++) {
       sprintf(filename,"/home/rjn/saltStuff/hockley2009/disk_%d/TEK%05d.txt",antFiles[i][0],antFiles[i][1]);
       grAnt[i]=getWaveFromTxt(filename);
       grAnt[i]=cropLikeOpen(grAnt[i],whichOpen[i]);
       if(i<3) {
	 //	  grAntPower[i]=FFTtools::makePowerSpectrumVoltsSecondsPadded(grAnt[i],8);
	  grAntPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grAnt[i],32,1024);
       }
       else {
	  grAntPower[i]=FFTtools::makePSVSBartlettPaddedOverlap(grAnt[i],32,1024);
	  //	  TGraph *grTemp=FFTtools::makePowerSpectrumVoltsSecondsPadded(grAnt[i],8);
	  //	  grAntPower[i]=FFTtools::smoothFFT(grTemp,2);
       }
       grAntRatio[i]=FFTtools::dbGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
       grAntTrans[i]=FFTtools::ratioSubtractOneGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
       //       Double_t *freq=grAntPower[i]->GetX();
       //       cout << i << "\t" << freq[1]-freq[0] << "\n";
       //       cout << i << "\t" << grAntPower[i]->GetN() << "\t" <<  freq[1]-freq[0] << "\n";
       //       cout << grAntRatio[i] << "\t" << grAntTrans[i] << endl;
   }


  TCanvas *canWave = new TCanvas("canBabyWave","canBabyWave",800,800);
  canWave->Divide(1,numFiles);
  for(int i=0;i<numFiles;i++) {
     canWave->cd(i+1);
     grAnt[i]->SetLineColor(getNiceColour(i));     
     grAnt[i]->Draw("al");
     grAnt[i]->SetTitle(graphDesc[i]);
     grAnt[i]->GetXaxis()->SetTitle("Time (s)");
     grAnt[i]->GetYaxis()->SetTitle("Voltage (V)");

  }

  TCanvas *can = new TCanvas("canBaby","canBaby",900,600);
  can->Divide(1,2);
  for(int subPad=1;subPad<=2;subPad++) {
    can->cd(subPad);
    TMultiGraph *mg = new TMultiGraph();
    
    TLegend *leggy = new TLegend(0.2,0.2,0.4,0.5);
    leggy->SetBorderSize(0);
    leggy->SetFillColor(0);
    leggy->SetFillStyle(0);  
    for(int i=0;i<numFiles;i++) {
      if(plotThis[i]!=subPad) continue;
      grAntRatio[i]->SetLineColor(getNiceColour(i));
      grAntRatio[i]->SetLineWidth(2);
      mg->Add(grAntRatio[i],"l");
      leggy->AddEntry(grAntRatio[i],graphDesc[i],"l");
    }
    mg->Draw("al");
    mg->SetTitle("S11 Measurements Baby");
    mg->SetMaximum(0);
    mg->SetMinimum(-20);
    mg->GetXaxis()->SetRangeUser(0,1000);
    mg->GetXaxis()->SetTitle("Frequency (MHz)");
    mg->GetYaxis()->SetTitle("Reflected Power (dB)");
    leggy->Draw();
  }


  TCanvas *canTrans = new TCanvas("canTransBaby","canTransBaby",900,600);
  canTrans->Divide(1,2);
  for(int subPad=1;subPad<=2;subPad++) {
    canTrans->cd(subPad);
    TMultiGraph *mg = new TMultiGraph();
    TLegend *leggy2 = new TLegend(0.2,0.55,0.4,0.85);
    leggy2->SetBorderSize(0);
    leggy2->SetFillColor(0);
    leggy2->SetFillStyle(0);  
    for(int i=0;i<numFiles;i++) {
      if(plotThis[i]!=subPad) continue;
      grAntTrans[i]->SetLineColor(getNiceColour(i));
      grAntTrans[i]->SetLineWidth(2);
      mg->Add(grAntTrans[i],"l");
      leggy2->AddEntry(grAntTrans[i],graphDesc[i],"l");
    }
    mg->Draw("al");
    mg->SetTitle("S11 Measurements Baby");
    mg->SetMaximum(1);
    mg->SetMinimum(0);
    mg->GetXaxis()->SetRangeUser(0,1000);
    mg->GetXaxis()->SetTitle("Frequency (MHz)");
    mg->GetYaxis()->SetTitle("Transmission (Ratio)");    
    leggy2->Draw();
  }


}


TGraph *getWaveFromCsv(char *fileName)
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


TGraph *getWaveFromTxt(char *fileName)
{
  ifstream TxtFile(fileName);  
  int numLines=0;
  int numLines2=0;
  TString line;
  Double_t xVals[100000];
  Double_t yVals[100000];
  Double_t xVals2[100000];
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
  
  
TGraph *getScopePeriodiogram(char *fileName)
{
  ifstream CsvFile(fileName);
  
  int numLines=0;
  int numLines2=0;
  TString line;
  Double_t xVals[100000];
  Double_t yVals[100000];
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
  Double_t xVals[100000];
  Double_t yVals[100000];
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

