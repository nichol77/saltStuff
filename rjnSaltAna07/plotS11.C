
TGraph *getWave(char *fileName);

void plotS11()
{
  gROOT->SetStyle("BABAR");
  gSystem->Load("/sw/lib/libfftw3.so");
  gSystem->Load("libMathMore.so");
  gSystem->Load("libRootFftwWrapper.so");
  //  Gstyle->SetHistLineWidth(2);
  //  gStyle->SetLineWidth(2);
    plotS11Christian();
  //  plotS11Raven();
  //  plotS11LongPole();
  //  plotS11Baby();
  //    plotOpens();
}


void plotOpens() {
   char filename[180];
   int openFiles[2][2]={{1,0},{5,9}};
   TGraph *grOpen[2]; 
   TGraph *grOpenPower[2];
   TGraph *grOpenPowerBart[2];
   TGraph *grOpenPowerBartSeg[2];
   for(int i=0;i<2;i++) {
      sprintf(filename,"/home/rjn/saltStuff/salt_trip3_data/disk_%d/TEK%05d.CSV",openFiles[i][0],openFiles[i][1]);
      grOpen[i]=getWave(filename); 
      cout << grOpen[i]->GetN() << endl;
      grOpenPower[i]=FFTtools::makePowerSpectrumVoltsSeconds(grOpen[i]);
      grOpenPowerBart[i]=FFTtools::makePowerSpectrumVoltsSecondsBartlett(grOpen[i]);
      grOpenPowerBartSeg[i]=FFTtools::makePSVSBartlettPaddedOverlap(grOpen[i],32,1024);
   }
   TCanvas *canOpen  = new TCanvas("canOpen","canOpen",600,600);
   canOpen->Divide(1,4);
   canOpen->cd(1);
   grOpen[0]->SetTitle("Open -- Above & Below Ground");
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
   grOpenPower[0]->SetTitle("Open -- Above & Below Ground");
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

   canOpen->cd(3);
   for(int i=0;i<2;i++) {
      gPad->SetLogy();
      grOpenPowerBart[i]->SetLineColor(getNiceColour(i));
      grOpenPowerBart[i]->SetLineWidth(2);
      //      grOpenPowerBart[i]->SetTitle("Open");
      if(i==0) {
	 grOpenPowerBart[i]->Draw("al");
	 grOpenPowerBart[i]->SetMaximum(1e-12);
	 grOpenPowerBart[i]->SetMinimum(1e-16);
	 grOpenPowerBart[i]->GetXaxis()->SetRangeUser(0,1000);
	 grOpenPowerBart[i]->GetXaxis()->SetTitle("Frequency (MHz)");
	 grOpenPowerBart[i]->GetYaxis()->SetTitle("Power (normalised to area)");
      }
      else
	 grOpenPowerBart[i]->Draw("l");
   }

   canOpen->cd(4);
   for(int i=0;i<2;i++) {
      gPad->SetLogy();
      grOpenPowerBartSeg[i]->SetLineColor(getNiceColour(i));
      grOpenPowerBartSeg[i]->SetLineWidth(2);
      //      grOpenPowerBartSeg[i]->SetTitle("Open");
      if(i==0) {
	 grOpenPowerBartSeg[i]->Draw("al");
	 grOpenPowerBartSeg[i]->SetMaximum(1e-12);
	 grOpenPowerBartSeg[i]->SetMinimum(1e-16);
	 grOpenPowerBartSeg[i]->GetXaxis()->SetRangeUser(0,1000);
	 grOpenPowerBartSeg[i]->GetXaxis()->SetTitle("Frequency (MHz)");
	 grOpenPowerBartSeg[i]->GetYaxis()->SetTitle("Power (normalised to area)");
      }
      else
	 grOpenPowerBartSeg[i]->Draw("l");
   }



}


void plotS11Raven() {

   const int numFiles=3;
   char filename[180];
   int openFiles[2][2]={{1,0},{5,9}};
   int antFiles[numFiles][2]={{1,2},{1,3},{5,10}};
   int whichOpen[numFiles]={0,0,1};
   int plotThis[numFiles]={1,1,1};
   char *graphDesc[numFiles]={"Air (1)","Air (2)","Salt 50ft (1)"};
   TGraph *grOpen[2]; 
   TGraph *grOpenPower[2];
   for(int i=0;i<2;i++) {
      sprintf(filename,"/home/rjn/saltStuff/salt_trip3_data/disk_%d/TEK%05d.CSV",openFiles[i][0],openFiles[i][1]);
      grOpen[i]=getWave(filename); 
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
       sprintf(filename,"/home/rjn/saltStuff/salt_trip3_data/disk_%d/TEK%05d.CSV",antFiles[i][0],antFiles[i][1]);
       grAnt[i]=getWave(filename);
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
  canWave->Divide(1,3);
  for(int i=0;i<numFiles;i++) {
     canWave->cd(i+1);
     grAnt[i]->SetLineColor(getNiceColour(i));     
     grAnt[i]->Draw("al");
     grAnt[i]->SetTitle(graphDesc[i]);
     grAnt[i]->GetXaxis()->SetTitle("Time (s)");
     grAnt[i]->GetYaxis()->SetTitle("Voltage (V)");

  }

  TCanvas *can = new TCanvas("canRaven","canRaven",900,600);
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
	grAntRatio[i]->SetTitle("S11 Measurements Short Raven");
	grAntRatio[i]->SetMaximum(0);
	grAntRatio[i]->SetMinimum(-30);
	grAntRatio[i]->GetXaxis()->SetRangeUser(0,1000);
	grAntRatio[i]->GetXaxis()->SetTitle("Frequency (MHz)");
	grAntRatio[i]->GetYaxis()->SetTitle("Reflected Power (dB)");
     }
     grAntRatio[i]->Draw("l");
     leggy->AddEntry(grAntRatio[i],graphDesc[i],"l");
  }
  leggy->Draw();


  TCanvas *canTrans = new TCanvas("canTransRaven","canTransRaven",900,600);
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
	grAntTrans[i]->SetTitle("S11 Measurements Short Raven");
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

void plotS11Christian() {
   const int numFiles=8;
   char filename[180];
   int openFiles[2][2]={{1,0},{5,9}};
   int antFiles[numFiles][2]={{1,8},{1,10},{2,1},{2,2},{5,2},{5,3},{5,4},{5,13}};
   int whichOpen[numFiles]={0,0,1,1,1,1,1,1};
   //   int whichOpen[numFiles]={0,0,0,0,0,0,0,0};
   int plotThis[numFiles]={1,0,0,0,1,0,0,1};
   char *graphDesc[numFiles]={"Air","Air (2)","Air (1)","On Salt (1)","Salt 90ft","Salt 20ft (1)","Salt 10ft (1)","Salt 50ft"};
   TGraph *grOpen[2]; 
   TGraph *grOpenPower[2];
   for(int i=0;i<2;i++) {
      sprintf(filename,"/home/rjn/saltStuff/salt_trip3_data/disk_%d/TEK%05d.CSV",openFiles[i][0],openFiles[i][1]);
      grOpen[i]=getWave(filename);
      grOpenPower[i]=FFTtools::makePowerSpectrumVoltsSecondsPadded(grOpen[i],16);
   }
   
   TGraph *grAnt[numFiles];
   TGraph *grAntPower[numFiles];
   TGraph *grAntRatio[numFiles];
   TGraph *grAntTransMhz[numFiles];
   TGraph *grAntTrans[numFiles];
   for(int i=0;i<numFiles;i++) {
       sprintf(filename,"/home/rjn/saltStuff/salt_trip3_data/disk_%d/TEK%05d.CSV",antFiles[i][0],antFiles[i][1]);
       grAnt[i]=getWave(filename);
       //       if(i<2) 
       grAntPower[i]=FFTtools::makePowerSpectrumVoltsSecondsPadded(grAnt[i],16);
       //       else {
       //       	  TGraph *grWave=FFTtools::makePowerSpectrumVoltsSecondsPadded(grAnt[i],8);
       //       	  grAntPower[i]=FFTtools::smoothFFT(grWave,2);
       //       }

       grAntRatio[i]=FFTtools::dbGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
       grAntTrans[i]=FFTtools::ratioSubtractOneGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
       // grAntTransMhz[i]=FFTtools::ratioSubtractOneGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
//        Double_t *tVals=grAntTransMhz[i]->GetY();
//        Double_t *mhzVals=grAntTransMhz[i]->GetX();
//        Int_t nPoints=grAntTransMhz[i]->GetN();
//        Double_t *ghzVals= new Double_t [nPoints];
//        for(int ind=0;ind<nPoints;ind++) {
// 	 ghzVals[ind]=mhzVals[ind]*1e-3;
//        }
//        grAntTrans[i] = new TGraph(nPoints,ghzVals,tVals);
   }

   TGraph *grWaveRat[numFiles];
   for(int i=4;i<8;i++) 
      grWaveRat[i]=FFTtools::subtractGraphs(grAnt[i],grAnt[6]);
   
   TCanvas *canRat = new TCanvas("canSaltyRat","canSaltyRat",900,600);
   for(int i=4;i<numFiles;i++) {
      grWaveRat[i]->SetLineColor(getNiceColour(i));
      if(i==4) {
	 grWaveRat[i]->Draw("al");
	 grWaveRat[i]->SetTitle("S11 Measurements Christian (Salty)");
	 //	 grWaveRat[i]->SetMaximum(10);
	 //	 grWaveRat[i]->SetMinimum(-10);
	 //	 grWaveRat[i]->GetXaxis()->SetRangeUser(0,1000);
	 grWaveRat[i]->GetXaxis()->SetTitle("Time (s)");
	 grWaveRat[i]->GetYaxis()->SetTitle("Ratio to 10ft");
      }
      grWaveRat[i]->Draw("l");
   }
   
   //   return;

  TCanvas *canWave = new TCanvas("canSaltyWave","canSaltyWave",800,800);
  canWave->Divide(2,4);
  for(int i=0;i<numFiles;i++) {
     canWave->cd(i+1);
     grAnt[i]->SetLineColor(getNiceColour(i));
     
     grAnt[i]->Draw("al");
     grAnt[i]->SetTitle(graphDesc[i]);
     grAnt[i]->GetXaxis()->SetTitle("Time (s)");
     grAnt[i]->GetYaxis()->SetTitle("Voltage (V)");

  }

  TCanvas *canPower = new TCanvas("canSaltyPower","canSaltyPower",800,800);
  canPower->Divide(1,2);
  canPower->cd(1);
  gPad->SetLogy();
  TLegend *legPower = new TLegend(0.2,0.2,0.5,0.5);
  legPower->SetBorderSize(0);
  legPower->SetFillColor(0);
  legPower->SetFillStyle(0);  
  for(int i=0;i<3;i++) {
     if(!plotThis[i]) continue;
     grAntPower[i]->SetLineColor(getNiceColour(i));
     grAntPower[i]->SetLineWidth(2);
     if(i==0) {
	grAntPower[i]->Draw("al");
	grAntPower[i]->SetTitle("Raw Power Christian (Salty)");
	grAntPower[i]->SetMaximum(1e-12);
	grAntPower[i]->SetMinimum(1e-16);
	grAntPower[i]->GetXaxis()->SetRangeUser(0,1000);
	grAntPower[i]->GetXaxis()->SetTitle("Frequency (MHz)");
	grAntPower[i]->GetYaxis()->SetTitle("FFT Power (arb.)");
     }
     grAntPower[i]->Draw("l");
     legPower->AddEntry(grAntPower[i],graphDesc[i],"l");
  } 
  canPower->cd(2);
  gPad->SetLogy();
  for(int i=4;i<numFiles;i++) {
     if(!plotThis[i]) continue;
     grAntPower[i]->SetLineColor(getNiceColour(i));
     grAntPower[i]->SetLineWidth(2);
     if(i==4) {
	grAntPower[i]->Draw("al");
	grAntPower[i]->SetTitle("Raw Power Christian (Salty)");
	//	grAntPower[i]->SetTitle("S11 Measurements Christian (Salty)");
	grAntPower[i]->SetMaximum(1e-12);
	grAntPower[i]->SetMinimum(1e-16);
	grAntPower[i]->GetXaxis()->SetRangeUser(0,1000);
	grAntPower[i]->GetXaxis()->SetTitle("Frequency (MHz)");
	grAntPower[i]->GetYaxis()->SetTitle("FFT Power (arb.)");
	//	grAntPower[i]->GetYaxis()->SetTitle("Reflected Power (dB)");
     }
     grAntPower[i]->Draw("l");
     legPower->AddEntry(grAntRatio[i],graphDesc[i],"l");
  }
  legPower->Draw();


  TCanvas *can = new TCanvas("canSalty","canSalty",800,800);
  can->Divide(1,2);
  can->cd(1);
  TLegend *leggy = new TLegend(0.2,0.2,0.5,0.5);
  leggy->SetBorderSize(0);
  leggy->SetFillColor(0);
  leggy->SetFillStyle(0);  
  for(int i=0;i<3;i++) {
     if(!plotThis[i]) continue;
     grAntRatio[i]->SetLineColor(getNiceColour(i));
     grAntRatio[i]->SetLineWidth(2);
     if(i==0) {
	grAntRatio[i]->Draw("al");
	grAntRatio[i]->SetTitle("S11 Measurements Christian (Salty)");
	grAntRatio[i]->SetMaximum(0);
	grAntRatio[i]->SetMinimum(-30);
	grAntRatio[i]->GetXaxis()->SetRangeUser(0,1000);
	grAntRatio[i]->GetXaxis()->SetTitle("Frequency (MHz)");
	grAntRatio[i]->GetYaxis()->SetTitle("Reflected Power (dB)");
     }
     grAntRatio[i]->Draw("l");
     leggy->AddEntry(grAntRatio[i],graphDesc[i],"l");
  } 
  can->cd(2);
  for(int i=4;i<numFiles;i++) {
     if(!plotThis[i]) continue;
     grAntRatio[i]->SetLineColor(getNiceColour(i));
     grAntRatio[i]->SetLineWidth(2);
     if(i==4) {
	grAntRatio[i]->Draw("al");
	grAntRatio[i]->SetTitle("S11 Measurements Christian (Salty)");
	grAntRatio[i]->SetMaximum(0);
	grAntRatio[i]->SetMinimum(-30);
	grAntRatio[i]->GetXaxis()->SetRangeUser(0,1000);
	grAntRatio[i]->GetXaxis()->SetTitle("Frequency (MHz)");
	grAntRatio[i]->GetYaxis()->SetTitle("Reflected Power (dB)");
     }
     grAntRatio[i]->Draw("l");
     leggy->AddEntry(grAntRatio[i],graphDesc[i],"l");
  }
  leggy->Draw();


  TCanvas *canTrans = new TCanvas("canTransSalty","canTransSalty",600,400);
  //  canTrans->Divide(1,2);
  //  canTrans->cd(1);
  TLegend *leggy2 = new TLegend(0.7,0.2,0.9,0.5);
  leggy2->SetBorderSize(0);
  leggy2->SetFillColor(0);
  leggy2->SetFillStyle(0);  
  {
    TH1F *framey = gPad->DrawFrame(0,0,600,1);
    framey->GetXaxis()->SetTitle("Frequency (MHz)");
    framey->GetYaxis()->SetTitle("Transmitted Power Fraction");
  }
  int count=0;
  //  SetPaletteGreys();
  
  for(int i=0;i<numFiles;i++) {
     if(!plotThis[i]) continue;
     grAntTrans[i]->SetLineColor(1);
     //     grAntTrans[i]->SetLineColor(kGray+count);
     grAntTrans[i]->SetLineWidth(2);
     grAntTrans[i]->SetLineStyle(2*count+1);
     count++;
     grAntTrans[i]->Draw("l");
     leggy2->AddEntry(grAntTrans[i],graphDesc[i],"l");
  //    char outName[180];
//      sprintf(outName,"christianAntennaTransmission_%d.txt",i);
//      ofstream OutFile(outName);
//      Double_t *freq = grAntTrans[i]->GetX();
//      Double_t *tran = grAntTrans[i]->GetY();
//      Int_t N=grAntTrans[i]->GetN();
//      for(int i=0;i<N;i++) {
//        OutFile << freq[i] << "\t" << tran[i] << "\n";
//      }
//      OutFile.close();
  }
//   canTrans->cd(2);  
//   {
//     TH1F *framey = gPad->DrawFrame(0,0,600,1);
//     framey->GetXaxis()->SetTitle("Frequency (MHz)");
//     framey->GetYaxis()->SetTitle("Transmision (Ratio)");
//   }
//   for(int i=4;i<numFiles;i++) {
//      if(!plotThis[i]) continue;
//      grAntTrans[i]->SetLineColor(getNiceColour(i));
//      grAntTrans[i]->SetLineWidth(2);
// //      if(i==4) {
// // 	grAntTrans[i]->Draw("al");
// // 	grAntTrans[i]->SetTitle("S11 Measurements Christian (Salty)");
// // 	grAntTrans[i]->SetMaximum(1);
// // 	grAntTrans[i]->SetMinimum(0);
// // 	grAntTrans[i]->GetXaxis()->SetRangeUser(0,1000);
// // 	grAntTrans[i]->GetXaxis()->SetTitle("Frequency (MHz)");
// // 	grAntTrans[i]->GetYaxis()->SetTitle("Transmission (Ratio)");
// //      }
//      grAntTrans[i]->Draw("l");
//      leggy2->AddEntry(grAntTrans[i],graphDesc[i],"l");

//      char outName[180];
//      sprintf(outName,"christianAntennaTransmission_%d.txt",i);
//      ofstream OutFile(outName);
//      Double_t *freq = grAntTrans[i]->GetX();
//      Double_t *tran = grAntTrans[i]->GetY();
//      Int_t N=grAntTrans[i]->GetN();
//      for(int ind=0;ind<N;ind++) {
//        OutFile << freq[ind] << "\t" << tran[ind] << "\n";
//      }
//      OutFile.close();
//   }
  leggy2->Draw();

}



void plotS11LongPole() {
   const int numFiles=4;
   char filename[180];
   int openFiles[2][2]={{1,0},{5,9}};
   int antFiles[numFiles][2]={{1,14},{1,16},{5,11},{5,12}};
   int whichOpen[numFiles]={0,0,1,1};
   int plotThis[numFiles]={1,1,1,1};
   char *graphDesc[numFiles]={"Air (1)","Air (2)","Salt 20ft (1)","Salt 75ft (1)"};
   TGraph *grOpen[2]; 
   TGraph *grOpenPower[2];
   for(int i=0;i<2;i++) {
      sprintf(filename,"/home/rjn/saltStuff/salt_trip3_data/disk_%d/TEK%05d.CSV",openFiles[i][0],openFiles[i][1]);
      grOpen[i]=getWave(filename);
      grOpenPower[i]=FFTtools::makePowerSpectrumVoltsSecondsPadded(grOpen[i],8);
      //      cout << grOpen[i]->GetN() << "\t" << grOpenPower[i]->GetN() << "\n";
   }
   
   TGraph *grAnt[numFiles];
   TGraph *grAntPower[numFiles];
   TGraph *grAntRatio[numFiles];
   TGraph *grAntTrans[numFiles];
   for(int i=0;i<numFiles;i++) {
       sprintf(filename,"/home/rjn/saltStuff/salt_trip3_data/disk_%d/TEK%05d.CSV",antFiles[i][0],antFiles[i][1]);
       grAnt[i]=getWave(filename);      
       if(i<2) {      
	  TGraph *grTemp=FFTtools::makePowerSpectrumVoltsSecondsPadded(grAnt[i],8);
	  grAntPower[i]=FFTtools::smoothFFT(grTemp,2);
	  delete grTemp;
       }
       else { 
	  grAntPower[i]=FFTtools::makePowerSpectrumVoltsSecondsPadded(grAnt[i],8);
       }
       grAntRatio[i]=FFTtools::dbGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
       grAntTrans[i]=FFTtools::ratioSubtractOneGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
       //       cout << grAnt[i]->GetN() << "\t" << grAntPower[i]->GetN() << "\t" << grAntRatio[i] << "\n";
   }


  TCanvas *canWave = new TCanvas("canLongPoleWave","canLongPoleWave",600,800);
  canWave->Divide(1,4);
  for(int i=0;i<numFiles;i++) {
     canWave->cd(i+1);
     grAnt[i]->SetLineColor(getNiceColour(i));   
     grAnt[i]->Draw("al");
     grAnt[i]->SetTitle(graphDesc[i]);
     grAnt[i]->GetXaxis()->SetTitle("Time (s)");
     grAnt[i]->GetYaxis()->SetTitle("Voltage (V)");

  }
  TCanvas *canPower = new TCanvas("canLongPolePower","canLongPolePower",800,800);
  canPower->Divide(1,4);
  for(int i=0;i<numFiles;i++) {
     canPower->cd(i+1);
     grAntPower[i]->SetLineColor(getNiceColour(i));     
     grAntPower[i]->Draw("al");
     grAntPower[i]->SetTitle(graphDesc[i]);
     grAntPower[i]->GetXaxis()->SetTitle("Frequency (MHz)");
     //     grAnt[i]->GetYaxis()->SetTitle("");

  }

  TCanvas *can = new TCanvas("canLongPole","canLongPole",900,600);
  TLegend *leggy = new TLegend(0.65,0.2,0.85,0.4);
  leggy->SetBorderSize(0);
  leggy->SetFillColor(0);
  leggy->SetFillStyle(0);  
  for(int i=0;i<numFiles;i++) {
     if(!plotThis[i]) continue;
     grAntRatio[i]->SetLineColor(getNiceColour(i));
     grAntRatio[i]->SetLineWidth(2);
     if(i==0) {
	grAntRatio[i]->Draw("al");
	grAntRatio[i]->SetTitle("S11 Measurements Long Pole");
	grAntRatio[i]->SetMaximum(0);
	grAntRatio[i]->SetMinimum(-30);
	grAntRatio[i]->GetXaxis()->SetRangeUser(0,500);
	grAntRatio[i]->GetXaxis()->SetTitle("Frequency (MHz)");
	grAntRatio[i]->GetYaxis()->SetTitle("Reflected Power (dB)");
     }
     grAntRatio[i]->Draw("l");
     leggy->AddEntry(grAntRatio[i],graphDesc[i],"l");
  }
  leggy->Draw();


  TCanvas *canTrans = new TCanvas("canTransLongPole","canTransLongPole",900,600);
  TLegend *leggy2 = new TLegend(0.23,0.2,0.43,0.4);
  leggy2->SetBorderSize(0);
  leggy2->SetFillColor(0);
  leggy2->SetFillStyle(0);  
  for(int i=0;i<numFiles;i++) {
     if(!plotThis[i]) continue;
     grAntTrans[i]->SetLineColor(getNiceColour(i));
     grAntTrans[i]->SetLineWidth(2);
     if(i==0) {
	grAntTrans[i]->Draw("al");
	grAntTrans[i]->SetTitle("S11 Measurements Long Pole");
	grAntTrans[i]->SetMaximum(1);
	grAntTrans[i]->SetMinimum(0);
	grAntTrans[i]->GetXaxis()->SetRangeUser(0,500);
	grAntTrans[i]->GetXaxis()->SetTitle("Frequency (MHz)");
	grAntTrans[i]->GetYaxis()->SetTitle("Transmission (Ratio)");
     }
     grAntTrans[i]->Draw("l");
     leggy2->AddEntry(grAntTrans[i],graphDesc[i],"l");
  }
  leggy2->Draw();

}





void plotS11Baby() {
   const int numFiles=4;
   char filename[180];
   int openFiles[2][2]={{1,0},{5,9}};
   int antFiles[numFiles][2]={{1,18},{1,20},{5,7},{5,8}};
   int whichOpen[numFiles]={0,0,1,1};
   int plotThis[numFiles]={1,1,1,1};
   char *graphDesc[numFiles]={"Air (1)","Air (2)","Salt 10ft (1)","Salt 90ft (1)"};
   TGraph *grOpen[2]; 
   TGraph *grOpenPower[2];
   for(int i=0;i<2;i++) {
      sprintf(filename,"/home/rjn/saltStuff/salt_trip3_data/disk_%d/TEK%05d.CSV",openFiles[i][0],openFiles[i][1]);
      grOpen[i]=getWave(filename);
      grOpenPower[i]=FFTtools::makePowerSpectrumVoltsSeconds(grOpen[i]);
      //      cout << grOpen[i]->GetN() << "\t" << grOpenPower[i]->GetN() << "\n";
   }
   
   TGraph *grAnt[numFiles];
   TGraph *grAntPower[numFiles];
   TGraph *grAntRatio[numFiles];
   TGraph *grAntTrans[numFiles];
   for(int i=0;i<numFiles;i++) {
       sprintf(filename,"/home/rjn/saltStuff/salt_trip3_data/disk_%d/TEK%05d.CSV",antFiles[i][0],antFiles[i][1]);
       if(i<2) {
	  grAnt[i]=getWave(filename);       
	  TGraph *grTemp=FFTtools::makePowerSpectrumVoltsSeconds(grAnt[i]);
	  grAntPower[i]=FFTtools::smoothFFT(grTemp,2);
	  delete grTemp;
       }
       else {
	  grAnt[i]=0;
	  grAntPower[i]=getScopeFFTNormToArea(filename);
       }
       grAntRatio[i]=FFTtools::dbGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
       grAntTrans[i]=FFTtools::ratioSubtractOneGraphs(grAntPower[i],grOpenPower[whichOpen[i]]);
   }


  TCanvas *canPower = new TCanvas("canBabyPower","canBabyPower",800,800);
  canPower->Divide(1,6);
  canPower->cd(5);
  grOpenPower[0]->Draw("al");
  canPower->cd(6);
  grOpenPower[1]->Draw("al");
  for(int i=0;i<numFiles;i++) {
     canPower->cd(i+1);
     grAntPower[i]->SetLineColor(getNiceColour(i));     
     grAntPower[i]->Draw("al");
     grAntPower[i]->SetTitle(graphDesc[i]);
     grAntPower[i]->GetXaxis()->SetTitle("Frequency (MHz)");
     //     grAnt[i]->GetYaxis()->SetTitle("");

  }

  TCanvas *can = new TCanvas("canBaby","canBaby",900,600);
  TLegend *leggy = new TLegend(0.65,0.2,0.85,0.4);
  leggy->SetBorderSize(0);
  leggy->SetFillColor(0);
  leggy->SetFillStyle(0);  
  for(int i=0;i<numFiles;i++) {
     if(!plotThis[i]) continue;
     grAntRatio[i]->SetLineColor(getNiceColour(i));
     grAntRatio[i]->SetLineWidth(2);
     if(i==0) {
	grAntRatio[i]->Draw("al");
	grAntRatio[i]->SetTitle("S11 Measurements Baby");
	grAntRatio[i]->SetMaximum(0);
	grAntRatio[i]->SetMinimum(-30);
	grAntRatio[i]->GetXaxis()->SetRangeUser(0,1000);
	grAntRatio[i]->GetXaxis()->SetTitle("Frequency (MHz)");
	grAntRatio[i]->GetYaxis()->SetTitle("Reflected Power (dB)");
     }
     grAntRatio[i]->Draw("l");
     leggy->AddEntry(grAntRatio[i],graphDesc[i],"l");
  }
  leggy->Draw();


  TCanvas *canTrans = new TCanvas("canTransBaby","canTransBaby",900,600);
  TLegend *leggy2 = new TLegend(0.23,0.2,0.43,0.4);
  leggy2->SetBorderSize(0);
  leggy2->SetFillColor(0);
  leggy2->SetFillStyle(0);  
  for(int i=0;i<numFiles;i++) {
     if(!plotThis[i]) continue;
     grAntTrans[i]->SetLineColor(getNiceColour(i));
     grAntTrans[i]->SetLineWidth(2);
     if(i==0) {
	grAntTrans[i]->Draw("al");
	grAntTrans[i]->SetTitle("S11 Measurements Baby");
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
