#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TBox.h"
#include "TMinuit.h"
#include "TF1.h"
#include "TF2.h"
#include "TTree.h"
#include "TFile.h"
#include "TLegend.h"

using namespace std;

TStyle* RootStyle();
TStyle *color=RootStyle();

float calcAtten(float power1, float power2, float dist1, float dist2);
                                                 //calculates the attenuation length

int main() {

  gStyle=color;

  //variables to read out of the rootfile
  const int HALF_MAX_ARRAY_SIZE=16384; //max array size in the time domain
  float f[HALF_MAX_ARRAY_SIZE]; //frequency
  float t[HALF_MAX_ARRAY_SIZE]; //time
  float vfreq[HALF_MAX_ARRAY_SIZE]; //v in the frequency domain (scaled by sqrt (dt*N))
  float vtime[HALF_MAX_ARRAY_SIZE]; //v in the time domain
  float power[HALF_MAX_ARRAY_SIZE]; //power in the frequency domain (scaled by dt*N)
  float maxvoltage,maxvoltage_freq, maxpower_freq;//maximum voltage in the time domain, 
                                    //frequency domain, and max power.  Used for plots.
  int idisk,ifile,itimebins,ifreqbins, filenumber; //disk number, nth file on disk, 
                                          //number of time bins, number of freq bins, 
                                          //and filenumber that appears in the filename.
  float v_max,v_min,vrms;//vmin and vmax in time domain for peak to peak

  float f_noise[HALF_MAX_ARRAY_SIZE];
  float t_noise[HALF_MAX_ARRAY_SIZE];
  float power_noise[HALF_MAX_ARRAY_SIZE];
  int idisk_noise,ifile_noise,itimebins_noise,ifreqbins_noise, filenumber_noise;
                                               //disk number, nth file on disk,      
                                             //number of time bins, number of freq bins, 
                                          //and filenumber that appears in the filename.
  //for the cut in frequency
  float christian_min=125e6;//low frequency edge in Hz for lowest band of antenna
  float christian_max=175e6;//high frequency edge in Hz for lowest band of antenna
  float long_min=100e6;
  float long_max=150e6;
  float baby_min=175e6;
  float baby_max=225e6;
  float christian_bandStep=50e6;
  float long_bandStep=50e6;
  float baby_bandStep=100e6;

  //variables for calculation of attenuation length
  int nantennas=3, ndepths=9, nbands=9, nmeasurements=2, nbandslong=2;
  float powerBand[nantennas][ndepths][nbands][nmeasurements];//power in band using FFT. 
                                    //[antennatype (0=long, 1=christian, 2=baby]
                                    //[depthindex, 0=10/20, 1=50, 2=90]
                                    //[frequency band number][0=close, 1=far hole]
  float powerBand_noise[nantennas][nbands][nmeasurements];
  float vpk2pk[nantennas][ndepths][nmeasurements];//vpeak2peak, with no index for frequency
  float hole1=163*.3048;//in meters
  float hole2=553*.3048;//in meters
  float len_atten[nantennas][ndepths][nbands];//attenuation length from FFT.  
                              //First index is for different antennas (long,christian,baby) 
                      //Second index is for depth (10/20,50,90), third is for frequency band
  float len_atten_avg[nantennas][nbands];
                              //avg attenuation length (for a given frequency bin)
  float len_atten_noise[nantennas][ndepths][nbands];
  float vpk2pk_atten[nantennas][ndepths];//attenuation length from Vpk2pk
  
  //for calculation of errors
  float vpk2pk_error_plus[nantennas][ndepths];//for error bars
  float vpk2pk_error_minus[nantennas][ndepths];
  float len_error_plus[nantennas][ndepths][nbands];
  float len_error_minus[nantennas][ndepths][nbands];
  float len_error_plus_avg[nantennas][nbands];//error bars for avg plot
  float len_error_minus_avg[nantennas][nbands];//error bars for average plot
  float sample_fluct=.054;//from spreadsheet, error in voltage 
                          //due to sampling of waveforms
  float volt_dev[ndepths][nmeasurements];//due to deviation in voltage
  float volt_dev_mean[nmeasurements], sum_volt_dev=0;
  float rms_dev_from_mean;
  float d12_uncer[9]={0.006,0.008,0.019,0.006,0.006,0.007,0.010,0.013,0.016};
                                                                 //uncertainty in d12, frac.
  float d13_uncer[9]={0.010,0.011,0.021,0.010,0.010,0.010,0.013,0.015,0.017};//fractional
  //in terms of pct of attenuation length
  float len_dev[nantennas][ndepths][nbands];//avg deviation for a frequency over all depths 
  float vpk2pk_dev[nantennas][ndepths];//avg deviation of attenlength
  float d12_uncer_atten[nantennas][ndepths][nbands];
  float d12_uncer_p2p[nantennas][ndepths];
  float d13_uncer_atten[nantennas][ndepths][nbands];
  float d13_uncer_p2p[nantennas][ndepths];
  float sample_fluct_atten[nantennas][ndepths][nbands];
  float sample_fluct_p2p[nantennas][ndepths];
  //for totals
  float total_uncer_plus[nantennas][ndepths][nbands], total_uncer_minus[nantennas][ndepths][nbands];
  float vpk2pk_total_uncer_plus[nantennas][ndepths],vpk2pk_total_uncer_minus[nantennas][ndepths];

  //toggles
  int includeshallow=0;// include 10 and 20' points in plots or not
  int logscale=1; //logscale on plots or not

  //initializing stuff
  for (int i=0; i<nantennas;i++){
    for (int j=0;j<ndepths;j++){
      for (int ctr=0;ctr<nbands;ctr++){
	for (int ctr1=0;ctr1<nmeasurements;ctr1++){
	  powerBand[i][j][ctr][ctr1]=0;
	  powerBand_noise[i][ctr][ctr1]=0;
	  vpk2pk[i][j][ctr1]=0;
	}
      }
    }
  }//end initialization
  
  //root variables
  TFile *f1=new TFile("waveforms.root"); //file with waveform information in  it
  TTree *t1=(TTree*)f1->Get("waveformtree");
  TFile *f2=new TFile("waveforms_noise.root"); //file with waveform information in  it
  TTree *t2=(TTree*)f2->Get("waveformtree");
  
  
  
  //read in waveform tree
  t1->SetBranchAddress("f",&f);//frequency
  t1->SetBranchAddress("t",&t);//time
  t1->SetBranchAddress("vfreq",&vfreq);//v in freq domain
  t1->SetBranchAddress("vtime",&vtime);//v in time domain
  t1->SetBranchAddress("idisk",&idisk); //disk number
  t1->SetBranchAddress("ifile",&ifile);//file number
  t1->SetBranchAddress("itimebins",&itimebins);//number of time bins
  t1->SetBranchAddress("ifreqbins",&ifreqbins);//number of frequency bins
  t1->SetBranchAddress("maxvoltage",&maxvoltage); //max voltage in time domain
  t1->SetBranchAddress("maxvoltage_freq",&maxvoltage_freq); 
                                               //max voltage in frequency domain
  t1->SetBranchAddress("maxpower_freq",&maxpower_freq); 
                                               //max power in frequency domain
  t1->SetBranchAddress("power",&power);//power in freq domain
  t1->SetBranchAddress("filenumber",&filenumber);//filenumber in the title of the file
  t1->SetBranchAddress("v_max",&v_max);
  t1->SetBranchAddress("v_min",&v_min);
  t1->SetBranchAddress("vrms",&vrms);
  
  int nentries=t1->GetEntries(); //get number of entries
  
  //LOOP THRU the entries in the tree
  for (int i=0;i<nentries;i++) {
    t1->GetEvent(i); //get an entry from the tree
    
    //pick out the files we are interested in, and for each one 
    //cut in frequency, add the power, and get the vpk2pk
    for (int ctr=0; ctr<nbands;ctr++){//loop through different frequency bands

      if (idisk==2 && ifile==10) {//50' christian antennas T12
	for (int i=0;i<ifreqbins;i++){//loop through all frequencies
	  if (f[i]>=christian_min+ctr*christian_bandStep && 
	      f[i]<=christian_max+ctr*christian_bandStep)//cut on freq band
	    powerBand[1][1][ctr][0]+=power[i];//add up the power in that band
	}
	vpk2pk[1][1][0]=v_max-v_min; //set the vpk2pk
      }// end 50' christian T12
      if (idisk==3 && ifile==6) {//50' christian antennas T13
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=christian_min+ctr*christian_bandStep && 
	      f[i]<=christian_max+ctr*christian_bandStep)
	    powerBand[1][1][ctr][1]+=power[i];
	}
	vpk2pk[1][1][1]=v_max-v_min;
      }//end 50' christain T13
      if (idisk==2 && ifile==31) {//50' longpole T12
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=long_min+ctr*long_bandStep && f[i]<=long_max+ctr*long_bandStep)
	    powerBand[0][1][ctr][0]+=power[i];
	}
	vpk2pk[0][1][0]=v_max-v_min;
      }// end 50' long T12
      if (idisk==3 && ifile==23) {//50' longpole T13
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=long_min+ctr*long_bandStep && f[i]<=long_max+ctr*long_bandStep)
	    powerBand[0][1][ctr][1]+=power[i];
	}
	vpk2pk[0][1][1]=v_max-v_min;
      }// end 50' long T13
      
      if (idisk==2 && ifile==26) {//50' baby T12
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=baby_min+ctr*baby_bandStep && f[i]<=baby_max+ctr*baby_bandStep)
	    powerBand[2][1][ctr][0]+=power[i];
	}
	vpk2pk[2][1][0]=v_max-v_min;
      }// end 50' baby T12
      if (idisk==3 && ifile==19) {//50' baby T13
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=baby_min+ctr*baby_bandStep && f[i]<=baby_max+ctr*baby_bandStep)
	    powerBand[2][1][ctr][1]+=power[i]*pow(10,(-34/10.));//account for amplifier
	}
	vpk2pk[2][1][1]=(v_max-v_min)*pow(10,(-34/20.));
      }// end 50' baby T13
      if (idisk==2 && ifile==19) {//90' christian antennas T12
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=christian_min+ctr*christian_bandStep && 
	      f[i]<=christian_max+ctr*christian_bandStep)
	    powerBand[1][2][ctr][0]+=power[i];
	}
	vpk2pk[1][2][0]=v_max-v_min;
      }// end 90' christian T12
      
      if (idisk==3 && ifile==14) {//90' christian antennas T13
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=christian_min+ctr*christian_bandStep 
	      && f[i]<=christian_max+ctr*christian_bandStep)
	    powerBand[1][2][ctr][1]+=power[i];
	}
	vpk2pk[1][2][1]=v_max-v_min;
      }//end 90' christain T13
      if (idisk==2 && ifile==32) {//75' long antennas T12
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=long_min+ctr*long_bandStep && f[i]<=long_max+ctr*long_bandStep)
	    powerBand[0][2][ctr][0]+=power[i];
	}
	vpk2pk[0][2][0]=v_max-v_min;
      }// end 75' long T12
      if (idisk==3 && ifile==24) {//75' long antennas T13
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=long_min+ctr*long_bandStep && f[i]<=long_max+ctr*long_bandStep)
	    powerBand[0][2][ctr][1]+=power[i];
	}
	vpk2pk[0][2][1]=v_max-v_min;
      }// end 75' long T13
      if (idisk==2 && ifile==28) {//90' baby T12
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=baby_min+ctr*baby_bandStep && f[i]<=baby_max+ctr*baby_bandStep)
	    powerBand[2][2][ctr][0]+=power[i];
	}
	vpk2pk[2][2][0]=v_max-v_min;
      }// end 90' baby T12
      if (idisk==3 && ifile==20) {//90' baby T13
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=baby_min+ctr*baby_bandStep && f[i]<=baby_max+ctr*baby_bandStep)
	    powerBand[2][2][ctr][1]+=power[i]*pow(10,(-34/10.));//account for amplifier
	}     
	vpk2pk[2][2][1]=(v_max-v_min)*pow(10,(-34/20.));
      }// end 90' baby T13
      if (idisk==2 && ifile==4) {//20' christian antennas T12
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=christian_min+ctr*christian_bandStep && 
	      f[i]<=christian_max+ctr*christian_bandStep) 
	    powerBand[1][0][ctr][0]+=power[i];
	} 
	vpk2pk[1][0][0]=v_max-v_min;
      }// end 20' christian T12    
      if (idisk==3 && ifile==1) {//20' christian antennas T13
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=christian_min+ctr*christian_bandStep && 
	      f[i]<=christian_max+ctr*christian_bandStep)
	    powerBand[1][0][ctr][1]+=power[i];
	}
	vpk2pk[1][0][1]=v_max-v_min;
    }//end 20' christain T13
      if (idisk==2 && ifile==33) {//20' long antennas T12
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=long_min+ctr*long_bandStep && f[i]<=long_max+ctr*long_bandStep)
	    powerBand[0][0][ctr][0]+=power[i];
	}
	vpk2pk[0][0][0]=v_max-v_min;
      }// end 20' long T12
      if (idisk==3 && ifile==22) {//20' long antennas T13
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=long_min+ctr*long_bandStep && f[i]<=long_max+ctr*long_bandStep)
	    powerBand[0][0][ctr][1]+=power[i];
	}
	vpk2pk[0][0][1]=v_max-v_min;
      }//end 20' long T13
      if (idisk==2 && ifile==25) {//10' baby T12
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=baby_min+ctr*baby_bandStep && f[i]<=baby_max+ctr*baby_bandStep)
	    powerBand[2][0][ctr][0]+=power[i];
	}
	vpk2pk[2][0][0]=v_max-v_min;
      }// end 10' baby T12
      if (idisk==3 && ifile==18) {//10' baby T13
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=baby_min+ctr*baby_bandStep && f[i]<=baby_max+ctr*baby_bandStep)
	    powerBand[2][0][ctr][1]+=power[i]*pow(10,(-34/10.));//account for amplifier
	}
	vpk2pk[2][0][1]=(v_max-v_min)*pow(10,(-34/20.));
      }// end 10' baby T13
      
      //other heights of christian antennas
      if (idisk==2 && ifile==2) {//10' christian antennas T12
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=christian_min+ctr*christian_bandStep && 
	      f[i]<=christian_max+ctr*christian_bandStep)
	    powerBand[1][3][ctr][0]+=power[i];
	} 
	vpk2pk[1][3][0]=v_max-v_min;
      }// end 10' christian T12
      if (idisk==3 && ifile==0) {//10' christian antennas T13
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=christian_min+ctr*christian_bandStep && 
	      f[i]<=christian_max+ctr*christian_bandStep) 
	    powerBand[1][3][ctr][1]+=power[i];
	} 
	vpk2pk[1][3][1]=v_max-v_min;
      }//end 10' christain T13
      if (idisk==2 && ifile==6) {//30' christian antennas T12
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=christian_min+ctr*christian_bandStep && 
	      f[i]<=christian_max+ctr*christian_bandStep)
	    powerBand[1][4][ctr][0]+=power[i];
	}
	vpk2pk[1][4][0]=v_max-v_min;
      }// end 30' christian T12
      if (idisk==3 && ifile==3) {//30' christian antennas T13
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=christian_min+ctr*christian_bandStep && 
	      f[i]<=christian_max+ctr*christian_bandStep)
	    powerBand[1][4][ctr][1]+=power[i];
	}
	vpk2pk[1][4][1]=v_max-v_min;
      }//end 30' christain T13
      if (idisk==2 && ifile==9) {//40' christian antennas T12
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=christian_min+ctr*christian_bandStep && 
	      f[i]<=christian_max+ctr*christian_bandStep)
	    powerBand[1][5][ctr][0]+=power[i];
	} 
	vpk2pk[1][5][0]=v_max-v_min;
      }// end 40' christian T12
      if (idisk==3 && ifile==4) {//40' christian antennas T13
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=christian_min+ctr*christian_bandStep && 
	      f[i]<=christian_max+ctr*christian_bandStep) 
	    powerBand[1][5][ctr][1]+=power[i];
	}
	vpk2pk[1][5][1]=v_max-v_min;
      }//end 40' christain T13
      if (idisk==2 && ifile==11) {//60' christian antennas T12
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=christian_min+ctr*christian_bandStep && 
	      f[i]<=christian_max+ctr*christian_bandStep)
	    powerBand[1][6][ctr][0]+=power[i];
	} 
	vpk2pk[1][6][0]=v_max-v_min;
      }// end 60' christian T12
      if (idisk==3 && ifile==8) {//60' christian antennas T13
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=christian_min+ctr*christian_bandStep && 
	      f[i]<=christian_max+ctr*christian_bandStep)
	    powerBand[1][6][ctr][1]+=power[i];
	}
	vpk2pk[1][6][1]=v_max-v_min;
      }//end 60' christain T13
      if (idisk==2 && ifile==14) {//70' christian antennas T12
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=christian_min+ctr*christian_bandStep && 
	      f[i]<=christian_max+ctr*christian_bandStep)
	    powerBand[1][7][ctr][0]+=power[i];
	} 
	vpk2pk[1][7][0]=v_max-v_min;
      }// end 70' christian T12
      if (idisk==3 && ifile==10) {//70' christian antennas T13
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=christian_min+ctr*christian_bandStep && 
	      f[i]<=christian_max+ctr*christian_bandStep)
	    powerBand[1][7][ctr][1]+=power[i];
	}
	vpk2pk[1][7][1]=v_max-v_min;
      }//end 70' christain T13
      if (idisk==2 && ifile==19) {//80' christian antennas T12
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=christian_min+ctr*christian_bandStep && 
	      f[i]<=christian_max+ctr*christian_bandStep) 
	    powerBand[1][8][ctr][0]+=power[i];
	} 
	vpk2pk[1][8][0]=v_max-v_min;
      }// end 80' christian T12
      if (idisk==3 && ifile==12) {//80' christian antennas T13
	for (int i=0;i<ifreqbins;i++){
	  if (f[i]>=christian_min+ctr*christian_bandStep && 
	      f[i]<=christian_max+ctr*christian_bandStep)
	  powerBand[1][8][ctr][1]+=power[i];
	}
	vpk2pk[1][8][1]=v_max-v_min;
      }//end 80' christain T13
    }//end loop of different frequencies for each measurement
  } // end loop over entries in the tree.  Done reading stuff!!!
  
  //read in noise waveforms
  t2->SetBranchAddress("f",&f_noise);//frequency
  t2->SetBranchAddress("t",&t_noise);//time
  t2->SetBranchAddress("idisk",&idisk_noise); //disk number
  t2->SetBranchAddress("ifile",&ifile_noise);//file number
  t2->SetBranchAddress("itimebins",&itimebins_noise);//number of time bins
  t2->SetBranchAddress("ifreqbins",&ifreqbins_noise);//number of frequency bins
  t2->SetBranchAddress("power",&power_noise);//power in freq domain
  t2->SetBranchAddress("filenumber",&filenumber_noise);//filenumber in the title of the file
  
  nentries=t2->GetEntries(); //get number of entries
  for (int i=0;i<nentries;i++) {
    t2->GetEvent(i); //get an entry from the tree
    for (int ctr=0; ctr<nbands;ctr++){//loop through different frequency bands
      if (idisk_noise==2 && ifile_noise==10) {
	for (int i=0;i<ifreqbins_noise;i++){//loop through all frequencies
	  if (f_noise[i]>=christian_min+ctr*christian_bandStep && 
	      f_noise[i]<=christian_max+ctr*christian_bandStep)//cut on freq band
	    powerBand_noise[1][ctr][0]+=power_noise[i];//add up the power in that band
	}
      }
      if (idisk_noise==2 && ifile_noise==31) {
	for (int i=0;i<ifreqbins_noise;i++){//loop through all frequencies
	  if (f_noise[i]>=long_min+ctr*long_bandStep && 
	      f_noise[i]<=long_max+ctr*long_bandStep)//cut on freq band
	    powerBand_noise[0][ctr][0]+=power_noise[i];//add up the power in that band
	}
      }
      if (idisk_noise==2 && ifile_noise==26) {
	for (int i=0;i<ifreqbins_noise;i++){//loop through all frequencies
	  if (f_noise[i]>=baby_min+ctr*baby_bandStep && 
	      f_noise[i]<=baby_max+ctr*baby_bandStep)//cut on freq band
	    powerBand_noise[2][ctr][0]+=power_noise[i]*12./38.;//add up the power in that band
	}
      }
      if (idisk_noise==3 && ifile_noise==04) {
	for (int i=0;i<ifreqbins_noise;i++){//loop through all frequencies
	  if (f_noise[i]>=christian_min+ctr*christian_bandStep && 
	      f_noise[i]<=christian_max+ctr*christian_bandStep)//cut on freq band
	    powerBand_noise[1][ctr][1]+=power_noise[i];//add up the power in that band
	}
      }
      if (idisk_noise==3 && ifile_noise==23) {
	for (int i=0;i<ifreqbins_noise;i++){//loop through all frequencies
	  if (f_noise[i]>=long_min+ctr*long_bandStep && 
	      f_noise[i]<=long_max+ctr*long_bandStep)//cut on freq band
	    powerBand_noise[0][ctr][1]+=power_noise[i];//add up the power in that band
	}
      }
      if (idisk_noise==3 && ifile_noise==19) {
	for (int i=0;i<ifreqbins_noise;i++){//loop through all frequencies
	  if (f_noise[i]>=baby_min+ctr*baby_bandStep && 
	      f_noise[i]<=baby_max+ctr*baby_bandStep)//cut on freq band
	    powerBand_noise[2][ctr][1]+=power_noise[i]*pow(10,(-34/10.));//add up power in that band
	}
      }
    }
  }//finish reading noise tree
  



  //calculate the attenuation lengths and then errors
  for (int i=0;i<nantennas;i++){
    for (int ctr=0;ctr<nbands;ctr++){
      for (int j=0;j<ndepths;j++){
	len_atten[i][j][ctr]=calcAtten(powerBand[i][j][ctr][0],
				       powerBand[i][j][ctr][1],hole1,hole2);
                                                //calc atten length w noise
	//subtract the noise
	powerBand[i][j][ctr][0]=powerBand[i][j][ctr][0]-powerBand_noise[i][ctr][0];
	powerBand[i][j][ctr][1]=powerBand[i][j][ctr][1]-powerBand_noise[i][ctr][1];
	
	len_atten_noise[i][j][ctr]=calcAtten(powerBand[i][j][ctr][0],//calc atten length wo noise
					     powerBand[i][j][ctr][1],hole1,hole2);
	len_atten_noise[i][j][ctr]=100.*(len_atten[i][j][ctr]
					 -len_atten_noise[i][j][ctr])/len_atten[i][j][ctr]; 
	                                            //calculate % difference with and wo noise
	//noise switch.  comment the line below out to not subtract noise.
	len_atten[i][j][ctr]=calcAtten(powerBand[i][j][ctr][0],
				       powerBand[i][j][ctr][1],hole1,hole2);
	//now make this attn length wo noise
	//calculate average attenuation length at a given frequency for plotting
	if (i==1 && j!=0 && j!=3) len_atten_avg[i][ctr]+=len_atten[i][j][ctr]; //christian's, don't avg in 10 and 20'
	if ((i==0 || i==2) && (j==1 || j==2)) len_atten_avg[i][ctr]+=len_atten[i][j][ctr];
	//don't avg shallow
	vpk2pk_atten[i][j]=calcAtten(vpk2pk[i][j][0]*vpk2pk[i][j][0],
				     vpk2pk[i][j][1]*vpk2pk[i][j][1],hole1,hole2);
	//calculate deviation error
	//if (j!=0 && i==1 && ctr==5) {
	if (j!=0 && i==1 && ctr==0) {
	  for (int k=0;k<nmeasurements;k++){
	    volt_dev[j][k]=(vpk2pk[i][j][k]-vpk2pk[i][j-1][k])/
	      (0.5*(vpk2pk[i][j][k]+vpk2pk[i][j-1][k]));//deviation for each point from next one
	    //volt_dev[j][k][ctr]=(sqrt(powerBand[i][j][ctr][k])-sqrt(powerBand[i][j-1][ctr][k]))/
	    //(0.5*(sqrt(powerBand[i][j][ctr][k])+sqrt(powerBand[i][j-1][ctr][k])));
	    volt_dev_mean[k]+=volt_dev[j][k];//start calculating mean deviation
	  }
	}
      }
      if (i==1) len_atten_avg[i][ctr]=len_atten_avg[i][ctr]/(ndepths-2);//calculate avg from sum
      if (i!=1) len_atten_avg[i][ctr]=len_atten_avg[i][ctr]/2.;
    }
  }
  //still doing deviation error
  for (int k=0;k<nmeasurements;k++){
      volt_dev_mean[k]=volt_dev_mean[k]/(ndepths-1.);//calculate the mean deviation
  }
  for (int k=0;k<nmeasurements;k++){
    for (int j=1;j<ndepths;j++){
      volt_dev[j][k]=volt_dev[j][k]-volt_dev_mean[k];//deviation from the mean
      volt_dev[j][k]=pow(volt_dev[j][k],2);//square deviation from the mean
      sum_volt_dev+=volt_dev[j][k];//sum the squares
    }
  }
  rms_dev_from_mean=sqrt(sum_volt_dev/(nmeasurements*(ndepths-1)));//take the rms
  
  //calculate total errors
  for (int i=0; i<nantennas;i++){
    for (int ctr=0;ctr<nbands;ctr++){
      for (int j=0;j<ndepths;j++){
	
	d12_uncer_atten[i][j][ctr]=fabs(((1-d12_uncer[j])*hole1-hole2)
					/len_atten[i][j][ctr]/
					(log(1-d12_uncer[j])+(hole1-hole2)/
					 len_atten[i][j][ctr])-1);//due hole1-2 distce
	d13_uncer_atten[i][j][ctr]=fabs((-1.*(1-d13_uncer[j])*hole2+hole1)/
					len_atten[i][j][ctr]/
					(log(1-d13_uncer[j])+(hole1-hole2)/
					 len_atten[i][j][ctr])-1);//hole 1-3
	sample_fluct_atten[i][j][ctr]=fabs(1./(len_atten[i][j][ctr]//due waveformsampling
					       *log(1+sample_fluct)/
					       (hole1-hole2)+1)-1);
	len_dev[i][j][ctr]=fabs(1./(len_atten[i][j][ctr]//due waveformsampling
				    *log(1+rms_dev_from_mean)/
				    (hole1-hole2)+1)-1);
	total_uncer_plus[i][j][ctr]=sqrt(pow(d12_uncer_atten[i][j][ctr],2)+
					 pow(len_dev[i][j][ctr],2)
					 +pow(d13_uncer_atten[i][j][ctr],2));
					 //+pow(sample_fluct_atten[i][j][ctr],2));//summed in quad
	d12_uncer_atten[i][j][ctr]=fabs(((1+d12_uncer[j])*hole1-hole2)
					/len_atten[i][j][ctr]/
					(log(1+d12_uncer[j])+(hole1-hole2)/
					 len_atten[i][j][ctr])-1);//due hole1-2 distce
	d13_uncer_atten[i][j][ctr]=fabs((-1.*(1+d13_uncer[j])*hole2+hole1)/
					len_atten[i][j][ctr]/
					(log(1+d13_uncer[j])+(hole1-hole2)/
					 len_atten[i][j][ctr])-1);//hole 1-3
	sample_fluct_atten[i][j][ctr]=fabs(1./(len_atten[i][j][ctr]//due waveformsampling
					       *log(1-sample_fluct)/
					       (hole1-hole2)+1)-1);
	len_dev[i][j][ctr]=fabs(1./(len_atten[i][j][ctr]//due waveformsampling
				    *log(1-rms_dev_from_mean)/
				    (hole1-hole2)+1)-1);
	total_uncer_minus[i][j][ctr]=sqrt(pow(d12_uncer_atten[i][j][ctr],2)+
					  pow(len_dev[i][j][ctr],2)
					  +pow(d13_uncer_atten[i][j][ctr],2));
	                                  //+pow(sample_fluct_atten[i][j][ctr],2));//summed in quad
	

	if (ctr==0){//calc total errors on atten length for vpk2pk
	  d12_uncer_p2p[i][j]=fabs(((1-d12_uncer[j])*hole1-hole2)/
				   vpk2pk_atten[i][j]/
				   (log(1-d12_uncer[j])+(hole1-hole2)/
				    vpk2pk_atten[i][j])-1);//to hole 1-2 distce
	  d13_uncer_p2p[i][j]=fabs((-1.*(1-d13_uncer[j])*hole2+hole1)/
				   vpk2pk_atten[i][j]/
				   (log(1-d13_uncer[j])+(hole1-hole2)/
				    vpk2pk_atten[i][j])-1);//hole 1-3
	  sample_fluct_p2p[i][j]=fabs(1./(vpk2pk_atten[i][j]//wavefo sampling
					  *log(1+sample_fluct)/(hole1-hole2)+1)-1);
	  vpk2pk_dev[i][j]=fabs(1./(vpk2pk_atten[i][j]//due waveformsampling
				    *log(1+rms_dev_from_mean)/
				    (hole1-hole2)+1)-1);
	  vpk2pk_total_uncer_plus[i][j]=sqrt(pow(vpk2pk_dev[i][j],2)+pow(d12_uncer_p2p[i][j],2)+
					     pow(d13_uncer_p2p[i][j],2));
	                                      //+pow(sample_fluct_p2p[i][j],2));//sum in quad
	  d12_uncer_p2p[i][j]=fabs(((1+d12_uncer[j])*hole1-hole2)/
				   vpk2pk_atten[i][j]/
				   (log(1+d12_uncer[j])+(hole1-hole2)/
				    vpk2pk_atten[i][j])-1);//to hole 1-2 distce
	  d13_uncer_p2p[i][j]=fabs((-1.*(1+d13_uncer[j])*hole2+hole1)/
				   vpk2pk_atten[i][j]/
				   (log(1+d13_uncer[j])+(hole1-hole2)/
				    vpk2pk_atten[i][j])-1);//hole 1-3
	  sample_fluct_p2p[i][j]=fabs(1./(vpk2pk_atten[i][j]//wavefo sampling
					  *log(1-sample_fluct)/(hole1-hole2)+1)-1);
	  vpk2pk_dev[i][j]=fabs(1./(vpk2pk_atten[i][j]//due waveformsampling
				    *log(1-rms_dev_from_mean)/
				    (hole1-hole2)+1)-1);
	  vpk2pk_total_uncer_minus[i][j]=sqrt(pow(vpk2pk_dev[i][j],2)+pow(d12_uncer_p2p[i][j],2)+
					      pow(d13_uncer_p2p[i][j],2));
					      //+pow(sample_fluct_p2p[i][j],2));
	  
	  
	  vpk2pk_error_plus[i][j]=vpk2pk_atten[i][j]*vpk2pk_total_uncer_plus[i][j];
	  vpk2pk_error_minus[i][j]=vpk2pk_atten[i][j]*vpk2pk_total_uncer_minus[i][j]; 
	}
	//calculate actual error from percentage
	len_error_plus[i][j][ctr]=len_atten[i][j][ctr]*total_uncer_plus[i][j][ctr];
	len_error_minus[i][j][ctr]=len_atten[i][j][ctr]*total_uncer_minus[i][j][ctr];
	if (i==1 && j!=0 && j!=3){//christian's, no shallow
	  len_error_plus_avg[i][ctr]+=len_atten[i][j][ctr]*total_uncer_plus[i][j][ctr];
	  len_error_minus_avg[i][ctr]+=len_atten[i][j][ctr]*total_uncer_minus[i][j][ctr];
	}
	if (i!=1 && (j==1 || j==2)){//baby.long, no shallow
	  len_error_plus_avg[i][ctr]+=len_atten[i][j][ctr]*total_uncer_plus[i][j][ctr];
	  len_error_minus_avg[i][ctr]+=len_atten[i][j][ctr]*total_uncer_minus[i][j][ctr];
	}	
      }
      if (i==1){
      	len_error_plus_avg[i][ctr]=len_error_plus_avg[i][ctr]/(ndepths-2);
      	len_error_minus_avg[i][ctr]=len_error_minus_avg[i][ctr]/(ndepths-2);
      }
      if (i!=1){
      	len_error_plus_avg[i][ctr]=len_error_plus_avg[i][ctr]/2.;
      	len_error_minus_avg[i][ctr]=len_error_minus_avg[i][ctr]/2.;
      }
    }
  }

  //print out what we've learned.
  cout<<endl;
  //print out attenuation lengths at all frequency bands per antenna.
  cout<<"Attenuation lengths for shallow,middle,deep depths: "<<endl;
  cout<<"Christian: "<<endl;
  cout<<"Vpk2pk: "<<vpk2pk_atten[1][3]<<", "<<vpk2pk_atten[1][0]<<", "
      <<vpk2pk_atten[1][4]<<", "<<vpk2pk_atten[1][5]<<", "
      <<vpk2pk_atten[1][1]<<", "<<vpk2pk_atten[1][6]<<", "
      <<vpk2pk_atten[1][7]<<", "<<vpk2pk_atten[1][8]<<", "
      <<vpk2pk_atten[1][2]<<endl;
  for (int ctr=0;ctr<nbands;ctr++){
    cout<<((christian_max-christian_min)/2+christian_min+ctr*christian_bandStep)/1.e6
	<<" MHz: "<<len_atten[1][3][ctr]<<", "<<len_atten[1][0][ctr]<<", "
	<<len_atten[1][4][ctr]<<", "<<len_atten[1][5][ctr]<<", "
	<<len_atten[1][1][ctr]<<", "<<len_atten[1][6][ctr]<<", "
	<<len_atten[1][7][ctr]<<", "<<len_atten[1][8][ctr]<<", "
	<<len_atten[1][2][ctr]<<endl;
  }
  cout<<"Longpole: "<<endl;
  cout<<"Vpk2pk: "<<vpk2pk_atten[0][0]<<", "<<vpk2pk_atten[0][1]<<", "<<
    vpk2pk_atten[0][2]<<endl;
  for (int ctr=0;ctr<nbands;ctr++){
    cout<<((long_max-long_min)/2+long_min+ctr*long_bandStep)/1.e6<<" MHz: "<<
      len_atten[0][0][ctr]<<", "<<len_atten[0][1][ctr]<<", "<<len_atten[0][2][ctr]<<endl;
  }
  cout<<"Baby: "<<endl;
  cout<<"Vpk2pk: "<<vpk2pk_atten[2][0]<<", "<<vpk2pk_atten[2][1]<<", "<<
    vpk2pk_atten[2][2]<<endl;
  for (int ctr=0;ctr<nbands;ctr++){
    cout<<((baby_max-baby_min)/2+baby_min+ctr*baby_bandStep)/1.e6<<" MHz: "<<
      len_atten[2][0][ctr]<<", "<<len_atten[2][1][ctr]<<", "<<len_atten[2][2][ctr]<<endl;
  }
  cout<<endl;
  //print out uncertainty values
  cout<<"rms_dev_from_mean (voltage): "<<rms_dev_from_mean<<endl;
  cout<<"sample_fluct (voltage): "<<sample_fluct<<endl;
  for (int i=0; i<ndepths;i++)
    cout<<"i: "<<i<<" d12_uncer (distance): "<<d12_uncer[i]<<" d13_uncer (distance): "
	<<d13_uncer[i]<<endl;
  
  

  //MAKE SOME PLOTS
  TGraphAsymmErrors *g1, *g1b, *g1c, *g1d, *g1e, *g1f, *g1g; //define graphs for depth plot
  TGraphAsymmErrors *g2, *g2b, *g2c, *g2d, *g2e, *g2f, *g2g, *g2h, *g2i; 
                                                        //define graphs for frequency plot
  TGraphAsymmErrors *g3,*g3b, *g3c, *g3d, *g3e, *g3f,*g3g, *g3h,*g3i; //define graphs for 
                                                        //christian antennas frequency plot
  TGraphAsymmErrors *gdepth; //depth plot for paper
  TGraph *gfreq, *gavgchristian, *gavgbaby, *gavglong, *gavg;//freq plot for paper

  TCanvas *c1=new TCanvas("c1","c1",1000,1600);
  TCanvas *c2=new TCanvas("c2","c2",1000,1600);
  TCanvas *cpaper1=new TCanvas("cpaper1","cpaper1",1000,800);
  TCanvas *cpaper2=new TCanvas("cpaper2","cpaper2",1000,800);
  TCanvas *cpaper3=new TCanvas("cpaper3","cpaper3",1000,800);
  TLegend *leg1, *leg2;//legends
  TLegend *leg3, *legavg, *legdepth;
  c1->Divide(1,2);
  c2->Divide(1,2);
  
  //depth plot...rearranging the arrays so I can plot them
  float depth[9]={10,20,30,40,50,60,70,80,90};
                                          //depth in feet corresponding to array indeces
  float atten_forplot_depth[nbands][ndepths];
  float atten_p2p_depth[ndepths];
  float plus_forplot_depth[nbands][ndepths], minus_forplot_depth[nbands][ndepths]; 
  float plus_vp2p_depth[ndepths], minus_vp2p_depth[ndepths];
  for (int ctr=0;ctr<nbands;ctr++){//rearrange so the arrays are in depth order for plots
    atten_forplot_depth[ctr][0]=len_atten[1][3][ctr];
    atten_p2p_depth[0]=vpk2pk_atten[1][3];
    atten_forplot_depth[ctr][1]=len_atten[1][0][ctr];
    atten_p2p_depth[1]=vpk2pk_atten[1][0];
    atten_forplot_depth[ctr][2]=len_atten[1][4][ctr];
    atten_p2p_depth[2]=vpk2pk_atten[1][4];
    atten_forplot_depth[ctr][3]=len_atten[1][5][ctr];
    atten_p2p_depth[3]=vpk2pk_atten[1][5];
    atten_forplot_depth[ctr][4]=len_atten[1][1][ctr];
    atten_p2p_depth[4]=vpk2pk_atten[1][1];
    atten_forplot_depth[ctr][5]=len_atten[1][6][ctr];
    atten_p2p_depth[5]=vpk2pk_atten[1][6];
    atten_forplot_depth[ctr][6]=len_atten[1][7][ctr];
    atten_p2p_depth[6]=vpk2pk_atten[1][7];
    atten_forplot_depth[ctr][7]=len_atten[1][8][ctr];
    atten_p2p_depth[7]=vpk2pk_atten[1][8];
    atten_forplot_depth[ctr][8]=len_atten[1][2][ctr];
    atten_p2p_depth[8]=vpk2pk_atten[1][2];

    plus_forplot_depth[ctr][0]=len_error_plus[1][3][ctr];
    minus_forplot_depth[ctr][0]=len_error_minus[1][3][ctr];
    plus_vp2p_depth[0]=vpk2pk_error_plus[1][3];
    minus_vp2p_depth[0]=vpk2pk_error_minus[1][3];
    plus_forplot_depth[ctr][1]=len_error_plus[1][0][ctr];
    minus_forplot_depth[ctr][1]=len_error_minus[1][0][ctr];
    plus_vp2p_depth[1]=vpk2pk_error_plus[1][0];
    minus_vp2p_depth[1]=vpk2pk_error_minus[1][0];
    plus_forplot_depth[ctr][2]=len_error_plus[1][4][ctr];
    minus_forplot_depth[ctr][2]=len_error_minus[1][4][ctr];
    plus_vp2p_depth[2]=vpk2pk_error_plus[1][4];
    minus_vp2p_depth[2]=vpk2pk_error_minus[1][4];
    plus_forplot_depth[ctr][3]=len_error_plus[1][5][ctr];
    minus_forplot_depth[ctr][3]=len_error_minus[1][5][ctr];
    plus_vp2p_depth[3]=vpk2pk_error_plus[1][5];
    minus_vp2p_depth[3]=vpk2pk_error_minus[1][5];
    plus_forplot_depth[ctr][4]=len_error_plus[1][1][ctr];
    minus_forplot_depth[ctr][4]=len_error_minus[1][1][ctr];
    plus_vp2p_depth[4]=vpk2pk_error_plus[1][1];
    minus_vp2p_depth[4]=vpk2pk_error_minus[1][1];
    plus_forplot_depth[ctr][5]=len_error_plus[1][6][ctr];
    minus_forplot_depth[ctr][5]=len_error_minus[1][6][ctr];
    plus_vp2p_depth[5]=vpk2pk_error_plus[1][6];
    minus_vp2p_depth[5]=vpk2pk_error_minus[1][6];
    plus_forplot_depth[ctr][6]=len_error_plus[1][7][ctr];
    minus_forplot_depth[ctr][6]=len_error_minus[1][7][ctr];
    plus_vp2p_depth[6]=vpk2pk_error_plus[1][7];
    minus_vp2p_depth[6]=vpk2pk_error_minus[1][7];
    plus_forplot_depth[ctr][7]=len_error_plus[1][8][ctr];
    minus_forplot_depth[ctr][7]=len_error_minus[1][8][ctr];
    plus_vp2p_depth[7]=vpk2pk_error_plus[1][8];
    minus_vp2p_depth[7]=vpk2pk_error_minus[1][8];
    plus_forplot_depth[ctr][8]=len_error_plus[1][2][ctr];
    minus_forplot_depth[ctr][8]=len_error_minus[1][2][ctr];
    plus_vp2p_depth[8]=vpk2pk_error_plus[1][2];
    minus_vp2p_depth[8]=vpk2pk_error_minus[1][2];
  }
  
  //define the depth graphs
  g1b=new TGraphAsymmErrors(ndepths,depth,atten_forplot_depth[1], 
			    0,0,minus_forplot_depth[1], plus_forplot_depth[1]);//bin0 freq
  g1=new TGraphAsymmErrors(ndepths,depth,atten_p2p_depth, 
			   0,0,minus_vp2p_depth, plus_vp2p_depth);//vpk2pk
  g1c=new TGraphAsymmErrors(ndepths,depth,atten_forplot_depth[2], 
			    0,0,minus_forplot_depth[2], plus_forplot_depth[2]);//bin1 freq
  g1d=new TGraphAsymmErrors(ndepths,depth,atten_forplot_depth[3], 
			    0,0,minus_forplot_depth[3], plus_forplot_depth[3]);//bin2 freq
  g1e=new TGraphAsymmErrors(ndepths,depth,atten_forplot_depth[4], 
			    0,0,minus_forplot_depth[4], plus_forplot_depth[4]);
  g1f=new TGraphAsymmErrors(ndepths,depth,atten_forplot_depth[6], 
			    0,0,minus_forplot_depth[6], plus_forplot_depth[6]);
  g1g=new TGraphAsymmErrors(ndepths,depth,atten_forplot_depth[8], 
			    0,0,minus_forplot_depth[8], plus_forplot_depth[8]);
  
  //now get rid of shallow for paper plot
  float atten_paperdepth[ndepths-2], plus_paperdepth[ndepths-2], minus_paperdepth[ndepths-2];
  float paper_depth[ndepths-2];
  for (int k=0;k<ndepths-2;k++){
    //choosing second index for paper plot here, 250 MHz
    atten_paperdepth[k]=atten_forplot_depth[2][k+2];
    plus_paperdepth[k]=plus_forplot_depth[2][k+2];
    minus_paperdepth[k]=minus_forplot_depth[2][k+2];
    paper_depth[k]=depth[k+2];
  }
  gdepth=new TGraphAsymmErrors(ndepths-2,paper_depth,atten_paperdepth,
			       0,0,minus_paperdepth, plus_paperdepth);
                               //use second index for paper plot: 250 MHz
  TH2F *h1=new TH2F("Field Attenuation length Versus Depth (ft) 1",
		    ""
		    ,100,0,100,100,0,300);
  h1->SetXTitle("Depth (ft)");
  h1->SetYTitle("Field Attenuation Length (m)");
  
  TH2F *hdepth=new TH2F("Field Attenuation length (m) versus Depth (ft) 2",
		    ""
		    ,100,20,100,100,0,200);
  hdepth->SetXTitle("Depth (ft)");
  hdepth->SetYTitle("Field Attenuation Length (m)");

  //add a legend
  leg1= new TLegend(0.65,0.6,.84,0.90);
  legdepth=new TLegend(0.2,0.75,0.4,0.85);
  leg1->AddEntry(g1,"Vpk2pk","l");
  leg1->AddEntry(g1b,"200 MHz","l");
  leg1->AddEntry(g1c,"250 MHz","l");
  leg1->AddEntry(g1d,"300 MHz","l");
  leg1->AddEntry(g1e,"350 MHz","l");
  leg1->AddEntry(g1f,"450 MHz","l");
  leg1->AddEntry(g1g,"550 MHz","l");
  legdepth->AddEntry(gdepth,"250 MHz","P");
  //frequency plot
  float freq_forplot[nantennas][nbands];
  for (int ctr=0;ctr<nbands;ctr++) freq_forplot[0][ctr]=((long_max-long_min)/2.
						    +long_min)+ctr*long_bandStep;
  for (int ctr=0;ctr<nbands;ctr++) freq_forplot[1][ctr]=((christian_max-christian_min)/2.
						    +christian_min)+ctr*christian_bandStep;
  for (int ctr=0;ctr<nbands;ctr++) freq_forplot[2][ctr]=((baby_max-baby_min)/2.
						    +baby_min)+ctr*baby_bandStep;

 
  //define the frequency graphs
  g2=new TGraphAsymmErrors(nbandslong,freq_forplot[0],len_atten[0][0],0,0,
			   len_error_minus[0][0],len_error_plus[0][0]);//long shallow
  g2b=new TGraphAsymmErrors(nbandslong,freq_forplot[0],len_atten[0][1],0,0,
			   len_error_minus[0][1],len_error_plus[0][1]);//long mid
  g2c=new TGraphAsymmErrors(nbandslong,freq_forplot[0],len_atten[0][2],0,0,
			   len_error_minus[0][2],len_error_plus[0][2]);//long deep
  g2d=new TGraphAsymmErrors(nbands,freq_forplot[1],len_atten[1][0],0,0,
			   len_error_minus[1][0],len_error_plus[1][0]);//christ shallow
  g2e=new TGraphAsymmErrors(nbands,freq_forplot[1],len_atten[1][1],0,0,
			   len_error_minus[1][1],len_error_plus[1][1]);//christ mid
  g2f=new TGraphAsymmErrors(nbands,freq_forplot[1],len_atten[1][2],0,0,
			   len_error_minus[1][2],len_error_plus[1][2]);//christ deep
  g2g=new TGraphAsymmErrors(nbands,freq_forplot[2],len_atten[2][0],0,0,
			   len_error_minus[2][0],len_error_plus[2][0]);//baby shallow
  g2h=new TGraphAsymmErrors(nbands,freq_forplot[2],len_atten[2][1],0,0,
			   len_error_minus[2][1],len_error_plus[2][1]);//baby mid
  g2i=new TGraphAsymmErrors(nbands,freq_forplot[2],len_atten[2][2],0,0,
			    len_error_minus[2][2],len_error_plus[2][2]);//baby deep
  TH2F *h2=new TH2F("Field Attenuation length (m) versus frequency",
		    "",
		    100,9e7,1.05e9,100,7,300);
                    //100,0,1.05e9,100,0,500);
  
  h2->SetXTitle("Frequency (Hz)");
  h2->SetYTitle("Field Attenuation Length (m)");
  
  //add a 1/f line
  int nx=100;
  float linex[nx],liney[nx];
  for (int i=0;i<nx;i++){
    linex[i]=(i+1)*10.e6;//10 MHz steps
    liney[i]=1./(i+1)*75.*30.;
  }
  TGraph *gline;
  gline=new TGraph(nx, linex,liney);

  //add a fit to all the data. (for attenuation length  vs. frequency plot)
  float atten_forfit[(nbands+nbandslong+nbands)*2];//total number of points on the plot.
  for (int i=0;i<(nbands+nbandslong+nbands)*2;i++) atten_forfit[i]=0;
  float freq_forfit[(nbands+nbandslong+nbands)*2]; //not including shallow points in fit
  for (int i=0;i<(nbands+nbandslong+nbands)*2;i++) freq_forfit[i]=0;
  for (int i=0;i<nantennas;i++){ //loop over antennas
    for (int j=0;j<2;j++){//loop over depths
      for (int ctr=0;ctr<nbands;ctr++){ //loop over bands
	if (i==0 && ctr<2){
	  atten_forfit[j*nbandslong+ctr]=len_atten[i][j+1][ctr];
	  freq_forfit[j*nbandslong+ctr]=freq_forplot[i][ctr];
	    }
	else if (i!=0){
	  atten_forfit[((i-1)*2+j)*nbands+ctr+nbandslong*2]=len_atten[i][j+1][ctr];
	  freq_forfit[((i-1)*2+j)*nbands+ctr+nbandslong*2]=freq_forplot[i][ctr];
	}
      }
    }
  }
  gfreq=new TGraph(nbands*2*2+nbandslong*2, freq_forfit,atten_forfit);//define graph
  TF1 *ffreq=new TF1("ffreq","[0]*pow(x,[1])",1e8,1e9);//define fit
  ffreq->SetParNames("Constant","Power Law");
  gfreq->Fit(ffreq);
  
  //add a legend
  if (logscale!=1) leg2= new TLegend(0.6,0.5,.80,0.90);
  if (logscale==1) leg2= new TLegend(0.18,0.18,.50,0.52);
  if (includeshallow==1) leg2->AddEntry(g2,"LongPoles, 20'","P");
  leg2->AddEntry(g2b,"LF Antennas, 50'","P");
  leg2->AddEntry(g2c,"LF Antennas, 75'","P");
  if (includeshallow==1) leg2->AddEntry(g2d,"Christian's, 20'","P");
  leg2->AddEntry(g2e,"Midband Antennas, 50'","P");
  leg2->AddEntry(g2f,"Midband Antennas, 90'","P");
  if (includeshallow==1) leg2->AddEntry(g2g,"Babies, 10'","P");
  leg2->AddEntry(g2h,"HF Antennas, 50'","P");
  leg2->AddEntry(g2i,"HF Antennas, 90'","P");
  //leg2->AddEntry(gline,"1/f, arbitrary scaling","l");
  
  //define the frequency christian graphs
  g3=new TGraphAsymmErrors(nbands,freq_forplot[1],len_atten[1][3],0,0,
			   len_error_minus[1][3],len_error_plus[1][3]);//christ 10'
  g3b=new TGraphAsymmErrors(nbands,freq_forplot[1],len_atten[1][0],0,0,
			    len_error_minus[1][0],len_error_plus[1][0]);//christ 20'
  g3c=new TGraphAsymmErrors(nbands,freq_forplot[1],len_atten[1][4],0,0,
			    len_error_minus[1][4],len_error_plus[1][4]);//christ 30'
  g3d=new TGraphAsymmErrors(nbands,freq_forplot[1],len_atten[1][5],0,0,
			    len_error_minus[1][5],len_error_plus[1][5]);//christ 40'
  g3e=new TGraphAsymmErrors(nbands,freq_forplot[1],len_atten[1][1],0,0,
			    len_error_minus[1][1],len_error_plus[1][1]);//christ 50'
  g3f=new TGraphAsymmErrors(nbands,freq_forplot[1],len_atten[1][6],0,0,
			    len_error_minus[1][6],len_error_plus[1][6]);//christ 60'
  g3g=new TGraphAsymmErrors(nbands,freq_forplot[1],len_atten[1][7],0,0,
			    len_error_minus[1][7],len_error_plus[1][7]);//christ 70'
  g3h=new TGraphAsymmErrors(nbands,freq_forplot[1],len_atten[1][8],0,0,
			    len_error_minus[1][8],len_error_plus[1][8]);//christ 80'
  g3i=new TGraphAsymmErrors(nbands,freq_forplot[1],len_atten[1][2],0,0,
			    len_error_minus[1][2],len_error_plus[1][2]);//christ 90'

  //make one array to plot for avg frequency plot.
  float avg_atten_forplot[16];
  float avg_freq_forplot[16];
  float avg_plus[16], avg_minus[16];
  avg_atten_forplot[0]=len_atten_avg[0][0];
  avg_freq_forplot[0]=freq_forplot[0][0];
  avg_plus[0]=len_error_plus_avg[0][0];
  avg_minus[0]=len_error_minus_avg[0][0];
  avg_atten_forplot[1]=len_atten_avg[1][0];
  avg_freq_forplot[1]=freq_forplot[1][0];
  avg_plus[1]=len_error_plus_avg[1][0];
  avg_minus[1]=len_error_minus_avg[1][0];
  avg_atten_forplot[2]=len_atten_avg[0][1];
  avg_freq_forplot[2]=freq_forplot[0][1];
  avg_plus[2]=len_error_plus_avg[0][1];
  avg_minus[2]=len_error_minus_avg[0][1];
  avg_atten_forplot[3]=(len_atten_avg[1][1]+len_atten_avg[2][0])/2.;
  avg_freq_forplot[3]=freq_forplot[1][1];
  avg_plus[3]=(len_error_plus_avg[1][1]+len_error_plus_avg[2][0])/2.;
  avg_minus[3]=(len_error_minus_avg[1][1]+len_error_minus_avg[2][0])/2.;
  avg_atten_forplot[4]=len_atten_avg[1][2];
  avg_freq_forplot[4]=freq_forplot[1][2];
  avg_plus[4]=len_error_plus_avg[1][2];
  avg_minus[4]=len_error_minus_avg[1][2];
  avg_atten_forplot[5]=(len_atten_avg[1][3]+len_atten_avg[2][1])/2.;
  avg_freq_forplot[5]=freq_forplot[1][3];
  avg_plus[5]=(len_error_plus_avg[1][3]+len_error_plus_avg[2][1])/2.;
  avg_minus[5]=(len_error_minus_avg[1][3]+len_error_minus_avg[2][1])/2.;
  avg_atten_forplot[6]=len_atten_avg[1][4];
  avg_freq_forplot[6]=freq_forplot[1][4];
  avg_plus[6]=len_error_plus_avg[1][4];
  avg_minus[6]=len_error_minus_avg[1][4];
  avg_atten_forplot[7]=(len_atten_avg[1][5]+len_atten_avg[2][2])/2.;
  avg_freq_forplot[7]=freq_forplot[1][5];
  avg_plus[7]=(len_error_plus_avg[1][5]+len_error_plus_avg[2][2])/2.;
  avg_minus[7]=(len_error_minus_avg[1][5]+len_error_minus_avg[2][2])/2.;
  avg_atten_forplot[8]=len_atten_avg[1][6];
  avg_freq_forplot[8]=freq_forplot[1][6];
  avg_plus[8]=len_error_plus_avg[1][6];
  avg_minus[8]=len_error_minus_avg[1][6];
  avg_atten_forplot[9]=(len_atten_avg[1][7]+len_atten_avg[2][3])/2.;
  avg_freq_forplot[9]=freq_forplot[1][7];
  avg_plus[9]=(len_error_plus_avg[1][7]+len_error_plus_avg[2][3])/2.;
  avg_minus[9]=(len_error_minus_avg[1][7]+len_error_minus_avg[2][3])/2.;
  avg_atten_forplot[10]=len_atten_avg[1][8];
  avg_freq_forplot[10]=freq_forplot[1][8];
  avg_plus[10]=len_error_plus_avg[1][8];
  avg_minus[10]=len_error_minus_avg[1][8];
  avg_atten_forplot[11]=len_atten_avg[2][4];
  avg_freq_forplot[11]=freq_forplot[2][4];
  avg_plus[11]=len_error_plus_avg[2][4];
  avg_minus[11]=len_error_minus_avg[2][4];
  avg_atten_forplot[12]=len_atten_avg[2][5];
  avg_freq_forplot[12]=freq_forplot[2][5];
  avg_plus[12]=len_error_plus_avg[2][5];
  avg_minus[12]=len_error_minus_avg[2][5];
  avg_atten_forplot[13]=len_atten_avg[2][6];
  avg_freq_forplot[13]=freq_forplot[2][6];
  avg_plus[13]=len_error_plus_avg[2][6];
  avg_minus[13]=len_error_minus_avg[2][6];
  avg_atten_forplot[14]=len_atten_avg[2][7];
  avg_freq_forplot[14]=freq_forplot[2][7];
  avg_plus[14]=len_error_plus_avg[2][7];
  avg_minus[14]=len_error_minus_avg[2][7];
  avg_atten_forplot[15]=len_atten_avg[2][8];
  avg_freq_forplot[15]=freq_forplot[2][8];
  avg_plus[15]=len_error_plus_avg[2][8];
  avg_minus[15]=len_error_minus_avg[2][8];

  //better yet... make error bars based on the scatter between same frequency points
  for (int i=0;i<16;i++){
    avg_plus[i]=avg_atten_forplot[i]*((len_atten_avg[1][1]-len_atten_avg[2][0])/len_atten_avg[1][1]+
				      (len_atten_avg[1][3]-len_atten_avg[2][1])/len_atten_avg[1][3]+
				      (len_atten_avg[1][5]-len_atten_avg[2][2])/len_atten_avg[1][5]+
				      (len_atten_avg[1][7]-len_atten_avg[2][3])/len_atten_avg[1][7])/4.;
    avg_minus[i]=avg_plus[i];
  }

  gavg=new TGraphAsymmErrors(16, avg_freq_forplot, avg_atten_forplot,0,0,avg_minus,avg_plus);
  
  gavgchristian=new TGraphAsymmErrors(nbands, freq_forplot[1], len_atten_avg[1],0,0,0,0);//plots for paper
  //   len_error_minus_avg[1],len_error_plus_avg[1]);
  gavgbaby=new TGraphAsymmErrors(nbands, freq_forplot[2], len_atten_avg[2],0,0,0,0);
  //	    len_error_minus_avg[2],len_error_plus_avg[2]);
  gavglong=new TGraphAsymmErrors(nbandslong, freq_forplot[0], len_atten_avg[0],0,0,0,0);
				 //len_error_minus_avg[0],len_error_plus_avg[0]);
  
  TF1 *ffreqavg=new TF1("ffreqavg","1.272e6*pow(x,[0])",1e8,1e9);//define fit
  ffreqavg->SetParNames("Power Law");
  gavg->Fit(ffreqavg);

  
  TH2F *h3=new TH2F("Field Atten. length (m) vs. freq., christian's", 
		    "Field Atten. length (m) vs. freq., christian's",
		    100,.95e8,1.05e9,100,10,300);
  //	 100,0,1.05e9,100,0,500);
  
  h3->SetXTitle("Frequency (Hz)");
  h3->SetYTitle("Field Attenuation Length (m)");
  
  if (logscale!=1) leg3= new TLegend(0.6,0.5,.80,0.90);
  if (logscale==1) leg3= new TLegend(0.18,0.18,.42,0.52);
  if (logscale!=1) legavg= new TLegend(0.6,0.5,0.80,0.90);
  if (logscale==1) legavg= new TLegend(0.2,0.2,0.70,0.30);
  if (includeshallow==1) leg3->AddEntry(g3,"Christian's, 10'","P");
  if (includeshallow==1) leg3->AddEntry(g3b,"Christian's, 20'","P");
  leg3->AddEntry(g3c,"Christian's, 30'","P");
  leg3->AddEntry(g3d,"Christian's, 40'","P");
  leg3->AddEntry(g3e,"Christian's, 50'","P");
  leg3->AddEntry(g3f,"Christian's, 60'","P");
  leg3->AddEntry(g3g,"Christian's, 70'","P");
  leg3->AddEntry(g3h,"Christian's, 80'","P");
  leg3->AddEntry(g3i,"Christian's, 90'","P");
  leg3->AddEntry(gline,"1/f, arbitrary scaling","l");
  legavg->AddEntry(gavgchristian,"Average Attenuation Length (Over All Depths)","P");
  //leg4->AddEntry(gline,"1/f, arbitrary scaling","l");
  
  
  TH2F *hnoise=new TH2F("(withnoise-withoutnoise)/withnoise",
			"(withnoise-withoutnoise)/withnoise",
			100,0,300,100,-20,20);
  hnoise->SetXTitle("Attenuation length (m)");
  hnoise->SetYTitle("% difference");
  for (int i=0;i<nantennas;i++){
    for (int ctr=0;ctr<nbands;ctr++){
      for (int j=1;j<ndepths;j++){
	if (j!=3)
	  hnoise->Fill(len_atten[i][j][ctr],len_atten_noise[i][j][ctr]);
      }
    }
  }

  //draw the stuff!
  c1->cd(1);//depth
  h1->Draw();
  g1->SetMarkerSize(.75);
  g1->SetLineColor(1);//vpk2pk=black
  g1->Draw("lP");
  g1b->SetMarkerSize(.75);
  g1b->SetLineColor(2);//fft bin0=red
  g1b->SetMarkerColor(2);//fft bin0=red
  g1b->Draw("lP");
  g1c->SetMarkerSize(.75);
  g1c->SetLineColor(6);//fft bin 1=magenta
  g1c->SetMarkerColor(6);//fft bin 1=magenta
  g1c->Draw("lP");
  g1d->SetMarkerSize(.75);
  g1d->SetLineColor(5);//bin2=yellow
  g1d->SetMarkerColor(5);//bin2=yellow
  g1d->Draw("lP");
  g1e->SetMarkerSize(.75);
  g1e->SetLineColor(3);//bin 3=green
  g1e->SetMarkerColor(3);//bin 3=green
  g1e->Draw("lP");
  g1f->SetMarkerSize(.75);
  g1f->SetLineColor(7);//bin4=cyan
  g1f->SetMarkerColor(7);//bin4=cyan
  g1f->Draw("lP");
  g1g->SetMarkerSize(.75);
  g1g->SetLineColor(4);//bin5=blue
  g1g->SetMarkerColor(4);//bin5=blue
  g1g->Draw("lP");
  leg1->Draw();
  
  //now for the paper plot
  cpaper1->cd(0);
  hdepth->Draw();
  gdepth->SetMarkerSize(.75);
  gdepth->SetLineColor(1);//250 MHz=bin1
  legdepth->Draw();
  gdepth->Draw("P");
  
  c1->cd(2);//frequency
  if (logscale==1) gPad->SetLogy();
  if (logscale==1) gPad->SetLogx();
  h2->Draw();
  g2->SetMarkerSize(.8);
  g2->SetMarkerColor(2);//long shallow=red circle
  if (includeshallow==1) g2->Draw("P");
  g2b->SetMarkerSize(.8);
  g2b->SetMarkerColor(2);//long mid=red triangle
  g2b->SetMarkerStyle(22);
  g2b->Draw("P");
  g2c->SetMarkerSize(.8);
  g2c->SetMarkerColor(2);//long deep=red square
  g2c->SetMarkerStyle(21);
  g2c->Draw("P");
  g2d->SetMarkerSize(.8);
  g2d->SetMarkerColor(3);//christ shallow=green circle
  if (includeshallow==1) g2d->Draw("P");
  g2e->SetMarkerSize(.8);
  g2e->SetMarkerColor(3);//christ mid=green triangle
  g2e->SetMarkerStyle(22);
  g2e->Draw("P");
  g2f->SetMarkerSize(.8);
  g2f->SetMarkerColor(3);//christ deep=green square
  g2f->SetMarkerStyle(21);
  g2f->Draw("P");
  g2g->SetMarkerSize(.8);
  g2g->SetMarkerColor(4);//baby shallow=blue circle
  if (includeshallow==1) g2g->Draw("P");
  g2h->SetMarkerSize(.8);
  g2h->SetMarkerColor(4);//baby mid=blue triangle
  g2h->SetMarkerStyle(22);
  g2h->Draw("P");
  g2i->SetMarkerSize(.8);
  g2i->SetMarkerColor(4);//baby deep=blue square
  g2i->SetMarkerStyle(21);
  g2i->Draw("P");
  gline->Draw();
  leg2->Draw();

  c1->Print("attenuation.eps");
  
  //now redraw for paper plots
  cpaper2->cd(0);//frequency
  if (logscale==1) gPad->SetLogy();
  if (logscale==1) gPad->SetLogx();
  h2->Draw();
  gfreq->SetMarkerSize(0);
  gStyle->SetOptFit(0000);//don't print chisquare
  gfreq->Draw("P");
  g2->SetMarkerSize(.8);
  g2->SetMarkerColor(2);//long shallow=red circle
  if (includeshallow==1) g2->Draw("P");
  g2b->SetMarkerSize(.8);
  g2b->SetMarkerColor(2);//long mid=red triangle
  g2b->SetMarkerStyle(22);
  g2b->Draw("P");
  g2c->SetMarkerSize(.8);
  g2c->SetMarkerColor(2);//long deep=red square
  g2c->SetMarkerStyle(21);
  g2c->Draw("P");
  g2d->SetMarkerSize(.8);
  g2d->SetMarkerColor(3);//christ shallow=green circle
  if (includeshallow==1) g2d->Draw("P");
  g2e->SetMarkerSize(.8);
  g2e->SetMarkerColor(3);//christ mid=green triangle
  g2e->SetMarkerStyle(22);
  g2e->Draw("P");
  g2f->SetMarkerSize(.8);
  g2f->SetMarkerColor(3);//christ deep=green square
  g2f->SetMarkerStyle(21);
  g2f->Draw("P");
  g2g->SetMarkerSize(.8);
  g2g->SetMarkerColor(4);//baby shallow=blue circle
  if (includeshallow==1) g2g->Draw("P");
  g2h->SetMarkerSize(.8);
  g2h->SetMarkerColor(4);//baby mid=blue triangle
  g2h->SetMarkerStyle(22);
  g2h->Draw("P");
  g2i->SetMarkerSize(.8);
  g2i->SetMarkerColor(4);//baby deep=blue square
  g2i->SetMarkerStyle(21);
  g2i->Draw("P");
  //gline->Draw();
  leg2->Draw();

  cpaper3->cd(0);
  if (logscale==1) gPad->SetLogy();
  if (logscale==1) gPad->SetLogx();
  h2->Draw();
  gavg->SetMarkerSize(0.8);
  legavg->Draw();
  gStyle->SetOptFit(0000);//don't print chisquare
  gavg->Draw("P");
  //gavgchristian->SetMarkerSize(.8);
  //gavgchristian->Draw("P");
  //gavgbaby->SetMarkerSize(.8);
  //gavgbaby->SetMarkerColor(2);
  //gavgbaby->Draw("P");
  //gavglong->SetMarkerSize(.8);
  //gavglong->SetMarkerColor(3);
  //gavglong->Draw("P");
  
  c2->cd(1);//christian's frequency
  if (logscale==1) gPad->SetLogy();
  if (logscale==1) gPad->SetLogx();
  h3->Draw();
  g3->SetMarkerSize(.75);
  g3->SetLineColor(1);//10'=black
  if (includeshallow==1) g3->Draw("lP");
  g3b->SetMarkerSize(.75);
  g3b->SetLineColor(2);//20'=red
  g3b->SetMarkerColor(2);//20'=red
  if (includeshallow==1) g3b->Draw("lP");
  g3c->SetMarkerSize(.75);
  g3c->SetLineColor(6);//30'=magenta
  g3c->SetMarkerColor(6);//30'=magenta
  g3c->Draw("lP");
  g3d->SetMarkerSize(.75);
  g3d->SetLineColor(5);//40'=yellow
  g3d->SetMarkerColor(5);//40'=yellow
  g3d->Draw("lP");
  g3e->SetMarkerSize(.75);
  g3e->SetLineColor(3);//50'=green
  g3e->SetMarkerColor(3);//50'=green
  g3e->Draw("lP");
  g3f->SetMarkerSize(.75);
  g3f->SetLineColor(7);//60'=cyan
  g3f->SetMarkerColor(7);//60'=cyan
  g3f->Draw("lP");
  g3g->SetMarkerSize(.75);
  g3g->SetLineColor(4);//70'=blue
  g3g->SetMarkerColor(4);//70'=blue
  g3g->Draw("lP");
  g3h->SetMarkerSize(.75);
  g3h->SetLineColor(40);//80'= grey
  g3h->SetMarkerColor(40);//80'= grey
  g3h->Draw("lP");
  g3i->SetMarkerSize(.75);
  g3i->SetLineColor(50);//90'= brick
  g3i->SetMarkerColor(50);//90'= brick
  g3i->Draw("lP");
  gline->Draw();
  leg3->Draw();
  c2->cd(2);
  hnoise->SetMarkerSize(.75);
  hnoise->Draw();
  c2->Print("attenuation_2.eps");
  cpaper1->Print("plots_for_paper1.eps");
  cpaper2->Print("plots_for_paper2.eps");
  cpaper3->Print("plots_for_paper3.eps");

    }

float calcAtten(float power1, float power2, float dist1, float dist2){
  float len_atten=(1./((log(sqrt(power1)*dist1)-log(sqrt(power2)*dist2))/(dist2-dist1)));
  return len_atten;
}

TStyle* RootStyle() {

  //  const char* modified = "Borrowed and adapted from paus et al";

  TStyle *RootStyle = new TStyle("Root-Style","The Perfect Style for Plots ;-)");

#ifdef __CINT__
  TStyle *GloStyle;
  GloStyle = gStyle;                          // save the global style reference

  gStyle = RootStyle;
#endif
// otherwise you need to call TROOT::SetStyle("Root-Style")

  // Paper size

  RootStyle->SetPaperSize(TStyle::kUSLetter);

  // Canvas

  RootStyle->SetCanvasColor     (0);
  RootStyle->SetCanvasBorderSize(10);
  RootStyle->SetCanvasBorderMode(0);
  RootStyle->SetCanvasDefH      (600);
  RootStyle->SetCanvasDefW      (600);
  RootStyle->SetCanvasDefX      (10);
  RootStyle->SetCanvasDefY      (10);

  // Pads

  RootStyle->SetPadColor       (0);
  RootStyle->SetPadBorderSize  (10);
  RootStyle->SetPadBorderMode  (0);
  //  RootStyle->SetPadBottomMargin(0.13);
  RootStyle->SetPadBottomMargin(0.16);
  RootStyle->SetPadTopMargin   (0.08);
  RootStyle->SetPadLeftMargin  (0.15);
  RootStyle->SetPadRightMargin (0.15);
  RootStyle->SetPadGridX       (0);
  RootStyle->SetPadGridY       (0);
  RootStyle->SetPadTickX       (1);
  RootStyle->SetPadTickY       (1);

  // Frames

  RootStyle->SetFrameFillStyle ( 0);
  RootStyle->SetFrameFillColor ( 0);
  RootStyle->SetFrameLineColor ( 1);
  RootStyle->SetFrameLineStyle ( 0);
  RootStyle->SetFrameLineWidth ( 2);
  RootStyle->SetFrameBorderSize(10);
  RootStyle->SetFrameBorderMode( 0);


  // Histograms

  RootStyle->SetHistFillColor(0);
  RootStyle->SetHistFillStyle(1);
  RootStyle->SetHistLineColor(1);
  RootStyle->SetHistLineStyle(0);
  RootStyle->SetHistLineWidth(2);

  // Functions

  RootStyle->SetFuncColor(1);
  RootStyle->SetFuncStyle(0);
  RootStyle->SetFuncWidth(1);

  //Legends 

  RootStyle->SetStatBorderSize(2);
  RootStyle->SetStatFont      (42);
  //  RootStyle->SetOptStat       (111111);
  RootStyle->SetOptStat       (0);
  RootStyle->SetStatColor     (0);
  RootStyle->SetStatX         (0.93);
  RootStyle->SetStatY         (0.90);
  RootStyle->SetStatFontSize  (0.07);
  //  RootStyle->SetStatW         (0.2);
  //  RootStyle->SetStatH         (0.15);

  // Labels, Ticks, and Titles

  RootStyle->SetTickLength ( 0.015,"X");
  RootStyle->SetTitleSize  ( 0.055,"X");
  RootStyle->SetTitleOffset( 1.00,"X");
  RootStyle->SetTitleBorderSize(0);
  //  RootStyle->SetTitleFontSize((float)3.);
  RootStyle->SetLabelOffset( 0.015,"X");
  RootStyle->SetLabelSize  ( 0.050,"X");
  RootStyle->SetLabelFont  ( 42   ,"X");

  RootStyle->SetTickLength ( 0.015,"Y");
  RootStyle->SetTitleSize  ( 0.055,"Y");
  RootStyle->SetTitleOffset( 1.300,"Y");
  RootStyle->SetLabelOffset( 0.015,"Y");
  RootStyle->SetLabelSize  ( 0.050,"Y");
  RootStyle->SetLabelFont  ( 42   ,"Y");

  RootStyle->SetTitleFont  (42,"XYZ");
  RootStyle->SetTitleColor  (1);

  // Options

  RootStyle->SetOptFit     (0001);

  RootStyle->SetMarkerStyle(20);
  RootStyle->SetMarkerSize(1.35);

  //  cout << ">> Style initialized with the Root Style!" << endl;
  //  cout << ">> " << modified << endl << endl;
  return RootStyle;
}
