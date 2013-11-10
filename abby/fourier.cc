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

using namespace std;

const double DEGRAD=57.2957795;    // degree/rad
const double MFT=0.3048;

const int HOWMANYDISKS=8;// how many days of measurements
const int MAXMEASUREMENTS=39; // how many files we're reading per day (how many mesurements)
const int MAXTIMESTEPS=15000; // record length

// which disk and file number to plot
int WHICHDISK=2;
int WHICHFILENUMBER=31;

ifstream infile; // this is the input file, different for every measurement
ofstream fout("out.txt"); // outputs average value of fourier transform for each frequency bin for each measurement
string whichfilename; // different for every measurement
string directory="/home/rjn/saltStuff/salt_trip3/data/";

//This does the work.  From numerical recipes.
// data is the array of data
// isign=1 from time->frequency, isign=-1 for frequency->time
// nsize is the size of the array
void four1(double *data, const int isign,int nsize);
void four1floats(float *data, const int isign,int nsize); // same, but acts of floats instead of doubles
float AverageVoltage(float *data,int ilower,int iupper); // averages the values in the array between ilower and iupper
float AveragePower(float *data,int ilower,int iupper); // averages the values in the array between ilower and iupper
float SumVoltage(float *data,int ilower,int iupper); // sums the values in the array between ilower and iupper
float SumPower(float *data,int ilower,int iupper); // sums the square of the values in the array between ilower and iupper



inline void SWAP(double &a, double &b) // swaps two numbers
	{double dum=a; a=b; b=dum;}
inline void SWAPfloats(float &a, float &b)
	{float dum=a; a=b; b=dum;}


int noise_switch=0;
//512=1/2.E-10/Deltaf
const int MAX_ARRAY_SIZE=16384*2; // the array size need to be a power of 2 and
// 500 is pretty close to 512, so use that
// we need another factor of 2 for the imaginary components

float voltage_timedomain[MAX_ARRAY_SIZE]; // voltages in time domain voltages
float voltage_firstpulse[MAX_ARRAY_SIZE]; //voltage cut in on first pulse only
float voltage_frequencydomain[MAX_ARRAY_SIZE]; // voltages in frequency domain from fourier transform + we add other factors

float vtime_forplot[MAX_ARRAY_SIZE/2]; // this is the same as voltage_timedomain with the 0's removed
float vfreq_forplot[MAX_ARRAY_SIZE/2]; // this is for plotting, contains sqrt of sum of squares of real and imaginary components for each step
float power_forplot[MAX_ARRAY_SIZE/2];
float pfreq_forplot[MAX_ARRAY_SIZE/2]; // powers
float xcount[MAX_ARRAY_SIZE/2]; // this is just the bin numbers from 0 to N
float t[MAX_ARRAY_SIZE/2]; // this is the time for each bin, in s
float f[MAX_ARRAY_SIZE/2]; // this is the freq for each bin, in Hz
float integral_timedomain=0.;
float integral_frequencydomain=0.;
float factor;  //accounts for attenuation, etc.
float TOTALENERGY=5.61E-5; //total energy in PCD pulse in Joules

const float LOWEREDGE1=50.E6; // lower edge of first frequency bin
const float UPPEREDGE1=150.E6; // upper edge of first freqency bin
const float LOWEREDGE2=150.E6; // lower edge 2nd frequency bin
const float UPPEREDGE2=250.E6; // upper edge 2nd freq bin
const float LOWEREDGE3=250.E6;
const float UPPEREDGE3=350.E6;

//int iloweredge1,iloweredge2,iloweredge3; // bin number that corresponds to lower edge of each frequency band
//int iupperedge1,iupperedge2,iupperedge3; // bin number that corresponds to upper edge of each frequency band

char ctest[3]; // temperary character string


// data from trip 1
float distances_trip1[3]={91.5*MFT,189.5*MFT,314.5*MFT};
float signal_300MHz[3]={1.,1./27.5,1./100.};
float signal_150MHz[3]={1.,1./17.0,1./64.6};

float S11POWERTRANS_TRIP1=0.01*0.899*0.858; // power transmitted by antenna in 250-350 MHz frequency band on first trip
float S11POWERTRANS_TRIP2=0.145*0.50*0.50; // power transmitted by antenna in 250-350 MHz frequency band on second trip

TStyle* RootStyle();
TStyle *color=RootStyle();

Float_t z[5],x[5],y[5],errorz[5];

Double_t func(Double_t *x,Double_t *par)
{

  float x_solid;
  float x_fractured;
  float atten;
 
  float r; // total distance traveled in salt

  // relevant distances in ft.
  float distance_fromclosewall=27.;
  float distance_straightacross=553.;
  float depth_of_cubby=12.;
  float distance_cubby_from_wall=20.;
  float width_cubby=21.;

  // d=par[0]
  // asolid=par[1]
  // afractured=par[2]

  // x[0]=theta
  // x[1]=r


  // find r
  if (x[0]<-1.*atan(distance_fromclosewall/distance_straightacross)) { // along the first wall
    r=distance_fromclosewall/sin(fabs(x[0]));
  }
  else if (x[0]>-1.*atan(distance_fromclosewall/distance_straightacross) &&
	   x[0]<-1.*atan((distance_fromclosewall-distance_cubby_from_wall)/distance_straightacross)) { // between wall and edge of cubby
    r=distance_straightacross/cos(x[0]);
  }
  else if (x[0]>-1.*atan((distance_fromclosewall-distance_cubby_from_wall)/distance_straightacross) && 
	   x[0]<atan((distance_cubby_from_wall+width_cubby-distance_fromclosewall)/distance_straightacross)) { // between edge of cubby and straight across
    r=(distance_straightacross-depth_of_cubby)/cos(x[0]);
  }
  else if (x[0]>atan((distance_cubby_from_wall+width_cubby-distance_fromclosewall)/distance_straightacross)) {// beyond edge of cubby
    r=distance_straightacross/cos(x[0]);
  }
  

  if (x[0]>-1.*atan(distance_fromclosewall/distance_straightacross))
    x_fractured=2.*par[0]/cos(x[0]);
  else {
    if(par[0]<27. && x[0]>-1.*atan((distance_fromclosewall-par[0])/par[0]))
      x_fractured=-1.*par[0]/sin(x[0])+par[0]/cos(x[0]);
    else 
      x_fractured=r;
  }
  //  cout << "borderline is " << -1.*atan((distance_fromclosewall-par[0])/par[0]) << "\n";

  x_solid=r-x_fractured;

  if (x_solid<0)
    cout << "x[0], r, x_fractured are " << x[0] << " " << r << " " << x_fractured << "\n";

  atten= log(par[3]*exp(-1.*x_fractured/par[2])*exp(-1.*x_solid/par[1]));
  
 
  //  cout << "atten is " << x_fractured << " " << par[2] << " " << x_solid << " " << par[1] << " " << atten << "\n";

 return atten;
}


//  void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
//  {
//    // this is the chisquare

//     const Int_t nbins = 5;
//     Int_t i;

//  //calculate chisquare
//     Double_t chisq = 0;
//     Double_t delta;
//     for (i=0;i<nbins; i++) {
//       delta  = (z[i]-func(x[i],y[i],par))/errorz[i];
//       chisq += delta*delta;
//     }
//     f = chisq;
//  }


int main() {

  gStyle=color;

  int howmanymeasurements[HOWMANYDISKS]={20,9,39,26,14,0,0,0}; // these are just time domain measurements

  int whichfiles[HOWMANYDISKS][MAXMEASUREMENTS]={
    {0,1,2,3,4,5,8,10,14,16,18,20,22,24,26,28,31,33,35,36},
    {0,1,2,3,4,5,6,7,8},
    {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,
     32,33,34,35,36,37,38},
    {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25},
    {2,3,4,9,10,11,12,13,14,15,16,17,18,19},
    {},
    {},
    {}};
  
  
  int howmanytimesteps[HOWMANYDISKS][MAXMEASUREMENTS];
  
  TGraph *g1;
  TGraph *g2;
  TH2F *htime;
  TH2F *hfreq;
  
  int idisk, ifile; //ith file on the ith disk.
  int filenumber=0; //number that is in the title of the filename
  int counttimesteps;
  int ifreqsteps=MAX_ARRAY_SIZE/2;
  float maxvoltage,maxvoltage_freq, maxpower_freq;
  float v_max,v_min,vrms,vsumsquare;//for peak to peak and  vrms
  float t_maxvoltage;//time of peak voltage


  TFile *waveformfile=new TFile("waveforms.root","RECREATE");
  TTree *waveformtree = new TTree("waveformtree","waveformtree"); // waveform tree filled for each file that we recorded
  
  waveformtree->Branch("idisk",&idisk,"idisk/I");
  waveformtree->Branch("ifile",&ifile,"ifile/I");
  waveformtree->Branch("itimebins",&counttimesteps,"counttimesteps/I");     
  waveformtree->Branch("ifreqbins",&ifreqsteps,"ifreqsteps/I");
  waveformtree->Branch("f",&f,"f[ifreqsteps]/F");
  waveformtree->Branch("t",&t,"t[ifreqsteps]/F");
  waveformtree->Branch("vtime",&vtime_forplot,"vtime_forplot[ifreqsteps]/F");
  waveformtree->Branch("vfreq",&vfreq_forplot,"vfreq_forplot[ifreqsteps]/F");
  waveformtree->Branch("maxvoltage",&maxvoltage,"maxvoltage/F");
  waveformtree->Branch("maxvoltage_freq",&maxvoltage_freq,"maxvoltage_freq/F");
  waveformtree->Branch("maxpower_freq",&maxpower_freq,"maxpower_freq/F");
  waveformtree->Branch("power",&power_forplot,"power_forplot[ifreqsteps]/F");
  waveformtree->Branch("filenumber",&filenumber,"filenumber/I");
  waveformtree->Branch("v_max",&v_max,"v_max/F");
  waveformtree->Branch("v_min",&v_min,"v_min/F");      
  waveformtree->Branch("vrms",&vrms,"vrms/F");
  
  for (idisk=0;idisk<HOWMANYDISKS;idisk++) {
    for (ifile=0;ifile<howmanymeasurements[idisk];ifile++) {
      
      integral_timedomain=0.;
      integral_frequencydomain=0.;
      
      
      
      for (int j=0;j<MAX_ARRAY_SIZE;j++) { // initialize voltage vs. time
	voltage_timedomain[j]=0.;
	t[j]=0;
      }
      
      whichfilename="/home/rjn/saltStuff/salt_trip3_data/";
      if (idisk==0)
	whichfilename+="above_ground_data/"; // all the file names start with this string
      if (idisk==1)
	whichfilename+="first_day/"; // all the file names start with this string
      if (idisk==2)
	whichfilename+="disk_3/"; // all the file names start with this string
      if (idisk==3)
	whichfilename+="disk_4/"; // all the file names start with this string
      if (idisk==4)
	whichfilename+="disk_5/"; // all the file names start with this string
      if (idisk==5)
	whichfilename+="disk_6/"; // all the file names start with this string
      if (idisk==6)
	whichfilename+="disk_7/"; // all the file names start with this string
      if (idisk==7)
	whichfilename+="disk_8/"; // all the file names start with this string
      
      
      whichfilename+="TEK000";

      if (whichfiles[idisk][ifile]<10)
	whichfilename+="0";

      sprintf(ctest,"%i",whichfiles[idisk][ifile]); // write number of measurement into character string ctest
	
      whichfilename+=(string)ctest; // append it to the file  name
	
      whichfilename+=".CSV"; // add extension
      


      //      cout << "input file is " << whichfilename << "\n";
      
      ifstream is(whichfilename.c_str()); // this is the input file 
      if(!is) {
	cerr << "Bollocks" << endl;
	exit(0);
      }

      char cjunk1[40]; // temporary character strings
      char cjunk2[40];
      
      //       for (int j=0;j<MAXTIMESTEPS;j++) {
      
      counttimesteps=0;
      maxvoltage=0.; // max voltage in the time domain
      v_max=0;
      v_min=0;
      vrms=0;
      vsumsquare=0;
      while (!is.eof()) {
	
 	is.getline(cjunk1,40,','); // go 40 characters at most, until you see a comma
	// 	// fin.getline(cjunk2,40,',');
 	is.getline(cjunk2,40,'\n'); // go another 40 characters until you see a return
	
 	voltage_timedomain[2*counttimesteps]=atof(cjunk2); // changes cjunk2 to a float
 	voltage_timedomain[2*counttimesteps+1]=0.; // fill in the imaginary components with 0's
	
	if (fabs(voltage_timedomain[2*counttimesteps])>maxvoltage)
	  maxvoltage=fabs(voltage_timedomain[2*counttimesteps]); // max voltage in the time domain

 	t[counttimesteps]=atof(cjunk1); // changes cjunk1 to float, it is the time
	counttimesteps++;
      } // end loop through time steps
      counttimesteps--; // the way we did it, it overcounts by 1, so subtract 1
      
      for (int j=counttimesteps*2;j<MAX_ARRAY_SIZE;j++) {
       	voltage_timedomain[j]=0.; // fills in the rest of the array with 0's
      }
      
      for (int j=counttimesteps;j<MAX_ARRAY_SIZE/2;j++) {
       	t[j]=t[j-1]+(t[j-1]-t[j-2]); // fills in times appropriately
      }
      
      for (int j=0; j<MAX_ARRAY_SIZE;j++) voltage_firstpulse[j]=voltage_timedomain[j];
      
      //abby's hack for cutting timewindow to the first pulse
      if (idisk==2 && whichfiles[idisk][ifile]==10){//50' christian antennas T12
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>0.035e-6 || t[j]<-0.010e-6 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==3 && whichfiles[idisk][ifile]==6){//50' christian antennas T13
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>0.035e-6 || t[j]<-0.010e-6 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==2 && whichfiles[idisk][ifile]==31){//50' longpole antennas T12
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>0.040e-6 ||t[j]<-0.020e-6 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==3 && whichfiles[idisk][ifile]==23){//50' longpole antennas T13
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>0.035e-6 || t[j]<-0.020e-6 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==2 && whichfiles[idisk][ifile]==26){//50' baby antennas T12
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>10.e-9 ||t[j]<-2.e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==3 && whichfiles[idisk][ifile]==19){//50' baby antennas T13
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>20.e-9 || t[j]<0.e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==2 && whichfiles[idisk][ifile]==19){//90' christian antennas T12
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>0.04e-6 || t[j]<-0.03e-6 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==3 && whichfiles[idisk][ifile]==14){//90' christian antennas T13
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>0.04e-6 || t[j]<-0.03e-6|| t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==2 && whichfiles[idisk][ifile]==32){//75' longpole antennas T12
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>0.035e-6 || t[j]<-0.02e-6 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==3 && whichfiles[idisk][ifile]==24){//75' longpole antennas T13
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>0.03e-6 || t[j]<-0.01e-6 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==3 && whichfiles[idisk][ifile]==25){//75' longpole antennas T13, test for Amy
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>910.e-9 || t[j]<870.e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==2 && whichfiles[idisk][ifile]==28){//90' baby antennas T12
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	   if(t[j]>15.e-9 || t[j]<-10e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==3 && whichfiles[idisk][ifile]==20){//90' baby antennas T13
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>8.e-9 || t[j]<-3e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==2 && whichfiles[idisk][ifile]==4){//20' christian antennas T12
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>0.04e-6 || t[j]<-0.005e-6 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==3 && whichfiles[idisk][ifile]==1){//20' christian antennas T13
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>0.045e-6 || t[j]<-0.005e-6 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==2 && whichfiles[idisk][ifile]==33){//20' longpole antennas T12
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>0.055e-6 || t[j]<-0.01e-6 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==3 && whichfiles[idisk][ifile]==22){//20' longpole antennas T13
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>50.e-9 || t[j]<10e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==2 && whichfiles[idisk][ifile]==25){//10' baby antennas T12
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>10.e-9 || t[j]<-5e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==3 && whichfiles[idisk][ifile]==18){//10' baby antennas T13
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>0.027e-6 || t[j]<0.e-6 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      
      //other heights of christian antennas.
      if (idisk==2 && whichfiles[idisk][ifile]==2){//10' christian antennas T12
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>40.e-9 || t[j]<-10.e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==3 && whichfiles[idisk][ifile]==0){//10' christian antennas T13
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>40.e-9 || t[j]<-5.e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==2 && whichfiles[idisk][ifile]==6){//30' christian antennas T12
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>25.e-9 || t[j]<-10.e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==3 && whichfiles[idisk][ifile]==3){//30' christian antennas T13
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>-40.e-9 || t[j]<-70.e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==2 && whichfiles[idisk][ifile]==9){//40' christian antennas T12
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>40e-9 || t[j]<-10e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==3 && whichfiles[idisk][ifile]==4){//40' christian antennas T13
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>45e-9 || t[j]<0e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==2 && whichfiles[idisk][ifile]==11){//60' christian antennas T12
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>45e-9 || t[j]<-10e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==3 && whichfiles[idisk][ifile]==8){//60' christian antennas T13
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>12e-9 || t[j]<-25e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==2 && whichfiles[idisk][ifile]==14){//70' christian antennas T12
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>30e-9 || t[j]<-25e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==3 && whichfiles[idisk][ifile]==10){//70' christian antennas T13
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>20e-9 || t[j]<-25e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==2 && whichfiles[idisk][ifile]==16){//80' christian antennas T12
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>40e-9 || t[j]<-30e-9|| t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      if (idisk==3 && whichfiles[idisk][ifile]==12){//80' christian antennas T13
	for (int j=0; j<MAX_ARRAY_SIZE;j++){
	  if(t[j]>20e-9 || t[j]<-25e-9 || t[j]==0) voltage_firstpulse[2*j]=0.; 
	}
      }
      //end abby's timewindow cuts

      for (int j=0; j<MAX_ARRAY_SIZE;j++){
	if (voltage_firstpulse[2*j]>v_max){
	  v_max=voltage_firstpulse[2*j];
	  if (fabs(v_max)>fabs(v_min)) t_maxvoltage=t[j];
	}
	if (voltage_firstpulse[2*j]<v_min){
	  v_min=voltage_firstpulse[2*j];
	  if (fabs(v_min)>fabs(v_max)) t_maxvoltage=t[j];
	}
	vsumsquare+=voltage_firstpulse[2*j];
      }
      vrms=sqrt(vsumsquare);
      
      //hack for making all the time windows the same for all antennas.
      if (noise_switch!=1){
	if ((idisk==2 && whichfiles[idisk][ifile]<=24) || 
	    (idisk==3 && whichfiles[idisk][ifile]<=17)){ //christian antennas
	  //	  cout << idisk << "\t" << whichfiles[idisk][ifile]  << "\t" << t_maxvoltage << endl;
	  for (int j=0; j<MAX_ARRAY_SIZE;j++){
	    if(t[j]<-10e-9+t_maxvoltage || t[j]>27e-9+t_maxvoltage || t[j]==0) 
	      //RJN change
	      voltage_timedomain[j]=0.; //cut to 10 ns before and 27 ns after peak.
	  }


	  int firstVal=0;
	  int lastVal=0;
	  for (int j=0; j<MAX_ARRAY_SIZE;j++){
	    if(voltage_timedomain[j]!=0) {
	      if(!firstVal)
		firstVal=j;
	      lastVal=j;
	    }
	  }
	  //	  cout << idisk << "\t" << whichfiles[idisk][ifile] << "\t" 
	  //	       << t[firstVal] << "\t" << t[lastVal] << "\t"
	  //	       << t[lastVal]-t[firstVal] << "\n";

	}
	if ((idisk==2 && whichfiles[idisk][ifile]>=30 && whichfiles[idisk][ifile]<=33) ||
	    (idisk==3 && whichfiles[idisk][ifile]>=22 && whichfiles[idisk][ifile]<=24)){//longpoles
	  cout << "Long\t" << idisk << "\t" << whichfiles[idisk][ifile]  << "\t" << t_maxvoltage << endl;
	  for (int j=0; j<MAX_ARRAY_SIZE;j++){
	    if(t[j]<-15e-9+t_maxvoltage || t[j]>23e-9+t_maxvoltage || t[j]==0) 
	    voltage_timedomain[j]=0.; //cut to 15 ns before and 23 ns after peak.
	  }
	}
	if ((idisk==2 && whichfiles[idisk][ifile]>=25 && whichfiles[idisk][ifile]<=29) ||
	    (idisk==3 && whichfiles[idisk][ifile]>=18 && whichfiles[idisk][ifile]<=20)){//babies
	  cout << "Baby\t" << idisk << "\t" << whichfiles[idisk][ifile]  << "\t" << t_maxvoltage << endl;
	  for (int j=0; j<MAX_ARRAY_SIZE;j++){
	    if(t[j]<-5e-9+t_maxvoltage || t[j]>7e-9+t_maxvoltage || t[j]==0) 
	      voltage_timedomain[j]=0.; //cut to 5 ns before and 7 ns after peak.
	  }
	}
      }
      if (noise_switch==1){
	if ((idisk==2 && whichfiles[idisk][ifile]<=24) || 
	    (idisk==3 && whichfiles[idisk][ifile]<=17)){ //christian antennas
	  for (int j=0; j<MAX_ARRAY_SIZE;j++){
	    if(t[j]>t[0]+37e-9 || t[j]==0) 
	      voltage_timedomain[j]=0.; //cut to 37 ns at beginning of the trace
	  }
	}
	if ((idisk==2 && whichfiles[idisk][ifile]>=30 && whichfiles[idisk][ifile]<=33) ||
	    (idisk==3 && whichfiles[idisk][ifile]>=22 && whichfiles[idisk][ifile]<=24)){//longpoles
	  for (int j=0; j<MAX_ARRAY_SIZE;j++){
	    if(t[j]>38e-9+t[0] || t[j]==0) 
	      voltage_timedomain[j]=0.; //cut to 38 ns at beginning of the trace
	  }
	}
	if ((idisk==2 && whichfiles[idisk][ifile]>=25 && whichfiles[idisk][ifile]<=29) ||
	    (idisk==3 && whichfiles[idisk][ifile]>=18 && whichfiles[idisk][ifile]<=20)){//babies
	  for (int j=0; j<MAX_ARRAY_SIZE;j++){
	    if(t[j]>12e-9+t[0] || t[j]==0) 
	      voltage_timedomain[j]=0.; //cut to 12 ns at beginning of the trace
	  }
	}
      }

      // N=1/(Delta f * Delta t)
      f[0]=0.;
      double fstep=1./((float)MAX_ARRAY_SIZE/2*(t[1]-t[0])); // f[1]=Delta f
      

      for (int j=1;j<MAX_ARRAY_SIZE/4;j++) { // fill in the rest of the frequencies appropriately
 	f[j]=f[j-1]+fstep;

      }
      f[MAX_ARRAY_SIZE/4]=-1.*f[MAX_ARRAY_SIZE/4-1];

      for (int j=MAX_ARRAY_SIZE/4+1;j<MAX_ARRAY_SIZE/2;j++) {
	f[j]=f[j-1]+fstep;
      }
      
   
      
      for (int j=0;j<MAX_ARRAY_SIZE;j++) { // copy over time domain signal to frequency domain signal before putting it into the fourier transform program
	voltage_frequencydomain[j]=voltage_timedomain[j]; // initialization of frequency domain voltages for inputting into fourier transform program
	integral_timedomain+=voltage_timedomain[j]*voltage_timedomain[j]*(t[1]-t[0]); // integral of power in time domain -> total energy of the pulse
	// (in units where 50 Ohms=1)
      }

      // take the fourier transform
      four1floats(voltage_frequencydomain,1,MAX_ARRAY_SIZE); // 1 says that we're going from time to frequency
      
      maxvoltage_freq=0; // maximum field in the frequency domain
      for (int j=0;j<MAX_ARRAY_SIZE/2;j++) {

 	voltage_frequencydomain[2*j]*=sqrt(t[1]-t[0]);//*sqrt(counttimesteps);  // scale them for different width of time bins (scale voltage by sqrt Hz).
 	voltage_frequencydomain[2*j+1]*=sqrt(t[1]-t[0]);//*sqrt(counttimesteps);//decided not to scale by Npoints in trace because we have zero-padded to the same N for all waveforms.
	
	vtime_forplot[j]=voltage_timedomain[2*j]; // this just picks out the real components for plotting
	vfreq_forplot[j]=sqrt(voltage_frequencydomain[2*j]*voltage_frequencydomain[2*j]+voltage_frequencydomain[2*j+1]*voltage_frequencydomain[2*j+1]); // sum the real and imaginary --the 2 accounts for the negative frequencies
	
	power_forplot[j]=(voltage_frequencydomain[2*j]*voltage_frequencydomain[2*j]+voltage_frequencydomain[2*j+1]*voltage_frequencydomain[2*j+1]);
	if (vfreq_forplot[j]>maxvoltage_freq){
	  maxvoltage_freq=vfreq_forplot[j]; 
	  maxpower_freq=power_forplot[j];
	}
      }
      
      filenumber=whichfiles[idisk][ifile];

      if (idisk==WHICHDISK && whichfiles[idisk][ifile]==WHICHFILENUMBER) { // picks out which measurement you want to plot

	htime=new TH2F("time","time",counttimesteps,t[0]*1.1,t[counttimesteps]*1.1,1000,-1.1*maxvoltage,1.1*maxvoltage);
	hfreq=new TH2F("freq","freq",MAX_ARRAY_SIZE/2,-1.1*f[MAX_ARRAY_SIZE/4-1],1.1*f[MAX_ARRAY_SIZE/4-1],1000,0,1.1*maxvoltage_freq);

 	g1=new TGraph(MAX_ARRAY_SIZE/2,t,vtime_forplot);
 	g2=new TGraph(MAX_ARRAY_SIZE/2,f,vfreq_forplot);
      }
  
      waveformtree->Fill();

    } // end loop over measurements


  }// end loop over days

  waveformfile->Write();




  TCanvas *c1 = new TCanvas("saltmcplots1","saltmcplots1",1000,880);
  
  htime->Draw();
  g1->Draw();
  
  c1->Print("time.eps");
  
  hfreq->Draw();
  g2->Draw();
  
  c1->Print("freq.eps");

 




}







// this takes Fourier transforms, from numerical recipes

void four1(double *data, const int isign,int nsize) {
	int n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;

	int nn=nsize/2;
	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j-1],data[i-1]);
			SWAP(data[j],data[i]);
		}
		m=nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j-1]-wi*data[j];
				tempi=wr*data[j]+wi*data[j-1];
				data[j-1]=data[i-1]-tempr;
				data[j]=data[i]-tempi;
				data[i-1] += tempr;
				data[i] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
void four1floats(float *data, const int isign,int nsize) {
	int n,mmax,m,j,istep,i;
	float wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;

	int nn=nsize/2;
	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAPfloats(data[j-1],data[i-1]);
			SWAPfloats(data[j],data[i]);
		}
		m=nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j-1]-wi*data[j];
				tempi=wr*data[j]+wi*data[j-1];
				data[j-1]=data[i-1]-tempr;
				data[j]=data[i]-tempi;
				data[i-1] += tempr;
				data[i] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
float AveragePower(float *data,int ilower,int iupper) {

  float ave=0.;
  for (int i=ilower;i<=iupper;i++) {
    ave+=data[i]*data[i];
  }
  ave/=(float)(iupper-ilower+1);

  return sqrt(ave);

}
float AverageVoltage(float *data,int ilower,int iupper) {

  float ave=0.;
  for (int i=ilower;i<=iupper;i++) {
    ave+=data[i];
  }
  ave/=(float)(iupper-ilower+1);

  return ave;

}
float SumVoltage(float *data,int ilower,int iupper) {

  float ave=0.;
  for (int i=ilower;i<=iupper;i++) {
    ave+=data[i];
  }

  return ave;

}
float SumPower(float *data,int ilower,int iupper) {

  float ave=0.;
  for (int i=ilower;i<=iupper;i++) {
    ave+=data[i]*data[i];
  }

  return ave;

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

  RootStyle->SetOptFit     (1);

  RootStyle->SetMarkerStyle(20);
  RootStyle->SetMarkerSize(1.35);

  //  cout << ">> Style initialized with the Root Style!" << endl;
  //  cout << ">> " << modified << endl << endl;
  return RootStyle;
}
