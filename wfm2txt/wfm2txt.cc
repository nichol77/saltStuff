/*
wfm2txt: a utility to convert Tektronix .WFM format binary files to ascii format.
Uses .WFM definition given in "TDS Waveform File Format", August 27, 1999, found on 
Tektronix website.  (Filename is "obsolete_scopes_WFMformat.pdf".)

Stephen Hoover 15 Dec 2006

Stephen Hoover 23 Aug 2007 : Increased array size to 300000, put in check that data size is not longer than array size.

** Known Issues / To Be Examined **
There appears to be a constant offset of 7.2e-11 seconds between time points from the scope's ascii output and from the .WFM output.  Verify?  Shouldn't be a problem, as the offset is constant and small.

Do I need to take vertGain / 256 instead of vertGain to find YOFF for a regular waveform? - Probably no, should be okay as is.

The program needs more testing, particularly with extended acquisitions.

*/

#include <cstdio>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>


//Define Tektronix structures used to store waveform header data.  
//Structures, with comments, copied directly from Tektronix document.
typedef struct tristarWfmHeader
{
  short cHeaderIndex;      /* c header index */
  short interpType;        /* various acq modes, interp type */
  short ttOffset;          /* ttOffset */
  short acqModeAndDigInfo; /* acqMode, interleave dig info, clip */
  short realPointOffset;   /* real point offset */
  short realPointSpacing;  /* real point spacing */
} __attribute__((packed)) TristarWfmHeader; /* 12 bytes */

typedef struct ealWfmHeader
{
  unsigned long acqLength;     /* Extended Acq Length, matches record length
				  if EAL mode is Off */
  double acqZoom;              /* Zoom factor used to produce the waveform
				  from the acquisition. */
  double acqTrigPos;           /* Trigger position percent within the 2, 4,
				  or 8 Meg record. Range is 0.0 to 1.0 */
  double recordStartInAcq;     /* Start point of recordLength in extended
				  record length in percent. */
  double acqTrigPosRelRecLen;  /* Trig pos of EAL as a percentage of the short
				  record length relative to beginning of 
				  the short record length. Units in percent. */
  short acqCharge;             /* Indicates the amount of pre and post charge data
				  in the acquisition. Typically 320 samples. */
  short acqState;              /* Stores whether the current acquisition has been
				  obtained by acquiring, recalling, is null, or
				  eal is off. States are ACQacquired, ACQrecalled,
				  ACQnull, and ACQealOff. */
  long expansionSpace[8];      /* 32 bytes of future expansion space. */
} __attribute__((packed)) EalWfmHeader; /* 72 bytes */

typedef struct wfmHeader
{
  short acqMode;            /* acquisition mode */
                            /* TOKEN - see SYMBOL table */
  short minMaxFormat;       /* minMax data format: ON | OFF */
                            /* TOKEN - see SYMBOL table */
  double duration;          /* duration of 1 acquisition in secs */
                            /* recordLength * horzScalePerPoint */
  short vertCpl;            /* vertical coupling ID */
                            /* TOKEN - see SYMBOL table */
  short horzUnit;           /* horizontal unit ID */
                            /* TOKEN - see SYMBOL table */
  double horzScalePerPoint; /* horizontal scale per point secs - XINCR */
                            /* horzScalePerPoint * 50 = time per Div in secs */
  short vertUnit;           /* vertical unit ID */
                            /* TOKEN - see SYMBOL table */
  double vertOffset;        /* vertical offset in volts - YZERO */
                            /* Gets added into digitized level to get real Volts */
                            /* Y = YZERO + YMULT * (Ydig - YOFF) */
  double vertPos;           /* vertical position in divs */
                            /* YOFF (in dig levels) = vertPos * 25 pts/div */
  double vertGain;          /* vertical gain in Volts per Div */
                            /* YMULT = vertGain/25 for single byte data */
                            /* YMULT = vertGain/25/256 for two byte data */
  unsigned long recordLength; /* record length in points */
                              /* NOT bytes - but points acquired */
                              /* NR_PT = recordLength */
  short trigPos;            /* trigger position in percent */
                            /* = (short) 100 * ptno / recordLength */
  short wfmHeaderVersion;   /* waveform header version: msbyte.lsbyte.
			       Was invertState in previous versions */
  short sampleDensity;      /* sample density */
  short burstSegmentLength; /* length of burst segment 0 = inactive */
  short sourceWfm;          /* Source waveform of saveref or first
			       source of math waveform */
  short videoTrigLineNumber;/* line number that video 
			       trigger was set to */
  short field;              /* video trig field number */
  short frameSize;          /* more video trig */
  double lineWidth;         /* yet more video */
  short standard;           /* video: NTSC, etc */
} __attribute__((packed)) WfmHeader; /* 80 bytes */

typedef struct referenceWfmHead
{
  double vertPos;
  double horzPos;
  double vertZoom;
  double horzZoom;
  WfmHeader cHeader;
  TristarWfmHeader tsHeader;
} __attribute__((packed)) ReferenceWfmHead;  /* 32 + 80 + 12 = 124 bytes */

typedef struct extendedReferenceWfmHead
{
  double vertPos;
  double horzPos;
  double vertZoom;
  double horzZoom;
  WfmHeader cHeader;
  EalWfmHeader ealHeader;
  TristarWfmHeader tsHeader;
} __attribute__((packed)) ExtendedReferenceWfmHead; /* 32 + 80 + 72 + 12 = 196 bytes */


void usage(char *execname) {
  printf("wfm2txt: This program converts TDS Waveform files (.WFM format) to ascii format.\n");
  printf("The program follows format guidelines given in a 27 August 1999 document.\n");
  printf("usage: %s [options] (input.WFM)\n", execname);
  printf("\tValid options are as follows:\n");
  printf("\t-o [outputfile] : Write output to file [outputfile].  Default is input.txt\n");
  printf("\t-c\tOutput data as comma separated values (default is space separated).\n");
  printf("\t-p\tPrint header information at top of output file.\n");
  printf("\t-q\tDo not print header information to screen.\n");
  printf("\t-n\tDo not print to file. (Use with -q or -d to examine a file.)\n");
  printf("\t-d\tdebug mode\n");
  printf("\t-h\tthis help menu\n");
  printf("\t-v\tversion\n");
  exit(1);
} //method usage

//Tektronix stores data with opposite Endianness from Intel processors.  Flip!
//Generic method to flip Endianness.
//WARNING: Flips byte order of anything put in - do not use on things like stuctures or classes!
template <class thing>
inline void flipEndian(thing &in) {
  int size = sizeof(thing);

  thing out;

  char* p_in = (char *) &in;
  char* p_out = (char *) &out;

  for(int i=0;i<size;i++) {
    p_out[i] = p_in[size-1-i];
  }

  in = out;

  return;
} //method flipEndian

void flipReferenceWfmHeadEndian(ReferenceWfmHead &in);
void flipExtReferenceWfmHeadEndian(ExtendedReferenceWfmHead &in);

std::string lookupSymbol(const int &token);

int main(int argc, char *argv[]) {
  int debug = 0; //Debug flag
  int print_header = 0;
  int screen_header = 1;
  int print_to_file = 1;
  std::string outputfilename = "";
  std::string output_extension = "txt";
  std::string output_spacer = " "; //Default is space separated values
  std::string inputfilename = "";
  std::string wfm_extension = "WFM";
  FILE *inputfile;
  std::ofstream output;
  long filesize = 0;
  int extended_wfm=-1; //Flag signifying that we're reading an extended waveform file
  int i=0;

  char clswitch;
  
  if(argc > 1) {//read arguments with getopts
    while((clswitch = getopt(argc, argv, "dhvo:cpqn")) != EOF) {
      switch(clswitch) {
      case 'd'://debug mode
	debug = 1;
	printf("DEBUG MODE\n\n");
	break;
      case 'o'://Use user specified output file
	outputfilename = optarg;
	break;
      case 'c'://Use comma separated values in output file
	output_spacer = ",";
	break;
      case 'p'://Print header info in output file
	print_header = 1;
	break;
      case 'q'://Do not print header info on screen
	screen_header = 0;
	break;
      case 'n'://Do not print to file
	print_to_file = 0;
	break;
      case 'h'://help menu
	usage(argv[0]);
	break;
      case 'v'://version
	printf("wfm2txt version 1.1.1 beta, 23 Aug 2007\n");
	exit(1);
	break;
      }//switch
    }//while
  } else {
    usage(argv[0]);
    exit(1);
  }//else
  
  //Give the usage to the user if no input file was provided.
  if(optind == argc) usage(argv[0]);//must pass filename
  
  if(debug && outputfilename != "") std::cout<<"DEBUG: Using output file "<<outputfilename<<std::endl;

  //Separate out the base and extension of the input filename
  char *basename, *extension;
  if(!strchr(argv[argc-1], '.')) 
    {
      std::cout << "Error: file must be of type " << wfm_extension << std::endl;
      exit(1);
    }//if
  basename = strtok(argv[argc-1], ".");
  extension = strtok(NULL, "\n");

  //Check to make sure that the input file extension marks it as a Tektronix waveform file (.WFM)
  if(extension != wfm_extension) 
    {
      std::cout << "Error: file must be of type " << wfm_extension << std::endl;
      exit(1);
    }//if

  inputfilename = std::string(basename) + "." + std::string(extension);
  if (debug) std::cout<<"DEBUG: Reading file "<<inputfilename<<"\n";

  //If the user does not supply an output file, use the basename of the input file with a .txt extension
  if(outputfilename == "")
    {
      outputfilename = std::string(basename) + "." + output_extension;
    }//if

  //Open input file
  if(!(inputfile = fopen(inputfilename.c_str(), "rb"))) 
    {
      std::cout << "Error: unable to open " << inputfilename << std::endl;
      exit(1);
    }//if
  
  if(debug)//Find filesize - only needed for debugging purposes, really.
    {
      fseek(inputfile , 0 , SEEK_END);
      filesize = ftell(inputfile);
      rewind(inputfile);
  
      std::cout<<"DEBUG: filesize = "<<filesize<<std::endl;
    }//if

  //Declare variables to store WFM data
  const int array_size = 300000;
  char start_of_header[8];
  char tmp[9];
  int digits_in_byte_count = 0;
  long byte_count = 0; //in bytes, obviously
  long magic_number = 0; //Some secret Tektronix thing - we don't care what it is.
  long length_of_waveform = 0; //in bytes
  double x[array_size]; //for storing processed data
  double y[array_size]; //for storing processed data

  //Variables for regular waveforms
  referenceWfmHead header;
  short int preamble[16];
  short int data[array_size];
  short int postamble[16];
  unsigned short int checksum = 0; //used for both waveforms
  unsigned short int checksum_check = 0; //my calculation of checksum, used for both waveforms
  unsigned short int *cs_ptr = NULL; //Used to find checksum contribution from header, used for both waveforms

  //Variables for extended waveforms
  extendedReferenceWfmHead ext_header;
  char ext_preamble[1000];
  char ext_postamble[1000];
  char ext_data[array_size];

  //Zero all variables
  memset(start_of_header, 0, 8);
  memset(tmp, 0, 9);
  memset(preamble, 0, 16);
  memset(data, 0, array_size);
  memset(x, 0, array_size);
  memset(y, 0, array_size);
  memset(postamble, 0, 16);
  memset(ext_preamble, 0, 1000);
  memset(ext_postamble, 0, 1000);
  memset(ext_data, 0, array_size);

  //cout<<"sizeof(int) = "<<sizeof(int)<<", sizeof(long) = "<<sizeof(long)<<", sizeof(long int) = "<<sizeof(long int)<<", sizeof(unsigned long) = "<<sizeof(unsigned long)<<std::endl;
  //cout<<"sizeof(char) = "<<sizeof(char)<<", sizeof(referenceWfmHead) = "<<sizeof(referenceWfmHead)<<std::endl;
  //cout<<"sizeof(extendedReferenceWfmHead) = "<<sizeof(extendedReferenceWfmHead)<<", sizeof(tristarWfmHeader) = "<<sizeof(tristarWfmHeader)<<std::endl;
  //cout<<"sizeof(ealWfmHeader) = "<<sizeof(ealWfmHeader)<<", sizeof(wfmHeader) = "<<sizeof(wfmHeader)<<std::endl;


  //Read start of header.  It is either 7 or 8 chars, and always ends with a '#'.
  fread(start_of_header, sizeof(char), 7, inputfile);
  if(start_of_header[6] != '#') 
    {
      rewind(inputfile);
      fread(start_of_header, sizeof(char), 8, inputfile);
    }
  if(debug) std::cout<<"DEBUG: start_of_header = "<<start_of_header<<std::endl;

  //The start of the header tells us if this is a regular or extended waveform.  
  //"LLWFM #" for regular, "LLWF1 #" for extended
  if (std::string(start_of_header) == "LLWFM #" || std::string(start_of_header) == ":LLWFM #") extended_wfm = 0;
  else if (std::string(start_of_header) == "LLWF1 #" || std::string(start_of_header) == ":LLWF1 #") extended_wfm = 1;
  if(debug) std::cout<<"DEBUG: extended_wfm = "<<extended_wfm<<std::endl;
  
  //The byte count is ascii, so we first read the number of digits in the byte count...
  fread(&tmp, sizeof(char), 1, inputfile);
  digits_in_byte_count = atoi(tmp);
  memset(tmp, 0, 9);
  if(debug) std::cout<<"DEBUG: digits_in_byte_count = "<<digits_in_byte_count<<std::endl;

  //...and then read that many digits worth of chars.
  //Byte count is the number of bytes following the byte count
  fread(tmp, sizeof(char), digits_in_byte_count, inputfile);
  byte_count = atoi(tmp);
  memset(tmp, 0, 9);
  if(debug) std::cout<<"DEBUG: byte_count = "<<byte_count<<std::endl;

  //This is the magic number that we don't need.  Ignore.
  fread(&magic_number, sizeof(magic_number), 1, inputfile);
  flipEndian(magic_number);
  if(debug) std::cout<<"DEBUG: magic_number = "<<magic_number<<std::endl;

  //This is supposed to be the number of bytes following.  I'm sceptical.
  fread(&length_of_waveform, sizeof(length_of_waveform), 1, inputfile);
  flipEndian(length_of_waveform);
  if(debug) std::cout<<"DEBUG: length_of_waveform = "<<length_of_waveform<<std::endl;

  //Read the header structure
  if (!extended_wfm) 
    {
      fread(&header, sizeof(header), 1, inputfile);
      flipReferenceWfmHeadEndian(header);
      if(debug) std::cout<<"DEBUG: vertPos = "<<header.vertPos<<std::endl;
      if(debug) std::cout<<"DEBUG: horzPos = "<<header.horzPos<<std::endl;
      if(debug) std::cout<<"DEBUG: vertZoom = "<<header.vertZoom<<std::endl;
      if(debug) std::cout<<"DEBUG: horzZoom = "<<header.horzZoom<<std::endl;
      if(debug) std::cout<<"DEBUG: acqMode = "<<lookupSymbol(header.cHeader.acqMode)<<std::endl;
      if(debug) std::cout<<"DEBUG: mimMaxFormat = "<<lookupSymbol(header.cHeader.minMaxFormat)<<std::endl;
      if(debug) std::cout<<"DEBUG: recordLength = "<<header.cHeader.recordLength<<std::endl;
      if(debug) std::cout<<"DEBUG: horzUnit = "<<lookupSymbol(header.cHeader.horzUnit)<<std::endl;
      if(debug) std::cout<<"DEBUG: horzScalePerPoint = "<<header.cHeader.horzScalePerPoint <<std::endl;
      if(debug) std::cout<<"DEBUG: duration = "<<header.cHeader.duration <<std::endl;
      if(debug) std::cout<<"DEBUG: trigPos = "<<header.cHeader.trigPos<<std::endl;
      if(debug) std::cout<<"DEBUG: realPointOffset = "<<header.tsHeader.realPointOffset <<std::endl;
      if(debug) std::cout<<"DEBUG: realPointSpacing = "<<header.tsHeader.realPointSpacing <<std::endl;
      if(debug) std::cout<<"DEBUG: sourceWfm = "<<lookupSymbol(header.cHeader.sourceWfm)<<std::endl;
      if(debug) std::cout<<"DEBUG: interpType = "<<header.tsHeader.interpType<<std::endl;
      
      if (int(header.cHeader.recordLength) > array_size)
	{
	  std::cout<<"Error!  WFM file has "<<header.cHeader.recordLength<<" data points, but my arrays will only hold "<<array_size<<"!  Exiting.\n";
	  exit(1);
	}
    } //if(!extended_wfm)
  else 
    {
      fread(&ext_header, sizeof(ext_header), 1, inputfile);
      flipExtReferenceWfmHeadEndian(ext_header);
      if(debug) std::cout<<"DEBUG: vertPos = "<<ext_header.vertPos<<std::endl;
      if(debug) std::cout<<"DEBUG: horzPos = "<<ext_header.horzPos<<std::endl;
      if(debug) std::cout<<"DEBUG: acqMode = "<<ext_header.cHeader.acqMode<<std::endl;
      if(debug) std::cout<<"DEBUG: recordLength = "<<ext_header.cHeader.recordLength<<std::endl;
    }//else

  //Now read the data, with a preamble before and postamble after.  I have no idea what
  //the preamble and postamble are, and they don't seem to be important.
  if (!extended_wfm) 
    {
      fread(preamble, sizeof(short int), 16, inputfile);
      for(i=0;i<16;i++) flipEndian(preamble[i]);

      fread(data, sizeof(short int), header.cHeader.recordLength, inputfile);
      for(i=0;i<int(header.cHeader.recordLength);i++) flipEndian(data[i]);

      fread(postamble, sizeof(short int), 16, inputfile);
      for(i=0;i<16;i++) flipEndian(postamble[i]);

      fread(&checksum, sizeof(checksum), 1, inputfile);
      flipEndian(checksum);
      if(debug) std::cout<<"DEBUG: checksum = "<<checksum<<std::endl;
      
      //Calculate our own checksum value.
      cs_ptr = (unsigned short int *)&header;
      for(i=0; i<(124/2); i++) checksum_check += cs_ptr[i]; //124 = size of a referenceWfmHead in bytes, each unsigned short int is 2 bytes
      for(i=0;i<16;i++) checksum_check += (unsigned short)preamble[i];
      for(i=0;i<int(header.cHeader.recordLength);i++) checksum_check += (unsigned short)data[i];
      for(i=0;i<16;i++) checksum_check += (unsigned short)postamble[i];
      if(debug) std::cout<<"DEBUG: checksum_check = "<<checksum_check<<std::endl;

      //Calculate the actual data values from the binary data.
      for(i=0;i<int(header.cHeader.recordLength);i++) 
	{
	  x[i] = ((double) i - ((double) header.cHeader.recordLength*(double) header.cHeader.trigPos / 100.0)) * header.cHeader.horzScalePerPoint;
	  y[i] = ((double) data[i] * (header.cHeader.vertGain / 25.0 / 256.0)) + header.cHeader.vertOffset - (header.cHeader.vertPos * header.cHeader.vertGain);
	}
    } //if (!extended_wfm)

  else
    {
      fread(ext_preamble, sizeof(char), ext_header.ealHeader.acqCharge, inputfile);
      fread(ext_data, sizeof(char), ext_header.cHeader.recordLength, inputfile);
      fread(ext_postamble, sizeof(char), ext_header.ealHeader.acqCharge, inputfile);
      fread(&checksum, sizeof(checksum), 1, inputfile);
      flipEndian(checksum);
      if(debug) std::cout<<"DEBUG: checksum = "<<checksum<<std::endl;

      //Calculate our own checksum value.
      cs_ptr = (unsigned short int *)&ext_header;
      for(i=0; i<(196/2); i++) checksum_check += cs_ptr[i];//196 = size of an extendedReferenceWfmHead in bytes, each unsigned short int is 2 bytes

      cs_ptr = (unsigned short int *)ext_preamble;
      for(i=0; i<(ext_header.ealHeader.acqCharge/2); i++) checksum_check += cs_ptr[i];
  
      cs_ptr = (unsigned short int *)ext_data;
      for(i=0; i<int(ext_header.cHeader.recordLength/2); i++) checksum_check += cs_ptr[i];

      cs_ptr = (unsigned short int *)ext_postamble;
      for(i=0;i<(ext_header.ealHeader.acqCharge/2);i++) checksum_check += cs_ptr[i];
      if(debug) std::cout<<"DEBUG: checksum_check = "<<checksum_check<<std::endl;

      //Calculate the actual data values from the binary data.
      for(i=0;i<int(header.cHeader.recordLength);i++) 
	{
	  x[i] = ((double) i - ((double)ext_header.ealHeader.acqLength * ext_header.ealHeader.acqTrigPos)) * ext_header.cHeader.horzScalePerPoint;
	  y[i] = ((double)ext_data[i] * (ext_header.cHeader.vertGain / 25.0)) + ext_header.cHeader.vertOffset - (ext_header.cHeader.vertPos * ext_header.cHeader.vertGain);
	}
    } //else (extended_wfm)

  //Verify that the checksum is correct.  If not, assume file is corrupted and exit program.
  if(checksum != checksum_check)
    {
      std::cout<<"Error!  Invalid checksum!\n";
      exit(1);
    }

  //Open output file
  output.open(outputfilename.c_str());
  if(!output.is_open()) 
    {
      std::cout<<"Error!  Could not open output file "<<outputfilename<<std::endl;
      exit(1);
    }  

  //Print data to the output file
  if(!extended_wfm)
    {
      if(screen_header) //Print header info to screen if desired.
	{
	  std::cout<<"Waveform taken from "<<lookupSymbol(header.cHeader.sourceWfm)<<" with coupling set to "<<lookupSymbol(header.cHeader.vertCpl)<<"."<<std::endl;
	  std::cout<<header.cHeader.recordLength<<" points recorded in "<<lookupSymbol(header.cHeader.acqMode)<<" mode."<<std::endl;
	  std::cout<<"X-axis (left column) in units of "<<lookupSymbol(header.cHeader.horzUnit)<<", Y-axis (right column) in units of "<<lookupSymbol(header.cHeader.vertUnit)<<std::endl;
	  std::cout<<"X-axis scale was "<<header.cHeader.horzScalePerPoint*50.<<" "<<lookupSymbol(header.cHeader.horzUnit)<<"/Div, Y-axis scale was "<<header.cHeader.vertGain<<" "<<lookupSymbol(header.cHeader.vertUnit)<<"/Div.\n";
	  std::cout<<"Vertical zoom = "<<header.vertZoom<<", horizontal zoom = "<<header.horzZoom<<std::endl<<std::endl;
	} //if(screen_header)

      if(print_header) //Print header info to file if desired.
	{
	  output<<"Waveform taken from "<<lookupSymbol(header.cHeader.sourceWfm)<<" with coupling set to "<<lookupSymbol(header.cHeader.vertCpl)<<"."<<std::endl;
	  output<<header.cHeader.recordLength<<" points recorded in "<<lookupSymbol(header.cHeader.acqMode)<<" mode."<<std::endl;
	  output<<"X-axis (left column) in units of "<<lookupSymbol(header.cHeader.horzUnit)<<", Y-axis (right column) in units of "<<lookupSymbol(header.cHeader.vertUnit)<<std::endl;
	  output<<"X-axis scale was "<<header.cHeader.horzScalePerPoint*50.<<" "<<lookupSymbol(header.cHeader.horzUnit)<<"/Div, Y-axis scale was "<<header.cHeader.vertGain<<" "<<lookupSymbol(header.cHeader.vertUnit)<<"/Div.\n";
	  output<<"Vertical zoom = "<<header.vertZoom<<", horizontal zoom = "<<header.horzZoom<<std::endl;
	  output<<"###START WAVEFORM###"<<std::endl;
	} //if(print_header)

      if (print_to_file)  //Print data to file if desired
	{
	  for(i=0;i<int(header.cHeader.recordLength);i++) 
	    {
	      output<<x[i]<<output_spacer<<y[i]<<std::endl;
	    } //for (print waveform to file)
	}
    }//if(regular waveform)
  else //if it is an extended waveform
    {
      if(screen_header) //Print header info to screen if desired.
	{
	  std::cout<<"Extended waveform taken from "<<lookupSymbol(ext_header.cHeader.sourceWfm)<<" with coupling set to "<<lookupSymbol(ext_header.cHeader.vertCpl)<<"."<<std::endl;
	  std::cout<<ext_header.cHeader.recordLength<<" points recorded in "<<lookupSymbol(ext_header.cHeader.acqMode)<<" mode."<<std::endl;
	  std::cout<<"Extended acquisition length: "<<ext_header.ealHeader.acqLength<<", Zoom factor: "<<ext_header.ealHeader.acqZoom<<std::endl;
	  std::cout<<"Start point of recordLength in extended record: "<<ext_header.ealHeader.recordStartInAcq<<"%."<<std::endl;
	  std::cout<<"X-axis (left column) in units of "<<lookupSymbol(ext_header.cHeader.horzUnit)<<", Y-axis (right column) in units of "<<lookupSymbol(ext_header.cHeader.vertUnit)<<std::endl;
	  std::cout<<"X-axis scale was "<<ext_header.cHeader.horzScalePerPoint*50.<<" "<<lookupSymbol(ext_header.cHeader.horzUnit)<<"/Div, Y-axis scale was "<<ext_header.cHeader.vertGain<<" "<<lookupSymbol(ext_header.cHeader.vertUnit)<<"/Div.\n";
	  std::cout<<"Vertical zoom = "<<ext_header.vertZoom<<", horizontal zoom = "<<ext_header.horzZoom<<std::endl<<std::endl;
	} //if(screen_header)

      if(print_header) //Print header info to file if desired
	{
	  output<<"Waveform taken from "<<lookupSymbol(ext_header.cHeader.sourceWfm)<<" with coupling set to "<<lookupSymbol(ext_header.cHeader.vertCpl)<<"."<<std::endl;
	  output<<ext_header.cHeader.recordLength<<" points recorded in "<<lookupSymbol(ext_header.cHeader.acqMode)<<" mode."<<std::endl;
	  output<<"Extended acquisition length: "<<ext_header.ealHeader.acqLength<<", Zoom factor: "<<ext_header.ealHeader.acqZoom<<std::endl;
	  output<<"Start point of recordLength in extended record: "<<ext_header.ealHeader.recordStartInAcq<<"%."<<std::endl;
	  output<<"X-axis (left column) in units of "<<lookupSymbol(ext_header.cHeader.horzUnit)<<", Y-axis (right column) in units of "<<lookupSymbol(ext_header.cHeader.vertUnit)<<std::endl;
	  output<<"X-axis scale was "<<ext_header.cHeader.horzScalePerPoint*50.<<" "<<lookupSymbol(ext_header.cHeader.horzUnit)<<"/Div, Y-axis scale was "<<ext_header.cHeader.vertGain<<" "<<lookupSymbol(ext_header.cHeader.vertUnit)<<"/Div.\n";
	  output<<"Vertical zoom = "<<ext_header.vertZoom<<", horizontal zoom = "<<ext_header.horzZoom<<std::endl;
	  output<<"###START WAVEFORM###"<<std::endl;
	} //if(print_header)

      if (print_to_file) //Print data to file if desired
	{ 
	  for(i=0;i<int(ext_header.cHeader.recordLength);i++) 
	    {
	      output<<x[i]<<output_spacer<<y[i]<<std::endl;
	    } //for (print waveform to file)
	} //if(print to file
    } //else (if it is an extended waveform)

  fclose(inputfile);
  output.close();
}//main

//This method converts integer tokens used by Tektronix into their proper meaning "On", "V", etc.
std::string lookupSymbol(const int &token) {
  std::string symbol = "Unknown token";
  char ref[5];
  memset(ref, 0, 5);

  switch (token) {
  case 97:
    symbol = "On";
    break;
  case 98:
    symbol = "Off";
    break;
  case 285:
    symbol = "Sample";
    break;
  case 2:
    symbol = "Peak Detect";
    break;
  case 3:
    symbol = "Hi Res";
    break;
  case 4:
    symbol = "Average";
    break;
  case 187:
    symbol = "Envelope";
    break;
  case 669:
    symbol = "Normal";
    break;
  case 930:
    symbol = "RMS Average";
    break;
  case 931:
    symbol = "Swept Spectrum";
    break;
  case 1132:
    symbol = "Transient";
    break;
  case 565:
    symbol = "DC Coupling";
    break;
  case 566:
    symbol = "AC Coupling";
    break;
  case 25:
    symbol = "Ground";
    break;
  case 609:
    symbol = "Volts";
    break;
  case 610:
    symbol = "S";
    break;
  case 632:
    symbol = "VV";
    break;
  case 626:
    symbol = "Valid";
    break;
  case 627:
    symbol = "Invalid";
    break;
  case 631:
    symbol = "Unknown";
    break;
  case 736:
    symbol = "Hz";
    break;
  case 740:
    symbol = "dB";
    break;
  case 766:
    symbol = "Vs";
    break;
  case 767:
    symbol = "VVs";
    break;
  case 768:
    symbol = "V/s";
    break;
  case 920:
    symbol = "V/v";
    break;
  case 107:
    symbol = "Ch1";
    break;
  case 108:
    symbol = "Ch2";
    break;
  case 109:
    symbol = "Ch3";
    break;
  case 110:
    symbol = "Ch4";
    break;
  case 907:
    symbol = "Ch5";
    break;
  case 908:
    symbol = "Ch6";
    break;
  case 909:
    symbol = "Ch7";
    break;
  case 910:
    symbol = "Ch8";
    break;
  case 911:
    symbol = "Ch9";
    break;
  case 912:
    symbol = "Ch10";
    break;
  case 913:
    symbol = "Ch11";
    break;
  case 914:
    symbol = "Ch12";
    break;
  case 915:
    symbol = "Ch13";
    break;
  case 916:
    symbol = "Ch14";
    break;
  case 917:
    symbol = "Ch15";
    break;
  case 918:
    symbol = "Ch16";
    break;
  case 111:
    symbol = "Math1";
    break;
  case 112:
    symbol = "Math2";
    break;
  case 113:
    symbol = "Math3";
    break;
  case 972:
    symbol = "Math4";
    break;
  case 973:
    symbol = "Math5";
    break;
  case 974:
    symbol = "Math6";
    break;
  case 975:
    symbol = "Math7";
    break;
  case 114:
    symbol = "Ref1";
    break;
  case 115:
    symbol = "Ref2";
    break;
  case 116:
    symbol = "Ref3";
    break;
  case 117:
    symbol = "Ref4";
    break;
  case 686:
    symbol = "NTSC";
    break;
  case 687:
    symbol = "PAL";
    break;
  case 723:
    symbol = "SECAM";
    break;
  case 688:
    symbol = "Custom";
    break;
  case 998:
    symbol = "Flex Format";
    break;
  case 936:
    symbol = "525 NTSC";
    break;
  case 937:
    symbol = "525 Mono";
    break;
  case 938:
    symbol = "625 PAL";
    break;
  case 939:
    symbol = "625 Mono";
    break;
  case 983:
    symbol = "HDTV";
    break;
  default:
    symbol = "Unknown token";
    break;
  } //switch(token)

  if (token >= 1140 && token <= 1199) {
    sprintf(ref,"Ref%i",token - 1135);
    symbol = std::string(ref);
  }

  return symbol;
} //method lookupSymbol

//Flips Endianness of all values in a regular waveform header
void flipReferenceWfmHeadEndian(ReferenceWfmHead &in) {

//   flipEndian(in->vertPos);
//   flipEndian(in->horzPos);
//   flipEndian(in->vertZoom);
//   flipEndian(in->horzZoom);
//   flipEndian(in->cHeader.acqMode);
//   flipEndian(in->cHeader.minMaxFormat);
//   flipEndian(in->cHeader.duration);
//   flipEndian(in->cHeader.vertCpl);
//   flipEndian(in->cHeader.horzUnit);
//   flipEndian(in->cHeader.horzScalePerPoint);
//   flipEndian(in->cHeader.vertUnit);
//   flipEndian(in->cHeader.vertOffset);
//   flipEndian(in->cHeader.vertPos);
//   flipEndian(in->cHeader.vertGain);
//   flipEndian(in->cHeader.recordLength);
//   flipEndian(in->cHeader.trigPos);
//   flipEndian(in->cHeader.wfmHeaderVersion);
//   flipEndian(in->cHeader.sampleDensity);
//   flipEndian(in->cHeader.burstSegmentLength);
//   flipEndian(in->cHeader.sourceWfm);
//   flipEndian(in->cHeader.videoTrigLineNumber);
//   flipEndian(in->cHeader.field);
//   flipEndian(in->cHeader.frameSize);
//   flipEndian(in->cHeader.lineWidth);
//   flipEndian(in->cHeader.standard);
//   flipEndian(in->tsHeader.cHeaderIndex);
//   flipEndian(in->tsHeader.interpType);
//   flipEndian(in->tsHeader.ttOffset);
//   flipEndian(in->tsHeader.acqModeAndDigInfo);
//   flipEndian(in->tsHeader.realPointOffset);
//   flipEndian(in->tsHeader.realPointSpacing);

  return;
} //method flipReferenceWfmHeadEndian

//Flips Endianness of all values in an extended waveform header.
void flipExtReferenceWfmHeadEndian(ExtendedReferenceWfmHead &in) {

//   flipEndian(in.ealHeader.acqLength);
//   flipEndian(in.ealHeader.acqZoom);
//   flipEndian(in.ealHeader.acqTrigPos);
//   flipEndian(in.ealHeader.recordStartInAcq);
//   flipEndian(in.ealHeader.acqTrigPosRelRecLen);
//   flipEndian(in.ealHeader.acqCharge);
//   flipEndian(in.ealHeader.acqState);

//   flipEndian(in.vertPos);
//   flipEndian(in.horzPos);
//   flipEndian(in.vertZoom);
//   flipEndian(in.horzZoom);
//   flipEndian(in.cHeader.acqMode);
//   flipEndian(in.cHeader.minMaxFormat);
//   flipEndian(in.cHeader.duration);
//   flipEndian(in.cHeader.vertCpl);
//   flipEndian(in.cHeader.horzUnit);
//   flipEndian(in.cHeader.horzScalePerPoint);
//   flipEndian(in.cHeader.vertUnit);
//   flipEndian(in.cHeader.vertOffset);
//   flipEndian(in.cHeader.vertPos);
//   flipEndian(in.cHeader.vertGain);
//   flipEndian(in.cHeader.recordLength);
//   flipEndian(in.cHeader.trigPos);
//   flipEndian(in.cHeader.wfmHeaderVersion);
//   flipEndian(in.cHeader.sampleDensity);
//   flipEndian(in.cHeader.burstSegmentLength);
//   flipEndian(in.cHeader.sourceWfm);
//   flipEndian(in.cHeader.videoTrigLineNumber);
//   flipEndian(in.cHeader.field);
//   flipEndian(in.cHeader.frameSize);
//   flipEndian(in.cHeader.lineWidth);
//   flipEndian(in.cHeader.standard);
//   flipEndian(in.tsHeader.cHeaderIndex);
//   flipEndian(in.tsHeader.interpType);
//   flipEndian(in.tsHeader.ttOffset);
//   flipEndian(in.tsHeader.acqModeAndDigInfo);
//   flipEndian(in.tsHeader.realPointOffset);
//   flipEndian(in.tsHeader.realPointSpacing);

  return;
} //method flipExtReferenceWfmHeadEndian
