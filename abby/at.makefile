CCFILE      = attenlength.cc

# Define filename suffixes
ObjSuf        = .o
SrcSuf        = .cc
ExeSuf        = .exe
DllSuf        = .so
OutPutOpt     = -o

# Define the compile and link commands
CXX           = g++
CXXFLAGS      = -c -g -o $(OBJ) -Wall -fPIC -I$(ROOTSYS)/include/
LD            = g++
LDFLAGS       = -O

# Define root libraries
ROOTLIBS      = `${ROOTSYS}/bin/root-config --libs` 

# Define all my libraries	
LIBS          = $(ROOTLIBS) 

# Define shortcuts for compiling and linking
COMPILE = $(CXX) $(CXXFLAGS) $(SCRATCH)
LINK = 	$(LD) $(LDFLAGS) $(OBJ) $(LIBS) $(OutPutOpt) $(EXE)
#------------------------------------------------------------------------------

# Define the file names
OBJ       = attenlength.obj
SRC       = $(CCFILE)
EXE       = attenlength.exe
SCRATCH	  = tmp.cc

all:            $(EXE)

# Update the executable if the object file has changed
$(EXE): 	$(OBJ)
		$(LINK)


# Update the object file if the source, or include file changed
$(OBJ):
		cat $(CCFILE) > $(SCRATCH)
		$(COMPILE)

clean:
	rm -f $(OBJ)	
