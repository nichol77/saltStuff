#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################


#Site Specific  Flags
SYSINCLUDES	= -I/sw/include
SYSLIBS         = -L/sw/lib

#Global Thingies
ROOTCONFIG   := root-config
ARCH         := $(shell $(ROOTCONFIG) --arch)
PLATFORM     := $(shell $(ROOTCONFIG) --platform)

CC 	= g++
LD	= g++
SOFLAGS	= -shared
OBJSUF	= o
SRCSUF	= cxx
DLLSUF  = so
OPT	= -O2 -g --debug  -Wall -fPIC
OutPutOpt     = -o # keep whitespace after "-o"

#ROOT Compiler and Linker Flags
ROOTCFLAGS	:= $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS	:= $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS	:= $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS	:= $(shell $(ROOTCONFIG) --glibs)


#Architechture Dependent Flags
ifeq ($(ARCH),linuxx8664gcc)
# AMD Opteron and Intel EM64T (64 bit mode) Linux with gcc 3.x
CXX           = g++
CXXFLAGS      = $(OPT2) -Wall -fPIC
LD            = g++
LDFLAGS       = $(OPT2)
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linux)
# Linux with egcs, gcc 2.9x, gcc 3.x
CXX           = g++
CXXFLAGS      = $(OPT2) -Wall -fPIC
LD            = g++
LDFLAGS       = $(OPT2)
SOFLAGS       = -shared
endif

ifeq ($(ARCH),macosx)
# MacOS X with cc (GNU cc 2.95.2 and gcc 3.3)
MACOSX_MINOR := $(shell sw_vers | sed -n 's/ProductVersion://p' | cut -d . -f 2)
MACOSXTARGET := MACOSX_DEPLOYMENT_TARGET=10.$(MACOSX_MINOR)
ifeq ($(MACOSX_MINOR),5)
MACOSX_MINOR  = 4
endif
CXX           = g++
CXXFLAGS      = $(OPT2) -pipe -Wall -W -Woverloaded-virtual
LD            = g++
LDFLAGS       = $(OPT2) -bind_at_load
# The SOFLAGS will be used to create the .dylib,
# the .so will be created separately
DLLSUF       = dylib
UNDEFOPT      = dynamic_lookup
ifneq ($(MACOSX_MINOR),4)
ifneq ($(MACOSX_MINOR),3)
UNDEFOPT      = suppress
LD            = g++
endif
endif
SOFLAGS       = -dynamiclib -single_module -undefined $(UNDEFOPT)
endif

#Add ROOT Flags
CXXFLAGS     += $(ROOTCFLAGS) $(SYSINCLUDES)
LDFLAGS      += $(ROOTLDFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)





#ROOT stuff

ROOT_LIBRARY = libRyanFft.${DLLSUF} 
LIB_OBJS = FFTtools.o FFTWComplex.o  fftDict.o
CLASS_HEADERS =  FFTtools.h FFTWComplex.h

all : $(ROOT_LIBRARY) 

fftDict.C: $(CLASS_HEADERS)
	@echo "Generating dictionary ..."
	@ rm -f *Dict* 
	rootcint $@ -c $(INCLUDES) $(CLASS_HEADERS) LinkDef.h


$(ROOT_LIBRARY) : $(LIB_OBJS) 
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
	$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
ifeq ($(MACOSX_MINOR),4)
	ln -sf $@ $(subst .$(DLLSUF),.so,$@)
else
	$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
	 $(OutPutOpt) $(subst .$(DLLSUF),.so,$@)
endif
else
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIB_OBJS) -o $@
endif

%.$(OBJSUF) : %.$(SRCSUF)
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

%.$(OBJSUF) : %.C
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@


clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(ROOT_LIBRARY)
	@rm -f $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY))	
	@rm -f $(TEST)
