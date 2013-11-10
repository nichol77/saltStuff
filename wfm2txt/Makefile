INCLUDE		=
CXX		= g++
CXXFLAGS	+= $(INCLUDE) -O0 -g -Wall -fPIC -Wno-deprecated
SOFLAGS		= -shared -g
LD		= g++
LDFLAGS		= -O0 -g
TOPDIR		= /usr/local
BINDIR		= $(TOPDIR)/bin

OBJS = wfm2txt.o

EXES = $(OBJS:.o=)

all:	$(EXES)

%.o:	%.cc
	$(CXX) $(CXXFLAGS) -c $<

wfm2txt:	wfm2txt.o
			$(LD) $(LDFLAGS) $< -o $@

clean:
	@rm -f core
	@rm -f $(OBJS)
	@rm -f $(EXES)

install:
	@cp $(EXES) $(BINDIR)
	@echo $(EXES) installed to $(BINDIR)

print:
	@echo $(EXES)
