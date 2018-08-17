#
# General Makefile for the OfflineUser package
#
#

# Replace the wildcard expression with .cc file list if you do
# not want to compile all .cc files in this directory
#

TOPDIR= $(shell pwd)

BINDIR  = $(TOPDIR)/bin
LIBDIR  = $(TOPDIR)/lib
SRCDIR  = $(TOPDIR)/src
INCDIR  = $(TOPDIR)/include
OBJDIR  = $(TOPDIR)/obj
XMLDIR  = $(TOPDIR)/xml

USER_SRCS = $(wildcard *.cc $(SRCDIR)/*.cc)
USER_HEADERS = $(wildcard $(INCDIR)/*.h $(INCDIR)/*.hh)
USER_XMLS = $(patsubst %.xml.in,%.xml,$(wildcard *.xml.in $(XMLDIR)/*.xml.in))

files := $(foreach USER_SRCS,$(SRCDIR),$(USER_SRCS))


MYHEADERS   = $(INCDIR)/TScintHit.hh $(INCDIR)/TPMTHit.hh $(INCDIR)/TParticleSimData.hh $(INCDIR)/TStationSimData.hh $(INCDIR)/TEventSimData.hh $(INCDIR)/TrackPoint.hh $(INCDIR)/Track.hh $(INCDIR)/Cluster.hh $(INCDIR)/MuonDetector.hh
MYSOURCES   = $(SRCDIR)/TScintHit.cc $(SRCDIR)/TPMTHit.cc $(SRCDIR)/TParticleSimData.cc $(SRCDIR)/TStationSimData.cc $(SRCDIR)/TEventSimData.cc $(SRCDIR)/TrackPoint.cc $(SRCDIR)/Track.cc $(SRCDIR)/Cluster.cc $(SRCDIR)/MuonDetector.cc
MYDICTOBJ   = $(PWD)/TScintHitDict.o $(PWD)/TPMTHitDict.o $(PWD)/TParticleSimData.o $(PWD)/TStationSimData.o $(PWD)/TEventSimData.o $(PWD)/TrackPointDict.o $(PWD)/TrackDict.o $(PWD)/ClusterDict.o $(PWD)/MuonDetector.o
MYSOURCESOBJ= $(SRCDIR)/TScintHit.o $(SRCDIR)/TPMTHit.o $(SRCDIR)/TParticleSimData.o $(SRCDIR)/TStationSimData.o $(SRCDIR)/TEventSimData.o $(SRCDIR)/TrackPointDict.o $(SRCDIR)/TrackDict.o $(SRCDIR)/ClusterDict.o $(SRCDIR)/MuonDetector.o $(MYDICTOBJ)


## Get platform 32/64 bit
LBITS   = $(shell getconf LONG_BIT)

# Set executable a name
EXE = userAugerOffline
#
#############################################################

## You should not need to change anything below this line ###

.PHONY: all depend clean

ifdef AUGEROFFLINEROOT
  AUGEROFFLINECONFIG = $(AUGEROFFLINEROOT)/bin/auger-offline-config
else
  AUGEROFFLINECONFIG = auger-offline-config
endif

OBJS = $(USER_SRCS:.cc=.o)


######################################
###  CPPFLAGS & CXXFLAGS  ############
######################################
CPPFLAGS    = -I$(INCDIR)
CPPFLAGS   += $(shell $(AUGEROFFLINECONFIG) --cppflags)
#CPPFLAGS   += -I$(INCDIR)
CXXFLAGS    = $(shell $(AUGEROFFLINECONFIG) --cxxflags)
MAIN        = $(shell $(AUGEROFFLINECONFIG) --main)
CONFIGFILES = $(shell $(AUGEROFFLINECONFIG) --config)
XMLSCHEMALOCATION = $(shell $(AUGEROFFLINECONFIG) --schema-location)


###########################
###  LDFLAGS   ############
###########################
LDFLAGS     = $(shell $(AUGEROFFLINECONFIG) --ldflags)

ifeq ($(LBITS),64)
#  # do 64 bit stuff here
#	CPPFLAGS += -I/usr/include -pthread -m64
  CXXFLAGS = -O2 -Wall -Wunused -Wuninitialized -fPIC -pthread -m64
	SYSLIBDIR = /usr/lib64
else
#  # do 32 bit stuff here
#	CPPFLAGS += -I/usr/include -pthread -m32
  CXXFLAGS = -O2 -Wall -Wunused -Wuninitialized -fPIC -pthread -m32
	SYSLIBDIR = /usr/lib
endif

SOFLAGS = -fPIC -ggdb3 -Wall -shared

## additional system library
LDFLAGS += -L$(SYSLIBDIR) -lrt 

## additional ROOT library (TGeo)
LDFLAGS += -L$(ROOTSYS)/lib -lGeom -lGeomPainter -lMinuit

################################################################

##all: GETOBJS PRINTINFO ClassDict.cc $(EXE) libTHits.so PUTOBJS $(USER_XMLS)
all: GETOBJS PRINTINFO $(EXE) libTHits.so PUTOBJS $(USER_XMLS)


PRINTINFO: 
	@echo 'Compiling on a $(LBITS) bit machine' \

PUTOBJS:
	@echo "Moving objects to $(OBJDIR) dir"
	- mv -f $(SRCDIR)/*.o $(OBJDIR) 2>/dev/null
	- mv -f *.o $(OBJDIR) 2>/dev/null

GETOBJS:
	@echo "Put objects again to $(SRCDIR) dir"
	- mv -f $(OBJDIR)/*.o $(SRCDIR) 2>/dev/null
	- mv -f $(SRCDIR)/ClassDict.o $(TOPDIR) 2>/dev/null


##ClassDict.cc: $(MYHEADERS) LinkDef.h
##	rootcint -f $@ -c $(CXXFLAGS) -p $^

ClassDict.cc: $(MYHEADERS) LinkDef.h
	rootcint -f $@ -c $^
	
#libTHits.so: 
#	g++ $(CXXFLAGS) -shared -o $@ $(LDFLAGS) $(CXXFLAGS) $(ROOTCPPFLAGS) $^

#libTHits.so: ClassDict.cc $(OBJS)
#	@echo generating $@ shared library...
#	@g++ $(SOFLAGS) -o $@ $(LDFLAGS) $(CPPFLAGS) $^

libTHits.so: $(OBJS)
	@echo generating $@ shared library...
	@g++ $(SOFLAGS) -o $@ $(LDFLAGS) $(CPPFLAGS) $^
	

$(EXE): $(OBJS) ClassDict.o
	@echo "Compiler is $(CXX) or $(CC), options are $(CXXFLAGS)"
	@$(CXX) $(CXXFLAGS) -o $@ $^ $(MAIN) $(CXXFLAGS) $(LDFLAGS) -lMinuit


%: %.in
	@echo -n "Generating $@ file..." 
	@sed -e 's!@CONFIGDIR@!$(CONFIGFILES)!g;s!@SCHEMALOCATION@!$(XMLSCHEMALOCATION)!g' $< >$@
	@echo "done"

#############################################################
# gcc can generate the dependency list

depend: Make-depend

Make-depend: $(USER_SRCS)
	$(CPP) $(CPPFLAGS) -MM $^ > $@


clean:
	- rm -f *.o $(OBJDIR)/*.o $(SRCDIR)/*.o $(BINDIR)/*.o *.so $(MYDICTOBJ) ClassDict* *.ps core $(USER_XMLS) Make-depend

#############################################################
# 'make run' will run the thing

run: GETOBJS PRINTINFO ClassDict.cc $(EXE) libTHits.so PUTOBJS $(USER_XMLS) 
	./$(EXE) -b xml/bootstrap.xml && touch $@

#############################################################
# the lines below are for running with debugger 'make run_gdb'

.INTERMEDIATE: gdb.cmdl

# batch mode gdb needs a file with commands
gdb.cmdl:
	echo "r -b bootstrap.xml" > $@

run_gdb: gdb.cmdl $(EXE) $(USER_XMLS)
	gdb -batch -x $< ./$(EXE) && touch $@

include Make-depend
