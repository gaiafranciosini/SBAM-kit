#SYSLIB=-lpthread
#============== check platform ====================#
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
# do something Linux-y
CCC=g++
endif
ifeq ($(UNAME), Darwin)
# do something Darwin-y
CCC=g++
endif
#==================================================#
CFLAGS = -g 
OPTLEV = -O2

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX           = $(CCC) -O2 -fdiagnostics-color
CXXFLAGS      = -Wall -fPIC
CXXFLAGS      += $(ROOTCFLAGS)

LD            = $(CCC) -fdiagnostics-color
LDFLAGS       = 
SOFLAGS       = -shared $(ROOTGLIBS) -lSpectrum -dynamiclib

NGLIBB        += $(ROOTGLIBS) 
NGLIBB        += -lMinuit -lSpectrum
GLIBB          = $(filter-out -lNew, $(NGLIBB))


#default: rebin_Dij smooth_Dij resample
default: Dij2Vec 

#resample:
#	$(CCC) $(CFLAGS) $(OPTLEV) mhd_resample.cpp Voxels.cpp mhd_io.cpp $(SYSLIB) -o mhd_resample.x

Dij2Vec:
	$(CXX) $(CFLAGS) $(CXXFLAGS) $(OPTLEV) $(LDFLAGS) Dij2Vec.cpp Voxels.cpp mhd_io.cpp $(SYSLIB) -o Dij2Vec.x $(GLIBB)

#rebin_Dij:
#	$(CCC) $(CFLAGS) $(OPTLEV) rebin_Dij.cpp Voxels.cpp mhd_io.cpp $(SYSLIB) -o rebin_Dij.x

clean:
	@rm -f *.o *.x


#===================   ZIP   ====================#
today=`date +"%d%b%y"`

zip:
	zip -r ./Resample_$(today).zip makefile *.h *.cpp
