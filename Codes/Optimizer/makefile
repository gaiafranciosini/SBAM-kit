ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX           = g++ -O2 -std=c++11 -fdiagnostics-color
CXXFLAGS      = -g -Wall -fPIC
CXXFLAGS      += $(ROOTCFLAGS)

LD            =  g++ -fdiagnostics-color
LDFLAGS       = -g
SOFLAGS       = -shared $(ROOTGLIBS) -lSpectrum -dynamiclib

NGLIBB        += $(PWD)/libOptimizer.so 
NGLIBB        += $(ROOTGLIBS) 
NGLIBB        += -lMinuit -lSpectrum
GLIBB          = $(filter-out -lNew, $(NGLIBB))

optiRoot: optimizeDose.o Lib
	$(LD) $(LDFLAGS) -o optiRoot optimizeDose.o $(GLIBB)

optimizeDose.o: optimizeDose.cpp
	$(CXX) $(CXXFLAGS) -Icommon -c optimizeDose.cpp

Optimizer.o: Optimizer.cpp
	$(CXX) $(CXXFLAGS) -Icommon -c Optimizer.cpp

mhd_io.o:common/mhd_io.cpp
	$(CXX) $(CXXFLAGS) -c common/mhd_io.cpp

InputTools.o: common/InputTools.cpp
	$(CXX) $(CXXFLAGS) -c common/InputTools.cpp

Lib: Optimizer.o mhd_io.o InputTools.o
	$(CXX) $(SOFLAGS) Optimizer.o mhd_io.o InputTools.o -o libOptimizer.so  

clean:
	rm -f *.o  *.so optiRoot

