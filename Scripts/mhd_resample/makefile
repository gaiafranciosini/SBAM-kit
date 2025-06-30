#============== check platform ====================#
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
# do something Linux-y
SYSLIB=-lpthread
CCC=g++
endif
ifeq ($(UNAME), Darwin)
# do something Darwin-y
SYSLIB=-lpthread
CCC=clang++
endif
#==================================================#
CFLAGS = -g 
OPTLEV = -O2


default: resample

resample:
	$(CCC) $(CFLAGS) $(OPTLEV) mhd_resample.cpp Voxels.cpp mhd_io.cpp $(SYSLIB) -o mhd_resample.x



clean:
	@rm -f *.o *.x


#===================   ZIP   ====================#
today=`date +"%d%b%y"`

zip:
	zip -r ./Resample_$(today).zip makefile *.h *.cpp
