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

smooth:
	$(CCC) $(CFLAGS) $(OPTLEV) mhd_smooth.cpp Voxels.cpp mhd_io.cpp $(SYSLIB) -o mhd_smooth.x

clean:
	@rm -f *.o *.x
