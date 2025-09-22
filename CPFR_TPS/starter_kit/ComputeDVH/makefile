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

default: ComputeDVH


ComputeDVH:
	$(CCC) $(CFLAGS) $(OPTLEV) ComputeDVH.cpp mhd_io.cpp InputTools.cpp  $(SYSLIB) -o ComputeDVH.x

clean:
	@rm -f *.o *.x


