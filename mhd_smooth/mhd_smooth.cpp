 /*
 * 
 *  Fred project
 *
 *  Copyright 2011-2018 Università di Roma LA SAPIENZA. All rights reserved.
 *  A. Schiavi (2018)
 */

 
const char *version = "1.0";
/*			version history
1.0 : initial version (Sep 2020)
*/
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cstdarg>
#include <ctime>
#include <vector>
#include <algorithm>
using namespace std;
#include "types.h"
#include "vector3d.h"
#include <pthread.h>

#include "Voxels.h"

template<>
int vector3d<float32>::outMode=1;

template<>
int vector3d<float64>::outMode=1;

inline int el2d(int i,int j,int n,int m) {return i+j*n;}
inline int el3d(int i,int j,int k,int nn[3]) {return i+j*nn[0]+k*nn[0]*nn[1];} // column-major and contiguous 

ifstream fin;	

vec3dRT xmin_A, xmin_B;
i3d nn_A, nn_B;
vec3dRT hs_A, hs_B;
vec3dRT xmax_A, xmax_B;

double * map_A;
double * map_B;

void evaluationLooper();

typedef struct threadInfo_s {
	int thread_num;
	int num_threads;
} threadInfo;

void *kernelOnCPU(void *info);
int pthreads_max_num = 4;

vector<const char *> args;

int imode; // 1=hotspotcleaner ; 2=arithmeticmean
float hotspotThreshold = 5;

//=====================================================================
Voxels voxelsA,voxelsB;
double *mapA,*mapB;

//=====================================================================
void banner(){
	/* cout<<endl; */
	/* cout<<"+------------------------------------------------------------------------------+"<<endl; */
	/* cout<<"|                       Smooth 3D map (Fred project)                           |"<<endl; */
	/* cout<<"|                                                                              |"<<endl; */
	/* cout<<"|                             A.Schiavi -  2020                                |"<<endl; */
	/* cout<<"+------------------------------------------------------------------------------+"<<endl; */
	/* cout<<endl; */
	/* cout<<"                             version "<<version<<endl; */
	/* cout<<endl; */
}
//=====================================================================
void usage(){
	cout<<"usage: mapA -o newMap"<<endl;
	cout<<"applies a smoothing kernel to voxels of mapA"<<endl;
	// cout<<"options: "<<endl;
	// cout<<"\t-harmonic [default]: apply harmonic mean smooth "<<endl;
}
//=====================================================================
int checkOption(vector<const char *> args,const char * str){
  for (int i=0;i<args.size();i++){
    if(strncmp(args[i],str,strlen(str))==0) return i;
  }
  return 0;
}

//=====================================================================
const char normalcolor[]  = "\e[0m";
const char blackcolor[]   = "\e[0;30m";
const char redcolor[]     = "\e[0;31m";
const char greencolor[]   = "\e[0;32m";
const char yellowcolor[]  = "\e[0;33m";
const char bluecolor[]    = "\e[0;34m";
const char magentacolor[] = "\e[0;35m";
const char cyancolor[]    = "\e[0;36m";
const char whitecolor[]   = "\e[0;37m";

//=====================================================================
void logError(int n_args, ...){
  cerr<<redcolor;

  va_list ap;
  va_start(ap, n_args);
  for(int i = 1; i <= n_args; i++) { cerr<<va_arg(ap,char *); }
  va_end(ap);
  
  cerr<<normalcolor<<endl;
}
//=====================================================================
void getGrid(Voxels &Vxl,i3d &nn,vec3dRT &hs, vec3dRT &x0){
	 nn = Vxl.nn;
	 hs = Vxl.hs;
	 x0 = Vxl.x0;
}
//=====================================================================
void importVoxelsToBuffer(Voxels &Vxl,double **_dstBuffer){
	double *dstBuffer = new double[Vxl.N];
	if(!dstBuffer) {logError(1,"failed memory allocation"); exit(1);}
	*_dstBuffer = dstBuffer;

	if(Vxl.dtype == "MET_CHAR") { int8   *p = (int8  *) Vxl.data; for(size_t i =0; i<Vxl.N; i++) { dstBuffer[i]=p[i]; } ; return;}
	if(Vxl.dtype == "MET_SHORT") { int16  *p = (int16 *) Vxl.data; for(size_t i =0; i<Vxl.N; i++) { dstBuffer[i]=p[i]; } ; return;}
	if(Vxl.dtype == "MET_INT") { int32  *p = (int32 *) Vxl.data; for(size_t i =0; i<Vxl.N; i++) { dstBuffer[i]=p[i]; } ; return;}
	if(Vxl.dtype == "MET_LONG") { int32  *p = (int32 *) Vxl.data; for(size_t i =0; i<Vxl.N; i++) { dstBuffer[i]=p[i]; } ; return;}
	if(Vxl.dtype == "MET_UCHAR") { uint8  *p = (uint8  *)Vxl.data; for(size_t i =0; i<Vxl.N; i++) { dstBuffer[i]=p[i]; } ; return;}
	if(Vxl.dtype == "MET_USHORT") { uint16 *p = (uint16 *)Vxl.data; for(size_t i =0; i<Vxl.N; i++) { dstBuffer[i]=p[i]; } ; return;}
	if(Vxl.dtype == "MET_UINT") { uint32 *p = (uint32 *)Vxl.data; for(size_t i =0; i<Vxl.N; i++) { dstBuffer[i]=p[i]; } ; return;}
	if(Vxl.dtype == "MET_ULONG") { uint32 *p = (uint32 *)Vxl.data; for(size_t i =0; i<Vxl.N; i++) { dstBuffer[i]=p[i]; } ; return;}
	if(Vxl.dtype == "MET_FLOAT") { float32 *p = (float32 *)Vxl.data; for(size_t i =0; i<Vxl.N; i++) { dstBuffer[i]=p[i]; } ; return;}
	if(Vxl.dtype == "MET_DOUBLE") { float64 *p = (float64 *)Vxl.data; for(size_t i =0; i<Vxl.N; i++) { dstBuffer[i]=p[i]; } ; return;}

	logError(1,"map datatype not allowed");
	exit(1);
}
//=====================================================================
void exportBufferToVoxels(double *srcBuffer,Voxels &dstVxl){

	if(!srcBuffer) {logError(1,"NULL input buffer"); exit(1);}

	if(dstVxl.dtype == "MET_CHAR") { int8   *p = (int8  *) dstVxl.data; for(size_t i =0; i<dstVxl.N; i++) { p[i]=srcBuffer[i]; } ; return;}
	if(dstVxl.dtype == "MET_SHORT") { int16  *p = (int16 *) dstVxl.data; for(size_t i =0; i<dstVxl.N; i++) { p[i]=srcBuffer[i]; } ; return;}
	if(dstVxl.dtype == "MET_INT") { int32  *p = (int32 *) dstVxl.data; for(size_t i =0; i<dstVxl.N; i++) { p[i]=srcBuffer[i]; } ; return;}
	if(dstVxl.dtype == "MET_LONG") { int32  *p = (int32 *) dstVxl.data; for(size_t i =0; i<dstVxl.N; i++) { p[i]=srcBuffer[i]; } ; return;}
	if(dstVxl.dtype == "MET_UCHAR") { uint8  *p = (uint8  *)dstVxl.data; for(size_t i =0; i<dstVxl.N; i++) { p[i]=srcBuffer[i]; } ; return;}
	if(dstVxl.dtype == "MET_USHORT") { uint16 *p = (uint16 *)dstVxl.data; for(size_t i =0; i<dstVxl.N; i++) { p[i]=srcBuffer[i]; } ; return;}
	if(dstVxl.dtype == "MET_UINT") { uint32 *p = (uint32 *)dstVxl.data; for(size_t i =0; i<dstVxl.N; i++) { p[i]=srcBuffer[i]; } ; return;}
	if(dstVxl.dtype == "MET_ULONG") { uint32 *p = (uint32 *)dstVxl.data; for(size_t i =0; i<dstVxl.N; i++) { p[i]=srcBuffer[i]; } ; return;}
	if(dstVxl.dtype == "MET_FLOAT") { float32 *p = (float32 *)dstVxl.data; for(size_t i =0; i<dstVxl.N; i++) { p[i]=srcBuffer[i]; } ; return;}
	if(dstVxl.dtype == "MET_DOUBLE") { float64 *p = (float64 *)dstVxl.data; for(size_t i =0; i<dstVxl.N; i++) { p[i]=srcBuffer[i]; } ; return;}

	logError(1,"map datatype not allowed");
	exit(1);
	
}

//=====================================================================
double ArithmeticMean(double *map, int idx[3], int n[3]){
	const bool debug=false;
	int ix,iy,iz;
	ix=idx[0];iy=idx[1];iz=idx[2];

	int imin = ix==0 ? 0 : -1;
	int jmin = iy==0 ? 0 : -1;
	int kmin = iz==0 ? 0 : -1;

	int imax = ix==n[0]-1 ? 0 : +1;
	int jmax = iy==n[1]-1 ? 0 : +1;
	int kmax = iz==n[2]-1 ? 0 : +1;

	double sum = 0;
	int nval = 0;
	for(int i=ix+imin;i<=ix+imax;i++){
		for(int j=iy+jmin;j<=iy+jmax;j++){
			for(int k=iz+kmin;k<=iz+kmax;k++){
				double v = map[i+j*n[0]+k*n[0]*n[1]]; 
				sum+=v; nval++;
			}
		}
	}
	double mean = sum/nval;
	return mean;
	// return map[ix+iy*n[0]+iz*n[0]*n[1]]; 
}
//=====================================================================
double HarmonicMean(double *map, int idx[3], int n[3]){
	const bool debug=false;
	int ix,iy,iz;
	ix=idx[0];iy=idx[1];iz=idx[2];

	int imin = ix==0 ? 0 : -1;
	int jmin = iy==0 ? 0 : -1;
	int kmin = iz==0 ? 0 : -1;

	int imax = ix==n[0]-1 ? 0 : +1;
	int jmax = iy==n[1]-1 ? 0 : +1;
	int kmax = iz==n[2]-1 ? 0 : +1;

	double sum = 0;
	int nval = 0;
	for(int i=ix+imin;i<=ix+imax;i++){
		for(int j=iy+jmin;j<=iy+jmax;j++){
			for(int k=iz+kmin;k<=iz+kmax;k++){
				double v = map[i+j*n[0]+k*n[0]*n[1]]; 
				// sum+=v; nval++;
				if(v>0){
					sum+=1/v;
					nval++;
				}
			}
		}
	}
	// double mean = sum/nval;
	double hmean = nval>0? nval/sum : 0;
	return hmean;
}
//=====================================================================
double HotSpotCleaner(double *map, int idx[3], int n[3]){
	const bool debug=false;
	int ix,iy,iz;
	ix=idx[0];iy=idx[1];iz=idx[2];

	int imin = ix==0 ? 0 : -1;
	int jmin = iy==0 ? 0 : -1;
	int kmin = iz==0 ? 0 : -1;

	int imax = ix==n[0]-1 ? 0 : +1;
	int jmax = iy==n[1]-1 ? 0 : +1;
	int kmax = iz==n[2]-1 ? 0 : +1;

	double sum = 0;
	int nval = 0;
	for(int i=ix+imin;i<=ix+imax;i++){
		for(int j=iy+jmin;j<=iy+jmax;j++){
			for(int k=iz+kmin;k<=iz+kmax;k++){
				double v = map[i+j*n[0]+k*n[0]*n[1]]; 
				sum+=v; nval++;
			}
		}
	}
	double value = map[ix+iy*n[0]+iz*n[0]*n[1]]; 
	double surroundingAvg = (sum-value)/(nval-1);
	
	return value>hotspotThreshold*surroundingAvg? surroundingAvg : value;
	// return value>2*surroundingAvg? 1 : 0;
}

//=====================================================================
int main(int argc, char *argv[]){
	banner();
	for(int i=0;i<argc;i++) args.push_back(argv[i]);
  


  if(checkOption(args,"-h") || checkOption(args,"--help")) {
    usage();
    return 0;
  }

  if(checkOption(args,"-v") || checkOption(args,"--version")) {
  	cout<<version<<endl;
    return 0;
  }

  int iarg;
  if((iarg=checkOption(args,"-serial"))) {
	pthreads_max_num = 1;
	args.erase(args.begin()+iarg,args.begin()+iarg+1);
  }

  string outPath="smoothed";
	if((iarg=checkOption(args,"-o"))) {
    if(iarg+1>=args.size()) { 
    	logError(1,"output file is missing");
    	return 1;
    }
    outPath = args[iarg+1];
	args.erase(args.begin()+iarg,args.begin()+iarg+1+1);
  }

  if(args.size()==1){
    usage();
    return 1;
  }
  
  //================================     Reading maps    =========================================//

	voxelsA.read(args[1]);
  
	int   ierr;
	// ierr = map3d_read(args[1],mapA_struct);
	//cout<<"Loading original map3d : "<<(voxelsA.good? "OK": "ERROR")<<endl;
	if(voxelsA.good) {
		voxelsA.info();
	} else {
		exit(1);
	}

	voxelsB.clone(voxelsA);
	voxelsB.info();
	// exit(0);
  
	// convert to double precision
	importVoxelsToBuffer(voxelsA,&mapA);
	importVoxelsToBuffer(voxelsB,&mapB);

	
	size_t numNeg=0;
	for(size_t i=0;i<voxelsA.N;i++) if(mapA[i]<0) {
		mapA[i]=0;
		numNeg++;
	}
	if(numNeg>0){
	  //	cout<<"I've fixed negative voxels by setting them to zero"<<endl;
	  //	cout<<"Num of negative voxels: "<<numNeg<<" = "<<100.*numNeg/voxelsA.N<<"%"<<endl;
	}

	double sum=0;
	for(size_t i=0;i<voxelsA.N;i++) sum+=mapA[i];
	//cout<<"Sum of map A: "<<sum<<endl;

	// double positiveMin = 1e32;
	// for(size_t i=0;i<voxelsA.N;i++) if(mapA[i]>0 and positiveMin>mapA[i]) positiveMin=mapA[i];
	// cout<<"Minimum positive value: "<<positiveMin<<endl;
	// for(size_t i=0;i<voxelsA.N;i++) if(mapA[i]<positiveMin) mapA[i]=positiveMin;

	for(size_t i=0;i<voxelsB.N;i++) mapB[i]=0; // B = 0
	
	//	cout<<"Remove isolated hot spots: ";
	//cout<<"threshold >="<<hotspotThreshold<<" times";
	imode = 1;  // 1=hotspotcleaner
	evaluationLooper();
	//cout<<endl;
	for(size_t i=0;i<voxelsB.N;i++) mapA[i]=mapB[i]; // B => A

	//	cout<<"Apply arithmetic mean smoothing: ";
	imode = 2;  // 2=arithmetic mean
	evaluationLooper();
	//cout<<endl;


	sum=0;
	for(size_t i=0;i<voxelsB.N;i++) sum+=mapB[i];
	//cout<<"Sum of map B: "<<sum<<endl;
	//cout<<endl;

	exportBufferToVoxels(mapB,voxelsB);
	voxelsB.write(outPath,true);
	//cout<<"Smoothed map written to file "<<outPath<<endl;
	return 0;  
}

//================================================================================================

	void evaluationLooper(){
			
		// setup pthread environment
		vector<pthread_t> threads(pthreads_max_num);
		vector<threadInfo> tInfo(pthreads_max_num);
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		// setup thread info
		for (int i=0; i<pthreads_max_num; i++) {
			tInfo[i].thread_num = i;
			tInfo[i].num_threads = pthreads_max_num;
		}
		
		//  create threads
		for (int i=0; i<pthreads_max_num; i++) {
			pthread_create(&threads[i], &attr, kernelOnCPU, (void *) &tInfo[i]);
		}

		// wait for all threads to complete
		for (int i=0; i<pthreads_max_num; i++) {
			pthread_join(threads[i], NULL);
		}
	}


void *kernelOnCPU(void *info) 
{
	int tid = ((threadInfo *) info)->thread_num;
	int nthreads = ((threadInfo *) info)->num_threads;

	i3d nnA,nnB;
	vec3dRT hsA,hsB,x0A,x0B;

	getGrid(voxelsA,nnA,hsA,x0A);
	getGrid(voxelsB,nnB,hsB,x0B);

	int idx[3];

	for(size_t i=0;i<voxelsB.N;i++) {
		if(i%nthreads!=tid) continue;

		idx[2] = i/(nnB.x*nnB.y);
		idx[1] = (i%(nnB.x*nnB.y))/nnB.x;
		idx[0] = (i%(nnB.x*nnB.y))%nnB.x;

		switch(imode){
			case 1:
				mapB[i] = HotSpotCleaner(mapA, idx, nnB.v);
			break;
			case 2:
				mapB[i] = ArithmeticMean(mapA, idx, nnB.v);
			break;
			case 3:
				mapB[i] = HarmonicMean(mapA, idx, nnB.v);
			break;
			default:
			printf("undefined operation mode %d\n",imode); exit(1);
		}

	}

	pthread_exit(NULL);	
}

