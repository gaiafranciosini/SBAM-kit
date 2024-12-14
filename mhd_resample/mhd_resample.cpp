 /*
 *  resampling mhd files
 *  Fred project
 *
 *  Copyright 2011-2018 Università di Roma LA SAPIENZA. All rights reserved.
 *  A. Schiavi (2018)
 */

// Jun 2020 => porting from m3d to mhd
 
const char *version = "2.0";
/*			version history
2.0 : from m3d to mhd
1.6 : ghost voxel for boundary conditions; -noInterpolation option for sharp voxel edges
1.5 : extented allowed datatypes: int8,int16,int32,int64,,uint8,uint16,uint32,uint64,float32,float64
1.3 : resampling on arbitrary (even non overlapping) grids  (floating point only)
1.0 : initial version for Magda and Antoni (May 2018)
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
// inline int el3d(int i,int j,int k,int istride,int jstride,int kstride) {return i*istride+j*jstride+k*kstride;}
inline int el3d(int i,int j,int k,int nn[3]) {return i+j*nn[0]+k*nn[0]*nn[1];} // column-major and contiguous 

ifstream fin;	

vec3dRT xmin_A, xmin_B;
i3d nn_A, nn_B;
vec3dRT hs_A, hs_B;
vec3dRT xmax_A, xmax_B;

int maxlevel;

double * map_A;
double * map_B;

double TrilinearInterp(double *map, float position[3], float h[3], int n[3]);
double getVoxelValue  (double *map, float position[3], float h[3], int n[3]);
void evaluationLooper();

typedef struct threadInfo_s {
	int thread_num;
	int num_threads;
} threadInfo;

void *kernelOnCPU(void *info);
int pthreads_max_num = 8;

bool noInterpolation = false;
bool rebinning = false;
int newbin[3];

vector<const char *> args;



//=====================================================================
Voxels voxelsA,voxelsB;
double *mapA,*mapB;
float frameTransform[4][4];

//=====================================================================
void banner(){
	cout<<endl;
	cout<<"+------------------------------------------------------------------------------+"<<endl;
	cout<<"|                       Resample 3D map (Fred project)                         |"<<endl;
	cout<<"|                                                                              |"<<endl;
	cout<<"|                             A.Schiavi -  2020                                |"<<endl;
	cout<<"+------------------------------------------------------------------------------+"<<endl;
	cout<<endl;
	cout<<"                             version "<<version<<endl;
	cout<<endl;
}
//=====================================================================
void usage(){
	cout<<"usage: mapA mapB"<<endl;
	cout<<"resamples mapA on the grid of mapB"<<endl;
	cout<<"options: "<<endl;
	cout<<"\t-noInterpolation : do not interpolate on the reference mapA, use voxel value for all points inside the voxel"<<endl;
	cout<<"\t-rebin NX NY NZ : rebin mapA and resample it into mapB"<<endl;
	cout<<"\t-level lev [3]: level of grid oversampling: num of voxel subdivisions is level^3 [default = 27]"<<endl;
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
void getFrameTransform_B2A(Voxels &A,Voxels &B,float M[4][4])
{

	float G2A[4][4]   ={ // global to local A
		{A.left[0], A.left[1], A.left[2], -A.x0[0]}, 
		{A.up[0], A.up[1], A.up[2], -A.x0[1]}, 
		{A.front[0], A.front[1], A.front[2], -A.x0[2]}, 
		{0.0f ,  0.0f,  0.0f, 1.0f}
	};

	float B2G[4][4]   ={ // local B to global
		{B.left[0], B.up[0], B.front[0], +B.x0[0]}, 
		{B.left[1], B.up[1], B.front[1], +B.x0[1]}, 
		{B.left[2], B.up[2], B.front[2], +B.x0[2]}, 
		{0.0f ,  0.0f,  0.0f, 1.0f}
	};



	for(int i=0;i<4;i++) {
		for(int j=0;j<4;j++) {
			M[i][j]=0;
			for(int k=0;k<4;k++) {
				M[i][j]+=G2A[i][k]*B2G[k][j];
			}
		}
	}

}
//=====================================================================
void applyTransformToPoint(float M[4][4],vec3dRT &P)
{
	float homo[4],homo_new[4];
	homo[0]=P[0];homo[1]=P[1];homo[2]=P[2];homo[3]=1.0f;
	for(int i=0;i<4;i++){
		homo_new[i]=0.0f;
		for(int j=0;j<4;j++){
			homo_new[i]+=M[i][j]*homo[j];
		}
	}
	P[0]=homo_new[0];P[1]=homo_new[1];P[2]=homo_new[2];	
}

//=====================================================================
void addGhostVoxels(Voxels &Vxl,double **buffer){

	int32   nnold[3]={Vxl.nn[0],Vxl.nn[1],Vxl.nn[2]};
	int32   nnnew[3]={Vxl.nn[0]+2,Vxl.nn[1]+2,Vxl.nn[2]+2};
	float   x0new[3]={Vxl.x0[0]-Vxl.hs[0],Vxl.x0[1]-Vxl.hs[1],Vxl.x0[2]-Vxl.hs[2]};
	float   hsnew[3]={Vxl.hs[0],Vxl.hs[1],Vxl.hs[2]};

	size_t ghostN = 1UL*nnnew[0]*nnnew[1]*nnnew[2];
	double *p = new double[ghostN];
	for (int i = 0; i < ghostN; ++i) p[i]=0;

	// copy inner part (real map)
	for (int iz = 0; iz < Vxl.nn[2]; ++iz)
	{
		for (int iy = 0; iy < Vxl.nn[1]; ++iy)
		{
			for (int ix = 0; ix < Vxl.nn[0]; ++ix)
			{
				p[el3d(ix+1,iy+1,iz+1,nnnew)] = (*buffer)[el3d(ix,iy,iz,nnold)];
			}
		}
	}

	// set ghost voxels
	for (int iy = 1; iy < nnnew[1]-1; ++iy)
	{
		for (int ix = 1; ix < nnnew[0]-1; ++ix)
		{
			p[el3d(ix,iy,0,nnnew)] = (*buffer)[el3d(ix-1,iy-1,0,nnold)];
			p[el3d(ix,iy,nnnew[2]-1,nnnew)] = (*buffer)[el3d(ix-1,iy-1,nnold[2]-1,nnold)];
		}
	}
	for (int iz = 1; iz < nnnew[2]-1; ++iz)
	{
		for (int ix = 1; ix < nnnew[0]-1; ++ix)
		{
			p[el3d(ix,0,iz,nnnew)] = (*buffer)[el3d(ix-1,0,iz-1,nnold)];
			p[el3d(ix,nnnew[1]-1,iz,nnnew)] = (*buffer)[el3d(ix-1,nnold[1]-1,iz-1,nnold)];
		}
	}
	for (int iz = 1; iz < nnnew[2]-1; ++iz)
	{
		for (int iy = 1; iy < nnnew[1]-1; ++iy)
		{
			p[el3d(0,iy,iz,nnnew)] = (*buffer)[el3d(0,iy-1,iz-1,nnold)];
			p[el3d(nnnew[0]-1,iy,iz,nnnew)] = (*buffer)[el3d(nnold[0]-1,iy-1,iz-1,nnold)];
		}
	}

	*buffer = p;

	Vxl.nn[0]=nnnew[0];
	Vxl.nn[1]=nnnew[1];
	Vxl.nn[2]=nnnew[2];
	Vxl.x0[0]=x0new[0];
	Vxl.x0[1]=x0new[1];
	Vxl.x0[2]=x0new[2];

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

//============================== 3-linear interpolation ========================================//

double TrilinearInterp(double *map, float position[3], float h[3], int n[3]){
	const bool debug=false;
	float x_d, y_d, z_d;
	double c00,c10,c01,c11,c0,c1,val;
	int ix,iy,iz;
	float fx,fy,fz;

	// position in is local frame of reference

	// compute point logical coords
	fx=position[0]/h[0]; 
	fy=position[1]/h[1];
	fz=position[2]/h[2];

	if(fx<=0 || fy<=0 || fz<=0 ) return 0.0f; // OOB
	if(fx>=n[0] || fy>=n[1] || fz>=n[2]) return 0.0f; // OOB

	

	if( fx>=0.5f && fy>=0.5f && fz>=0.5f && 
		fx<=n[0]-0.5f &&  fy<=n[1]-0.5f && fz<=n[2]-0.5f) {

		// shift coords system from nodes to voxel centers
		fx -= 0.5f; fy -= 0.5f; fz -= 0.5f; 

		// compute voxel indices
		ix = fx;
		iy = fy;
		iz = fz;

		// check upper index range
		if (ix==n[0]-1) ix--;
		if (iy==n[1]-1) iy--;
		if (iz==n[2]-1) iz--;

		// compute normalized distance from voxel center (where the map value is known)
		x_d= fx - ix; 
		y_d= fy - iy; 
		z_d= fz - iz; 

		// interpolation
		c00 = map[ix+iy*n[0]+iz*n[0]*n[1]]*(1.0f-x_d) + map[(ix+1)+iy*n[0]+iz*n[0]*n[1]]*x_d;
		c10 = map[ix+(iy+1)*n[0]+iz*n[0]*n[1]]*(1.0f-x_d) + map[(ix+1)+(iy+1)*n[0]+iz*n[0]*n[1]]*x_d;
		c01 = map[ix+iy*n[0]+(iz+1)*n[0]*n[1]]*(1.0f-x_d) + map[(ix+1)+iy*n[0]+(iz+1)*n[0]*n[1]]*x_d;
		c11 = map[ix+(iy+1)*n[0]+(iz+1)*n[0]*n[1]]*(1.0f-x_d) + map[(ix+1)+(iy+1)*n[0]+(iz+1)*n[0]*n[1]]*x_d;

		c0 = c00*(1.0f-y_d) + c10*y_d;
		c1 = c01*(1.0f-y_d) + c11*y_d;

		val = (c0*(1.0f-z_d) + c1*z_d);

		if(debug) cout<<"trilininterp "<<endl;
		// if(debug) cout<<"x0  "<<x0[0]<<' '<<x0[1]<<' '<<x0[2]<<endl;
		if(debug) cout<<"h   "<< h[0]<<' '<< h[1]<<' '<< h[2]<<endl;
		if(debug) cout<<"n   "<< n[0]<<' '<< n[1]<<' '<< n[2]<<endl;
		if(debug) cout<<"pos "<<position[0]<<' '<<position[1]<<' '<<position[2]<<endl;
		if(debug) cout<<"fx  "<<fx<<' '<<fy<<' '<<fz<<endl;
		if(debug) cout<<"idx "<<ix<<' '<<iy<<' '<<iz<<endl;
		if(debug) cout<<"xd  "<<x_d<<' '<<y_d<<' '<<z_d<<" => "<<val<<endl<<endl;

		return val;
	} else{
		// boundary points => return non-interpolated dose value at the center of the voxel
		ix = fx;
		iy = fy;
		iz = fz;
		return map[ix+iy*n[0]+iz*n[0]*n[1]];
	}
}
//=====================================================================
double getVoxelValue(double *map, float position[3], float h[3], int n[3]){
	const bool debug=false;
	int ix,iy,iz;
	// position in is local frame of reference

	// compute voxel indices
	ix=floor(position[0]/h[0]);
	iy=floor(position[1]/h[1]);
	iz=floor(position[2]/h[2]);

	if(ix<0 || iy<0 || iz<0 ) return 0.0f; // OOB
	if(ix>=n[0] || iy>=n[1] || iz>=n[2]) return 0.0f; // OOB

	return map[ix+iy*n[0]+iz*n[0]*n[1]];
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
    return 0;
  }

  int iarg;
  if((iarg=checkOption(args,"-serial"))) {
	pthreads_max_num = 1;
	args.erase(args.begin()+iarg,args.begin()+iarg+1);
  }

  if((iarg=checkOption(args,"-noInterpolation"))) {
    noInterpolation = true;
	args.erase(args.begin()+iarg,args.begin()+iarg+1);
  }

  if((iarg=checkOption(args,"-rebin"))) {
    rebinning = true;
    if(iarg+3>=args.size()) { 
    	logError(1,"new bins are missing");
    	return 1;
    }
    newbin[0] = atoi(args[iarg+1]);
    newbin[1] = atoi(args[iarg+2]);
    newbin[2] = atoi(args[iarg+3]);
	args.erase(args.begin()+iarg,args.begin()+iarg+1+3);
  }

  maxlevel = 3;
	if((iarg=checkOption(args,"-level"))) {
    if(iarg+1>=args.size()) { 
    	logError(1,"subdivision level is missing");
    	return 1;
    }
    maxlevel = atoi(args[iarg+1]);
	args.erase(args.begin()+iarg,args.begin()+iarg+1+1);
  }

  // for(int i=0;i<args.size();i++) cout<<i<<' '<<args[i]<<endl;

  if(args.size()==1){
    usage();
    return 1;
  }
  
  if(rebinning){
  	if(args.size()<2){
	    logError(1,"original mhd file is missing");
	    return 1;  	
	}
  } else if(args.size()<2+1){
    logError(1,"mhd file are missing");
    return 1;
  }

  //================================     Reading maps    =========================================//

	voxelsA.read(args[1]);
  
	int   ierr;
	// ierr = map3d_read(args[1],mapA_struct);
	cout<<"Loading original map3d : "<<(voxelsA.good? "OK": "ERROR")<<endl;
	if(voxelsA.good) {
		voxelsA.info();
	} else {
		exit(1);
	}

	if(rebinning){
		voxelsB.clone(voxelsA);
		voxelsB.resize(newbin);
	} else {
		voxelsB.read(args[2]);
 		voxelsB.realloc(voxelsA.dtype,voxelsB.nn,voxelsB.hs,voxelsB.x0); // use grid of B, but datatype of A
	}
	voxelsB.info();

	getFrameTransform_B2A(voxelsA,voxelsB,frameTransform);
	for(int i=0;i<4;i++) {
		for(int j=0;j<4;j++) {
			cout<<frameTransform[i][j]<<' ';
		}
		cout<<endl;
	}
  
	// convert to double precision
	importVoxelsToBuffer(voxelsA,&mapA);
	importVoxelsToBuffer(voxelsB,&mapB);

	float VolA = voxelsA.hs[0]*voxelsA.hs[1]*voxelsA.hs[2];
	float VolB = voxelsB.hs[0]*voxelsB.hs[1]*voxelsB.hs[2];

	double sum=0;
	for(size_t i=0;i<voxelsA.N;i++) sum+=mapA[i];
	cout<<"Sum of map A: "<<sum<<endl;
	cout<<"Integral of map A: "<<sum*VolA<<endl;
	cout<<endl;
	for(size_t i=0;i<voxelsB.N;i++) mapB[i]=0;

	// for(size_t i=0;i<voxelsA.N;i++) cout<<i<<' '<<mapA[i]<<endl;

	addGhostVoxels(voxelsA,&mapA);

	getFrameTransform_B2A(voxelsA,voxelsB,frameTransform);
	for(int i=0;i<4;i++) {
		for(int j=0;j<4;j++) {
			cout<<frameTransform[i][j]<<' ';
		}
		cout<<endl;
	}

	if(noInterpolation){
		cout<<"Reference map sampling mode: use voxel value for all points inside the voxel (no interpolation)"<<endl;
	} else {
		cout<<"Reference map sampling mode: interpolate using neighboring voxels"<<endl;
	}
	cout<<endl;

	printf("Executing the kernel....\n");
  
	evaluationLooper();
	cout<<endl;

	sum=0;
	for(size_t i=0;i<voxelsB.N;i++) sum+=mapB[i];
	cout<<"Sum of map B: "<<sum<<endl;
	cout<<"Integral of map B: "<<sum*VolB<<endl;
	cout<<endl;

	// for(size_t i=0;i<voxelsB.N;i++) cout<<i<<' '<<mapA[i]<<' '<<mapB[i]<<endl;


	exportBufferToVoxels(mapB,voxelsB);
	if(rebinning){
		voxelsB.write("rebinned.mhd",false);	
	} else {
		voxelsB.write("resampled.mhd");	
	}
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

	i3d idx;
	vec3dRT nodeB,pos;


	int subdiv = maxlevel;
	int npt = subdiv*subdiv*subdiv;
	
	vector<vec3dRT> subpts(npt);
	
	int32 i,j,k;
	int32 n=0;
	double f = 1.0/subdiv, f2=f/2;
	for (i=0; i<subdiv; i++) {
		for (j=0; j<subdiv; j++) {
			for (k=0; k<subdiv; k++) {
				subpts[n].set(f2+i*f,f2+j*f,f2+k*f);
				n++;
			}
		}
	}

	if(tid==0){ 
		cout<<endl<<"Num of sub-voxels: "<<npt<<endl;
		cout<<"Intra-voxel sampling points: "<<endl;
		for (int i = 0; i < subpts.size(); ++i) cout<<i+1<<" : "<<subpts[i]<<endl;
		cout<<endl;
		cout<<"nnA,hsA,x0A: "<<nnA<<' '<<hsA<<' '<<x0A <<" ==> "<<x0A+hsA*nnA<<endl;
		cout<<"nnB,hsB,x0B: "<<nnB<<' '<<hsB<<' '<<x0B <<" ==> "<<x0B+hsB*nnB<<endl;
	}

	for (int i = 0; i < subpts.size(); ++i) subpts[i] *= hsB;
	
	for(size_t i=0;i<voxelsB.N;i++) {
		if(i%nthreads!=tid) continue;

		int iz = i/(nnB.x*nnB.y);
		int iy = (i%(nnB.x*nnB.y))/nnB.x;
		int ix = (i%(nnB.x*nnB.y))%nnB.x;

		idx.set(ix,iy,iz);

		nodeB = hsB*idx;
		mapB[i]=0;
		for (int j = 0; j < subpts.size(); ++j) {
			pos = nodeB + subpts[j];
			// pos is local frame of B => go to local frame of A
			// cout<<pos<<"  => ";
			applyTransformToPoint(frameTransform,pos);
			// cout<<pos<<endl;

			float valueMapA; 
			if(noInterpolation) {
				valueMapA = getVoxelValue(mapA, pos.v, hsA.v, nnA.v);
			} else {
				valueMapA = TrilinearInterp(mapA, pos.v, hsA.v, nnA.v); 
			}
			mapB[i] += valueMapA;
			// cout<<i<<' '<<valueMapA<<" -> "<<mapB[i]<<endl;
		}
		mapB[i] /= npt;
	}

	pthread_exit(NULL);	
}

