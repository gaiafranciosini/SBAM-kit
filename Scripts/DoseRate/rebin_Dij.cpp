/*
 *  rebin Dij binary files
 *  Fred project
 *
 *  Copyright 2011-2018 Università di Roma LA SAPIENZA. All rights reserved.
 *  A. Schiavi (2018)
 */

// Jun 2020
 
const char *version = "1.0";
/*			version history
1.0 : initial version for Marta (Jun 2020) based on previous resampling routine
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
int nnnew[3];
float hsnew[3];

vector<const char *> args;



//=====================================================================
Voxels voxelsA,voxelsB;
Voxels voxelsA_orig;
double *mapA,*mapB;
float frameTransform[4][4];

//=====================================================================
void banner(){
	cout<<endl;
	cout<<"+------------------------------------------------------------------------------+"<<endl;
	cout<<"|                  Rebin dose influence matrix Dij (Fred project)              |"<<endl;
	cout<<"|                                                                              |"<<endl;
	cout<<"|                             A.Schiavi -  2020                                |"<<endl;
	cout<<"+------------------------------------------------------------------------------+"<<endl;
	cout<<endl;
	cout<<"                             version "<<version<<endl;
	cout<<endl;
}
//=====================================================================
void usage(){
	cout<<"usage: Dij.bin NX NY NZ"<<endl;
	cout<<"rebin original Dij into Dijrebin with dimensions NX NY NZ"<<endl;
	cout<<"options: "<<endl;
	cout<<"\t-noInterpolation : do not interpolate on the reference mapA, use voxel value for all points inside the voxel"<<endl;
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

	int32   nnghost[3]={Vxl.nn[0]+2,Vxl.nn[1]+2,Vxl.nn[2]+2};
	int32   nnreal[3]={Vxl.nn[0],Vxl.nn[1],Vxl.nn[2]};

	float   x0new[3]={Vxl.x0[0]-Vxl.hs[0],Vxl.x0[1]-Vxl.hs[1],Vxl.x0[2]-Vxl.hs[2]};
	float   hsnew[3]={Vxl.hs[0],Vxl.hs[1],Vxl.hs[2]};

	size_t ghostN = 1UL*nnghost[0]*nnghost[1]*nnghost[2];
	double *p = new double[ghostN];
	for (int i = 0; i < ghostN; ++i) p[i]=0;

	// copy inner part (real map)
	for (int iz = 0; iz < Vxl.nn[2]; ++iz)
	{
		for (int iy = 0; iy < Vxl.nn[1]; ++iy)
		{
			for (int ix = 0; ix < Vxl.nn[0]; ++ix)
			{
				p[el3d(ix+1,iy+1,iz+1,nnghost)] = (*buffer)[el3d(ix,iy,iz,nnreal)];
			}
		}
	}

	// set ghost voxels
	for (int iy = 1; iy < nnghost[1]-1; ++iy)
	{
		for (int ix = 1; ix < nnghost[0]-1; ++ix)
		{
			p[el3d(ix,iy,0,nnghost)] = (*buffer)[el3d(ix-1,iy-1,0,nnreal)];
			p[el3d(ix,iy,nnghost[2]-1,nnghost)] = (*buffer)[el3d(ix-1,iy-1,nnreal[2]-1,nnreal)];
		}
	}
	for (int iz = 1; iz < nnghost[2]-1; ++iz)
	{
		for (int ix = 1; ix < nnghost[0]-1; ++ix)
		{
			p[el3d(ix,0,iz,nnghost)] = (*buffer)[el3d(ix-1,0,iz-1,nnreal)];
			p[el3d(ix,nnghost[1]-1,iz,nnghost)] = (*buffer)[el3d(ix-1,nnreal[1]-1,iz-1,nnreal)];
		}
	}
	for (int iz = 1; iz < nnghost[2]-1; ++iz)
	{
		for (int iy = 1; iy < nnghost[1]-1; ++iy)
		{
			p[el3d(0,iy,iz,nnghost)] = (*buffer)[el3d(0,iy-1,iz-1,nnreal)];
			p[el3d(nnghost[0]-1,iy,iz,nnghost)] = (*buffer)[el3d(nnreal[0]-1,iy-1,iz-1,nnreal)];
		}
	}

	*buffer = p;

	Vxl.nn[0]=nnghost[0];
	Vxl.nn[1]=nnghost[1];
	Vxl.nn[2]=nnghost[2];
	Vxl.x0[0]=x0new[0];
	Vxl.x0[1]=x0new[1];
	Vxl.x0[2]=x0new[2];

}
//=====================================================================
void fillBufferWithGhostVoxels(Voxels &Vxl,double *bufferGhost,double *bufferReal){

	int32   nnghost[3]={Vxl.nn[0],Vxl.nn[1],Vxl.nn[2]};
	int32   nnreal[3]={Vxl.nn[0]-2,Vxl.nn[1]-2,Vxl.nn[2]-2};

	// copy inner part (real map)
	for (int iz = 0; iz < nnreal[2]; ++iz)
	{
		for (int iy = 0; iy < nnreal[1]; ++iy)
		{
			for (int ix = 0; ix < nnreal[0]; ++ix)
			{
				bufferGhost[el3d(ix+1,iy+1,iz+1,nnghost)] = bufferReal[el3d(ix,iy,iz,nnreal)];
			}
		}
	}

	// set ghost voxels
	for (int iy = 1; iy < nnghost[1]-1; ++iy)
	{
		for (int ix = 1; ix < nnghost[0]-1; ++ix)
		{
			bufferGhost[el3d(ix,iy,0,nnghost)] = bufferReal[el3d(ix-1,iy-1,0,nnreal)];
			bufferGhost[el3d(ix,iy,nnghost[2]-1,nnghost)] = bufferReal[el3d(ix-1,iy-1,nnreal[2]-1,nnreal)];
		}
	}
	for (int iz = 1; iz < nnghost[2]-1; ++iz)
	{
		for (int ix = 1; ix < nnghost[0]-1; ++ix)
		{
			bufferGhost[el3d(ix,0,iz,nnghost)] = bufferReal[el3d(ix-1,0,iz-1,nnreal)];
			bufferGhost[el3d(ix,nnghost[1]-1,iz,nnghost)] = bufferReal[el3d(ix-1,nnreal[1]-1,iz-1,nnreal)];
		}
	}
	for (int iz = 1; iz < nnghost[2]-1; ++iz)
	{
		for (int iy = 1; iy < nnghost[1]-1; ++iy)
		{
			bufferGhost[el3d(0,iy,iz,nnghost)] = bufferReal[el3d(0,iy-1,iz-1,nnreal)];
			bufferGhost[el3d(nnghost[0]-1,iy,iz,nnghost)] = bufferReal[el3d(nnreal[0]-1,iy-1,iz-1,nnreal)];
		}
	}
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

  if(args.size()!=1+4){
    usage();
    return 1;
  }
  

    nnnew[0] = atoi(args[2]);
    nnnew[1] = atoi(args[3]);
    nnnew[2] = atoi(args[4]);


    ifstream fin(args[1],ios::in|ios::binary);
    if(!fin){cerr<<"Error: cannot open "<<args[1]<<" for input"<<endl;exit(1);}

    ofstream fout("Dijrebin.bin",ios::out|ios::binary);
    if(!fout){cerr<<"Error: cannot open "<<"Dijrebin.bin"<<" for output"<<endl;exit(1);}


    int nn[3];
    float hs[3];
    float x0[3];
    int npb;

    fin.read((char *)nn,3*sizeof(int));
    fin.read((char *)hs,3*sizeof(float));
    fin.read((char *)x0,3*sizeof(float));
    fin.read((char *)&npb,1*sizeof(int));

    hsnew[0]=hs[0]*nn[0]/nnnew[0];
    hsnew[1]=hs[1]*nn[1]/nnnew[1];
    hsnew[2]=hs[2]*nn[2]/nnnew[2];

    fout.write((char *)nnnew,3*sizeof(int));
    fout.write((char *)hsnew,3*sizeof(float));
    fout.write((char *)x0,3*sizeof(float));
    fout.write((char *)&npb,1*sizeof(int));

    voxelsA.realloc("MET_FLOAT",nn,hs,x0);
    cout<<"Original map:"<<endl;
	voxelsA.info();
	cout<<endl;
	int numVoxelA=voxelsA.N;

    cout<<"Rebinned map:"<<endl;
    voxelsB.realloc("MET_FLOAT",nnnew,hsnew,x0);
	voxelsB.info();
	cout<<endl;

	double *mapAreal = new double[voxelsA.N];
	mapA = new double[voxelsA.N];
	mapB = new double[voxelsB.N];


	// rebin
	addGhostVoxels(voxelsA,&mapA); 

	getFrameTransform_B2A(voxelsA,voxelsB,frameTransform);

	// exit(0);

	cout<<"Num of spots: "<<npb<<endl;

// loop on the stored spots
	int jpb,ifield,ni;
	vector<int> Vi;
	vector<float> di;

	for(int j=0;j<npb;j++){
		// read Vi and di
		fin.read((char *)&jpb,sizeof(int));
	    if(!fin){cerr<<"IO Error at line __LINE__"<<endl;exit(1);}
		fin.read((char *)&ifield,sizeof(int));
	    if(!fin){cerr<<"IO Error at line __LINE__"<<endl;exit(1);}
		fin.read((char *)&ni,sizeof(int));
	    if(!fin){cerr<<"IO Error at line __LINE__"<<endl;exit(1);}
		Vi.resize(ni);
		di.resize(ni);
		cout<<"\t"<<jpb<<' '<<ifield<<" : "<<ni<<endl;
		if(ni>0){
			fin.read((char *)Vi.data(),ni*sizeof(int)  );
		    if(!fin){cerr<<"IO Error at line __LINE__"<<endl;exit(1);}
			fin.read((char *)di.data(),ni*sizeof(float));
		    if(!fin){cerr<<"IO Error at line __LINE__"<<endl;exit(1);}
		}

		// fill dense map A
		for(size_t i=0;i<numVoxelA;i++) mapAreal[i]=0; // reset
		for(size_t i=0;i<ni;i++) mapAreal[Vi[i]]=di[i]; // fill

		// rebin
		fillBufferWithGhostVoxels(voxelsA,mapA,mapAreal);
		printf("Rebinning ...\n");
		evaluationLooper();
		// zero-suppression on map B
		ni=0;
		for(size_t i=0;i<voxelsB.N;i++) if(mapB[i]!=0) ni++;
		Vi.resize(ni);
		di.resize(ni);
		ni=0;
		for(size_t i=0;i<voxelsB.N;i++) if(mapB[i]!=0)
		{
			Vi[ni]=i;
			di[ni]=mapB[i];
			ni++;
		}
		cout<<" ni = "<<ni<<endl;
		// store to new Dij matrix
		fout.write((char *)&jpb,sizeof(int));
		fout.write((char *)&ifield,sizeof(int));
		fout.write((char *)&ni,sizeof(int));
		if(ni>0){
			fout.write((char *)Vi.data(),ni*sizeof(int)  );
			fout.write((char *)di.data(),ni*sizeof(float));
		}
	}
	cout<<"Done: output written to Dijrebin.bin"<<endl;
 
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

	// if(tid==0){ 
	// 	cout<<endl<<"Num of sub-voxels: "<<npt<<endl;
	// 	cout<<"Intra-voxel sampling points: "<<endl;
	// 	for (int i = 0; i < subpts.size(); ++i) cout<<i+1<<" : "<<subpts[i]<<endl;
	// 	cout<<endl;
	// 	cout<<"nnA,hsA,x0A: "<<nnA<<' '<<hsA<<' '<<x0A <<" ==> "<<x0A+hsA*nnA<<endl;
	// 	cout<<"nnB,hsB,x0B: "<<nnB<<' '<<hsB<<' '<<x0B <<" ==> "<<x0B+hsB*nnB<<endl;
	// }

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

