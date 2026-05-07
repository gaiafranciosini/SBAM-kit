 /*
 *  
 *  Fred project
 *
 *  
 *  A. Schiavi (2020)
 */

#include "Voxels.h"
#include "mhd_io.h"
#include <cstring>

Voxels::Voxels() : N(_N)
{ 
	nn[0]=nn[1]=nn[2]=-1;
	hs[0]=hs[1]=hs[2]=0;
	L[0]=L[1]=L[2]=0;
	dtype="UNKNOWN";
	good=false;
	data=nullptr;
	_N=0;
}

Voxels::Voxels(string _dtype,int _nn[3],float _hs[3], float _x0[3]): N(_N)
{
	data = nullptr;
	realloc(_dtype,_nn,_hs,_x0);
}

void Voxels::realloc(string _dtype,int _nn[3],float _hs[3], float _x0[3])
{
	dtype=_dtype;
	for(int i=0;i<3;i++) nn[i]=_nn[i];
	for(int i=0;i<3;i++) x0[i]=_x0[i];
	for(int i=0;i<3;i++) hs[i]=_hs[i];

	hs2L();
	_N = 1UL*nn[0]*nn[1]*nn[2];
	eye();


	size_t sz = dataSize();
	if(data!=nullptr) delete [] data; // release memory
	data = new char[sz];

	good= data!=nullptr;
}


Voxels::~Voxels(){
	if(data!=nullptr) delete [] data;
}

void Voxels::eye(){
	float _left[3] ={1,0,0};
	float _up[3]   ={0,1,0};
	float _front[3]={0,0,1};
	for(int i=0;i<3;i++) left[i]=_left[i];
	for(int i=0;i<3;i++) up[i]=_up[i];
	for(int i=0;i<3;i++) front[i]=_front[i];
}

int Voxels::read(string path){
	size_t pos;
	string dataFile;
	int ierr = read_mhd_header(path,dtype,nn,x0,L,left,up,front,dataFile,pos);
	if (!ierr){
		ierr = read_mhd_data(dataFile,dtype,nn,pos,&data);
	}
	if (!ierr){
		L2hs();
		_N = 1UL*nn[0]*nn[1]*nn[2];
		good = true;
	}
	return ierr;
}
int Voxels::write(string path,bool localData){
	return write_mhd(path,dtype,data,nn,x0,L,left,up,front,localData);
}

void Voxels::info(ostream &os){
	os<<"dtype = "<<dtype<<endl;
	os<<"nn = "<<nn[0]<<' '<<nn[1]<<' '<<nn[2]<<endl;
	os<<"hs = "<<hs[0]<<' '<<hs[1]<<' '<<hs[2]<<endl;
	os<<"x0 = "<<x0[0]<<' '<<x0[1]<<' '<<x0[2]<<endl;
	os<<"L  = "<<L[0]<<' '<<L[1]<<' '<<L[2]<<endl;
	os<<"left  = "<<left[0]<<' '<<left[1]<<' '<<left[2]<<endl;
	os<<"up    = "<<up[0]<<' '<<up[1]<<' '<<up[2]<<endl;
	os<<"front = "<<front[0]<<' '<<front[1]<<' '<<front[2]<<endl;
}

size_t Voxels::dataSize(){
	size_t sz=0;
	if(dtype == "MET_CHAR") { sz = 1;}
	if(dtype == "MET_SHORT") { sz = 2;}
	if(dtype == "MET_INT") { sz = 4;}
	if(dtype == "MET_LONG") { sz = 4;}
	if(dtype == "MET_UCHAR") { sz = 1;}
	if(dtype == "MET_USHORT") { sz = 2;}
	if(dtype == "MET_UINT") { sz = 4;}
	if(dtype == "MET_ULONG") { sz = 4;}
	if(dtype == "MET_FLOAT") { sz = 4;}
	if(dtype == "MET_DOUBLE") { sz = 8;}
	return sz*N;
}

void Voxels::clone(Voxels &B){
	dtype=B.dtype;
	for(int i=0;i<3;i++) nn[i]=B.nn[i];
	for(int i=0;i<3;i++) x0[i]=B.x0[i];
	for(int i=0;i<3;i++) hs[i]=B.hs[i];
	for(int i=0;i<3;i++) L[i]=B.L[i];
	for(int i=0;i<3;i++) left[i]=B.left[i];
	for(int i=0;i<3;i++) up[i]=B.up[i];
	for(int i=0;i<3;i++) front[i]=B.front[i];
	_N=B.N;
	good=B.good;

	size_t sz = dataSize();
	if(data!=nullptr) delete [] data;
	data = new char[sz];
	memcpy(data,B.data,sz);
}

void Voxels::resize(int nnnew[3]){
	if(data!=nullptr) delete [] data;
	for(int i=0;i<3;i++) nn[i]=nnnew[i];
	L2hs();
	_N = 1UL*nn[0]*nn[1]*nn[2];
	size_t sz = dataSize();
	data = new char[sz];
	memset(data,0,sz);
}

void Voxels::resize(int nx,int ny,int nz){
	int nnnew[3]={nx,ny,nz};
	resize(nnnew);
}
