/*
 *  mhd_io.h
 *
 *  Fred project
 *
 *  Created by A.Schiavi 
 *  Copyright 2020. All rights reserved.
 *
 */


// purpose: minimal IO routines for mhd file format

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

// explicit data types -->
typedef signed char           int8;       // 8 bit signed
typedef unsigned char         uint8;      // 8 bit unsigned
typedef short                 int16;      // 16 bit signed
typedef unsigned short        uint16;     // 16 bit unsigned
typedef int                   int32;      // 32 bit signed
typedef unsigned int          uint32;     // 32 bit unsigned
typedef long long             int64;      // 64 bit signed
typedef unsigned long long    uint64;     // 64 bit unsigned

typedef float                 float32;     // 32 bit floating point 
typedef double                float64;     // 64 bit floating point
// explicit data types <--

template<typename T>
string dataTypeString(){
	string dtype = "UNKNOWN";
	if (typeid(T)==typeid(int8))    dtype = "MET_CHAR";
	if (typeid(T)==typeid(int16))   dtype = "MET_SHORT";
	if (typeid(T)==typeid(int32))   dtype = "MET_INT";
	if (typeid(T)==typeid(int64))   dtype = "MET_LONG";
	if (typeid(T)==typeid(uint8))   dtype = "MET_UCHAR";
	if (typeid(T)==typeid(uint16))  dtype = "MET_USHORT";
	if (typeid(T)==typeid(uint32))  dtype = "MET_UINT";
	if (typeid(T)==typeid(uint64))  dtype = "MET_ULONG";
	if (typeid(T)==typeid(float32)) dtype = "MET_FLOAT";
	if (typeid(T)==typeid(float64)) dtype = "MET_DOUBLE";
	return dtype;
}

int write_mhd(string fname,string &dtype,char *data,int nn[3],float x0[3],float L[3],float left[3],float up[3],float front[3],bool localData=true);

template<typename T>
void write_mhd(string fname,vector<T> &A,int nn[3],float x0[3],float L[3],float left[3],float up[3],float front[3],bool localData=true){
	string dtype = dataTypeString<T>();
	write_mhd(fname,dtype,(char *)A.data(),nn,x0,L,left,up,front,localData);
}



string tolower(const string &str);
int icompare(std::string const& a, std::string const& b);

int read_mhd_header(string fname,string &dtype,int nn[3],float x0[3],float L[3],float left[3],float up[3],float front[3],string &dataFile, size_t &dataFileOffset);
int read_mhd_data(string fname,string &dtype,int nn[3],size_t &fileOffset,char **data);

template<typename T>
int read_mhd(string fname,vector<T> &A,int nn[3],float x0[3],float L[3],float left[3],float up[3],float front[3])
{
	string dtype_src;
	size_t dataFileOffset;
	string dataFile;
	int ierr = read_mhd_header(fname,dtype_src,nn,x0,L,left,up,front,dataFile,dataFileOffset);
	if(ierr) return ierr;
	string dtype_dst = dataTypeString<T>();
    if(icompare(dtype_src,dtype_dst)){cerr<<"Error: datatypes do not match"<<endl;return 1;}

    ifstream fin;
	fin.open(dataFile,ios_base::in|ios_base::binary);
	fin.seekg(dataFileOffset,fin.beg);
    size_t numVxl = 1UL*nn[0]*nn[1]*nn[2];
	A.resize(numVxl);
	fin.read((char  *)A.data(),sizeof(T)*A.size());
	if(!fin) return 1;
	fin.close();

	return 0;
}
