/*
 *  mhd_io.cpp
 *
 *  Fred project
 *
 *  Created by A.Schiavi 
 *  Copyright 2020. All rights reserved.
 *
 */


// purpose: minimal IO routines for mhd file format

#include "mhd_io.h"
#include "InputTools.h"
#include <iostream>
using namespace std;

string tolower(const string &str){
    string s = str;
    for(int i=0;i<s.length();i++) s[i] = ::tolower(s[i]);
    return s;
}

int icompare(std::string const& a, std::string const& b)
{
    if (a.length()<b.length()) return -1;
    if (a.length()>b.length()) return +1;
    for(int i=0;i<a.length();i++) if(std::tolower(a[i])!=std::tolower(b[i])) return std::tolower(a[i])<std::tolower(b[i]) ? -1 : +1;
    return 0;
}


int read_mhd_header(string fname,string &dtype,int nn[3],float x0[3],float L[3],float left[3],float up[3],float front[3],string &dataFile,size_t &dataFileOffset)
{
	const bool debug = false;

	// set default transform matrix = Identity
	for(int i=0;i<3;i++) {left[i]=up[i]=front[i]=0;}
	left[0]=up[1]=front[2]=1; // eye


    ifstream fin(fname);
    if(!fin){cerr<<"Error: could not open "<<fname<<endl; exit(-1);}
    string tok;
    int NDims;
    float hs[3],offset[3];
    while(1){
        if(fin.eof()) break;
        fin>>tok;
        if(!fin) break; // error
        if(icompare(tok,"ObjectType")==0){
            fin>>tok;
            if(tok!="=") {cerr<<"Error in reading mhd header"<<endl; return 1;}
            fin>>tok;
            if(icompare(tok,"Image")) {cerr<<"Error in reading mhd header"<<endl; return 1;}
            continue;
        }
        if(icompare(tok,"NDims")==0){
            fin>>tok;
            if(tok!="=") {cerr<<"Error in reading mhd header"<<endl; return 1;}
            fin>>NDims;
            if(debug) cout<<"Ndims = "<<NDims<<endl;
            continue;
        }
        if(icompare(tok,"DimSize")==0){
            fin>>tok;
            if(tok!="=") {cerr<<"Error in reading mhd header"<<endl; return 1;}
            if(debug) cout<<"nn = ";
            for(int i=0;i<NDims;i++) {fin>>nn[i];if(debug) cout<<nn[i]<<' ';}
            if(debug) cout<<endl;
            continue;
        }
        if(icompare(tok,"ElementType")==0){
            fin>>tok;
            if(tok!="=") {cerr<<"Error in reading mhd header"<<endl; return 1;}
            fin>>dtype;
            if(debug) cout<<"ElementType = "<<dtype<<endl;
            continue;
        }
        if(icompare(tok,"ElementSpacing")==0){
            fin>>tok;
            if(tok!="=") {cerr<<"Error in reading mhd header"<<endl; return 1;}
            if(debug) cout<<"hs = ";
            for(int i=0;i<NDims;i++) {fin>>hs[i];if(debug) cout<<hs[i]<<' ';}
            if(debug) cout<<endl;
            continue;
        }
        if(icompare(tok,"Offset")==0){
            fin>>tok;
            if(tok!="=") {cerr<<"Error in reading mhd header"<<endl; return 1;}
            if(debug) cout<<"offset = ";
            for(int i=0;i<NDims;i++) {fin>>offset[i];if(debug) cout<<offset[i]<<' ';}
            if(debug) cout<<endl;
            continue;
        }
        if(icompare(tok,"TransformMatrix")==0){
            fin>>tok;
            if(tok!="=") {cerr<<"Error in reading mhd header"<<endl; return 1;}
            if(debug) cout<<"TransformMatrix = ";
            for(int i=0;i<NDims;i++) {fin>>left[i];if(debug) cout<<left[i]<<' ';}
            for(int i=0;i<NDims;i++) {fin>>up[i];if(debug) cout<<up[i]<<' ';}
            for(int i=0;i<NDims;i++) {fin>>front[i];if(debug) cout<<front[i]<<' ';}
            if(debug) cout<<endl;
            continue;
        }
        if(icompare(tok,"BinaryData")==0){
            fin>>tok;
            if(tok!="=") {cerr<<"Error in reading mhd header"<<endl; return 1;}
            fin>>tok;
            if(icompare(tok,"true")){cerr<<"Error: expected BinaryData = True"<<endl;return 1;}
            continue;
        }
        if(icompare(tok,"BinaryDataByteOrderMSB")==0){
            fin>>tok;
            if(tok!="=") {cerr<<"Error in reading mhd header"<<endl; return 1;}
            fin>>tok;
            if(icompare(tok,"false")){cerr<<"Error: expected BinaryDataByteOrderMSB = False"<<endl;return 1;}
            continue;
        }
        if(icompare(tok,"CompressedData")==0){
            fin>>tok;
            if(tok!="=") {cerr<<"Error in reading mhd header"<<endl; return 1;}
            fin>>tok;
            if(icompare(tok,"false")){cerr<<"Error: expected CompressedData = False"<<endl;return 1;}
            continue;
        }
        if(icompare(tok,"ElementDataFile")==0){
            fin>>tok;
            if(tok!="=") {cerr<<"Error in reading mhd header"<<endl; return 1;}
            fin>>tok;
            if(icompare(tok,"LOCAL")){
                // get current directory and substitute header file with data file
                vector<string> pathtokens = strtokens(fname,"/");
                // cout<<"fname "<<fname<<endl;
                // for(int i=0;i<pathtokens.size();i++) cout<<i<<" "<<pathtokens[i]<<endl;
                dataFile="";
                if(fname[0]=='/') dataFile="/";
                for(int i=0;i<pathtokens.size()-1;i++) dataFile+=pathtokens[i]+"/";
                dataFile += tok;
                // cout<<"ElementDataFile "<<dataFile<<endl;
                // cout<<endl;
                dataFileOffset=0;
            } else {
                dataFile=fname; // is LOCAL => same file
                if(fin.peek()=='\r') fin.seekg(1,fin.cur);
                if(fin.peek()=='\n') fin.seekg(1,fin.cur);
                dataFileOffset = fin.tellg();                
            }
            break;
        }
        // cout<<tok<<endl;
        // if(tok.rfind("ElementDataFile")!=string::npos) break;
    }
    if(debug) cout<<"data offset = "<<dataFileOffset<<endl;
    fin.close();

	// from mhd to fred 
	x0[0]=offset[0]-0.5*hs[0]*left[0]-0.5*hs[1]*up[0]-0.5*hs[2]*front[0];
	x0[1]=offset[1]-0.5*hs[0]*left[1]-0.5*hs[1]*up[1]-0.5*hs[2]*front[1];
	x0[2]=offset[2]-0.5*hs[0]*left[2]-0.5*hs[1]*up[2]-0.5*hs[2]*front[2];

    if(debug) cout<<"x0 = "<<x0[0]<<' '<<x0[1]<<' '<<x0[2]<<endl;

	L[0]=nn[0]*hs[0];
	L[1]=nn[1]*hs[1];
	L[2]=nn[2]*hs[2];

	return 0;
}

int read_mhd_data(string fname,string &dtype,int nn[3],size_t &fileOffset,char **data)
{
    ifstream fin;
    fin.open(fname,ios_base::in|ios_base::binary);
    fin.seekg(fileOffset,fin.beg);
    size_t numVxl = 1UL*nn[0]*nn[1]*nn[2];
    size_t buffSize;
    if(dtype == "MET_CHAR") { buffSize = 1;}
    if(dtype == "MET_SHORT") { buffSize = 2;}
    if(dtype == "MET_INT") { buffSize = 4;}
    if(dtype == "MET_LONG") { buffSize = 4;}
    if(dtype == "MET_UCHAR") { buffSize = 1;}
    if(dtype == "MET_USHORT") { buffSize = 2;}
    if(dtype == "MET_UINT") { buffSize = 4;}
    if(dtype == "MET_ULONG") { buffSize = 4;}
    if(dtype == "MET_FLOAT") { buffSize = 4;}
    if(dtype == "MET_DOUBLE") { buffSize = 8;}
    buffSize *= numVxl;
    // cout<<"buffSize "<<buffSize<<endl;
    // cout<<"fileOffset "<<fileOffset<<endl;
    char *buff = new char[buffSize];
    fin.read(buff,buffSize);
    if(!fin) return 1;
    fin.close();
    *data = buff;
    return 0;
}


 
int write_mhd(string fname,string &dtype,char *data,int nn[3],float x0[3],float L[3],float left[3],float up[3],float front[3], bool localData){

    // check if extension is present
    string basename;
    if(fname.size()<4) {basename = fname;}
    else if(fname.substr(fname.size()-4)==".mhd"){
        basename=fname.substr(0,fname.size()-4);
    } else{
        basename=fname;
    }
    // cout<<"fname "<<fname<<endl;
    // cout<<"basename "<<basename<<endl;
    // cout<<"mhd "<<basename+".mhd"<<endl;
    // cout<<"raw "<<basename+".raw"<<endl;
    // exit(0);
    fname=basename+".mhd";
    ofstream fout(fname);
    if(!fout) return 1;

    fout<<"ObjectType = Image"<<endl;
    fout<<"NDims = 3"<<endl;
    fout<<"DimSize = "<<nn[0]<<' '<<nn[1]<<' '<<nn[2]<<endl;
    fout<<"BinaryData = True"<<endl;
    fout<<"BinaryDataByteOrderMSB = False"<<endl;
    fout<<"CompressedData = False"<<endl;
    // fout<<"AnatomicalOrientation = RAI"<<endl;
    
    fout<<"TransformMatrix = ";
    fout<<left [0]<<' '<<left [1]<<' '<<left [2]<<' ';
    fout<<up   [0]<<' '<<up   [1]<<' '<<up   [2]<<' ';
    fout<<front[0]<<' '<<front[1]<<' '<<front[2]<<' ';
    fout<<endl;
    
    float spacing[3],delta[3],offset[3];
    spacing[0]=L[0]/nn[0];
    spacing[1]=L[1]/nn[1];
    spacing[2]=L[2]/nn[2];

    delta[0]=0.5f*spacing[0];
    delta[1]=0.5f*spacing[1];
    delta[2]=0.5f*spacing[2];

    offset[0]=x0[0]+delta[0]*left[0]+delta[1]*up[0]+delta[2]*front[0];
    offset[1]=x0[1]+delta[0]*left[1]+delta[1]*up[1]+delta[2]*front[1];
    offset[2]=x0[2]+delta[0]*left[2]+delta[1]*up[2]+delta[2]*front[2];  
    
    fout<<"Offset = "<<offset[0]<<' '<<offset[1]<<' '<<offset[2]<<endl;
    fout<<"ElementSpacing = "<<spacing[0]<<' '<<spacing[1]<<' '<<spacing[2]<<endl;
    fout<<"ElementType = "<<dtype<<endl;
    if(localData){
        fout<<"ElementDataFile = LOCAL"<<endl;
    } else {
        fout<<"ElementDataFile = "<<basename+".raw"<<endl;
    }

    fout.close();
    if(!fout) return 1;

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

    size_t buffSize = sz*nn[0]*nn[1]*nn[2];

    if(localData){
        fout.open(fname,ios_base::out|ios_base::binary|ios_base::app);
    } else {
        fname = basename+".raw";
        fout.open(fname,ios_base::out|ios_base::binary);
    }
    
    if(!fout) return 1;
    fout.write(data,buffSize);
    if(!fout) return 1;
    fout.close();
    if(!fout) return 1;    
    return 0;
}