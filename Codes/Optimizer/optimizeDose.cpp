
//
//  optimizeDose.cpp
//  Fred
//
//  Created by A. Schiavi on 09/04/13.
//  Copyright (c) 2013 University of Rome La Sapienza. All rights reserved.
//

// revamped in 2020 for Fred flash project


#include <iostream>
#include <map>
using namespace std;

#include "Optimizer.h"
#include "mhd_io.h"
#include "InputTools.h"

//=====================================================================
int main(int argc, char *argv[]){
    if(argc!=2){
        cout<<"usage: "<<argv[0]<<" inputFile"<<endl;
        cout<<"where inputFile is a file containing a list of commands and directives as follows: "<<endl;
        cout<<endl;
        cout<<'\t'<<"ptv: path plannedDose [DMF]"<<endl;
        cout<<'\t'<<"oar: path maxDose [DMF]"<<endl;
        cout<<'\t'<<"dij: path"<<endl;
        cout<<'\t'<<"NoT: maxDose DMF weight red"<<endl;
        cout<<'\t'<<"DoseCTRL: maxDose DMF weight red"<<endl;
        return 1;
    }

    system("rm -fr ./out");
    mkdir("./out");

    try{
    int ierr;

    Optimizer myOpti;
    myOpti.pthreads_max_num = 8; // number of concurrent POSIX threads (1 = serial execution)
    
    //    myOpti.EnableDebug();

    vector<string> lines = readLines(argv[1]);
    // for(auto & l : lines) cout<<l<<endl;

    float idbflg(1), initflg(0);
    float MaxIter = 100000000;
    float ConvPar = 1.e-5;
    cout<<"Parsing setup: directives:"<<endl;
    for(auto & l : lines){
        if(l.rfind("stpFlag:")==0){
            vector<string> tok=strtokens(l.substr(8));

            if((ierr=getFloatParamRequired(tok,"debug",idbflg))) throw ierr;
            if(idbflg!=0 && idbflg!=1) throw 10;

            if((ierr=getFloatParamRequired(tok,"init",initflg))) throw ierr;
            if(initflg!=0 && initflg!=1) throw 20;

            if((ierr=getFloatParam(tok,"iteration",MaxIter,100000))) throw ierr;
            if(MaxIter<=0) throw 30;

            if((ierr=getFloatParam(tok,"converg",ConvPar,1.e-5))) throw ierr;
            if(ConvPar<=0) throw 10;

            if((ierr=getIntParam(tok,"threads",myOpti.pthreads_max_num,12))) throw ierr;
            if(myOpti.pthreads_max_num<=0) throw 10;

      }
    }
    cout<<"----------------------------------------------------------------------"<<endl;
    
    myOpti.setConv(ConvPar);
    myOpti.setIterMax((int)MaxIter);

    cout<<"----------------------------------------------------------------------"<<endl;
    string outFlag = "def";
    cout<<"Parsing output: directives:"<<endl;
    for(auto & l : lines){
        if(l.rfind("outFlag:")==0){
            vector<string> tok=strtokens(l.substr(8));
            outFlag = tok[0].substr(0,tok[0].size());
        }
    }
    cout<<" Flagging the output with following flag: "<<outFlag.data()<<endl;
    myOpti.bookHistos(outFlag);
    cout<<"----------------------------------------------------------------------"<<endl;

    cout<<"Parsing NT: directives:"<<endl;
    bool foundNoT = false;
    for(auto & l : lines){
        if(l.rfind("NoT:")==0){
            foundNoT = true;
            vector<string> tok=strtokens(l.substr(4));

            float DMF,maxDose,weight, threshold,ReductionFactor;
            if((ierr=getFloatParamRequired(tok,"maxDose",maxDose))) throw ierr;
            if(maxDose<=0) throw 10;
            if((ierr=getFloatParam(tok,"DMF",DMF,1.0))) throw ierr;
            if(DMF<=0) throw 10;
            if((ierr=getFloatParam(tok,"weight",weight,100.0))) throw ierr;
            if(weight<0 or weight>100) throw 20;
            if((ierr=getFloatParam(tok,"thre",threshold,5.))) throw ierr;
            if(threshold<0) throw 20;
            if((ierr=getFloatParam(tok,"red",ReductionFactor,1.))) throw ierr;
            if(threshold<0) throw 20;
        //            myOpti.loadPTV_mhd(name,tok[0],plannedDose,weight/100.,DMF);

        myOpti.DMF_NT     = DMF;
        myOpti.maxDose_NT = maxDose;
        myOpti.weight_NT  = weight/100.;
        myOpti.thre_NT    = threshold*1.e-11;
        myOpti.Reduction_NT= ((int)ReductionFactor);
        }
    }
    if(foundNoT){
        cout<<"Configuring NT with maxDose "<<myOpti.maxDose_NT<<" DMF "<<myOpti.DMF_NT<<" and weight "<<myOpti.weight_NT<<" threshold "<<myOpti.thre_NT<<" Reduction "<<myOpti.Reduction_NT<<endl;
    } else {
        cout<<'\t'<<"Normal Tissue not used"<<endl;
    }
    cout<<"----------------------------------------------------------------------"<<endl;

    cout<<"----------------------------------------------------------------------"<<endl;
    cout<<"Parsing DoseCTRL: directives:"<<endl;
    bool foundDoseCTRL = false;

    for(auto & l : lines){
        if(l.rfind("DoseCTRL:")==0){
            foundDoseCTRL = true;
            if(foundNoT){
                cout<<"error: NoT and DoseCTRL directives are mutually exclusive"<<endl;
                exit(-1);
            }
            vector<string> tok=strtokens(l.substr(4));

            float DMF,maxDose,weight,threshold,ReductionFactor;
            if((ierr=getFloatParamRequired(tok,"maxDose",maxDose))) throw ierr;
            if(maxDose<=0) throw 10;
            if((ierr=getFloatParam(tok,"DMF",DMF,1.0))) throw ierr;
            if(DMF<=0) throw 10;
            if((ierr=getFloatParam(tok,"weight",weight,100.0))) throw ierr;
            if(weight<0 or weight>100) throw 20;
            if((ierr=getFloatParam(tok,"thre",threshold,5.))) throw ierr;
            if(threshold<0) throw 20;
            if((ierr=getFloatParam(tok,"red",ReductionFactor,1.))) throw ierr;
            if(threshold<0) throw 20;

            myOpti.DMF_DCTRL     = DMF;
            myOpti.maxDose_DCTRL = maxDose;
            myOpti.weight_DCTRL  = weight/100.;
            myOpti.thre_DCTRL    = threshold/100.; // from percentage to fraction
            myOpti.Reduction_DCTRL= ReductionFactor;
        }
    }
    if(foundDoseCTRL){
        cout<<"Configuring DoseCTRL with maxDose "<<myOpti.maxDose_DCTRL<<" DMF "<<myOpti.DMF_DCTRL<<" and weight "<<myOpti.weight_DCTRL<<" threshold "<<myOpti.thre_DCTRL<<" Reduction "<<myOpti.Reduction_DCTRL<<endl;
    } else {
        cout<<'\t'<<"Dose Control Volume not used"<<endl;
    }
    cout<<"----------------------------------------------------------------------"<<endl;
    // exit(0);

    cout<<"Parsing ptv: directives:"<<endl;
    for(auto & l : lines){
        if(l.rfind("ptv:")==0){
            vector<string> tok=strtokens(l.substr(4));
            string name = tok[0].substr(0,tok[0].size()-4);
            float DMF,plannedDose,weight;
            if((ierr=getFloatParamRequired(tok,"Dgoal",plannedDose))) throw ierr;
            if(plannedDose<=0) throw 10;
            if((ierr=getFloatParam(tok,"DMF",DMF,1.0))) throw ierr;
            if(DMF<=0) throw 10;
            if((ierr=getFloatParam(tok,"weight",weight,100.0))) throw ierr;
            if(weight<0 or weight>100) throw 20;
            myOpti.loadPTV_mhd(name,tok[0],plannedDose,weight/100.,DMF);
        }
    }
    
    if(idbflg) cout<<"OUT size red2denseIdx = "<<myOpti.red2denseIdx.size()<<"  numVxl "<<myOpti.numVxl<<" myOpti.Von.size() "<<myOpti.Von.size()<<endl;
        
    for(auto roi : myOpti.PTVs) roi.info();
    if(myOpti.PTVs.size()==0){cout<<"error: at least a PTV is needed"<<endl; exit(1);}
    cout<<"----------------------------------------------------------------------"<<endl;

    cout<<"Parsing oar: directives:"<<endl;
    for(auto & l : lines){
        if(l.rfind("oar:")==0){
            vector<string> tok=strtokens(l.substr(4));
            string name = tok[0].substr(0,tok[0].size()-4);
            float DMF,maxDose,weight;
            if((ierr=getFloatParamRequired(tok,"maxDose",maxDose))) throw ierr;
            if(maxDose<=0) throw 10;
            if((ierr=getFloatParam(tok,"DMF",DMF,1.0))) throw ierr;
            if(DMF<=0) throw 10;
            if((ierr=getFloatParam(tok,"weight",weight,100.0))) throw ierr;
            if(weight<0 or weight>100) throw 20;
            myOpti.loadOAR_mhd(name,tok[0],maxDose,weight/100.,DMF);
        }
    }

    if(idbflg) cout<<"size red2denseIdx = "<<myOpti.red2denseIdx.size()<<"  numVxl "<<myOpti.numVxl<<" myOpti.Von.size() "<<myOpti.Von.size()<<endl;
    
    for(auto roi : myOpti.OARs) roi.info();
    cout<<"----------------------------------------------------------------------"<<endl;
    // exit(0);

    if(idbflg) {
        cout<<endl<<"Entering builDgoal size red2denseIdx = "<<myOpti.red2denseIdx.size()<<endl;;
        for(int ii = 0; ii<(int)myOpti.red2denseIdx.size();ii++){
            if(myOpti.red2denseIdx[ii]<0 ){
                cout <<ii<< "red2denseIdx[i] = "<< myOpti.red2denseIdx[ii]<<endl;
            }
        }
    }

    // read common info on maps ---->
    vector<float32> map;
    int nn[3];
    float x0[3],L[3],l[3],u[3],f[3];
    string dtype;
    size_t pos;
    string dataFile;
    read_mhd_header(myOpti.PTV_paths[0],dtype,nn,x0,L,l,u,f,dataFile,pos);
    size_t N = 1UL*nn[0]*nn[1]*nn[2];
    // read common info on maps <----

    myOpti.buildDgoal();
    cout<<"numVxl    = "<<myOpti.numVxl<<endl;    
    cout<<"numVxlPTV = "<<myOpti.numVxlPTV<<endl;    
    cout<<"numVxlOAR = "<<myOpti.numVxlOAR<<endl;    
    myOpti.init_DMF();
    if(idbflg)
      cout<<"IN size red2denseIdx = "<<myOpti.red2denseIdx.size()<<"  numVxl "<<myOpti.numVxl<<" myOpti.Von.size() "<<myOpti.Von.size()<<endl;

    cout<<"----------------------------------------------------------------------"<<endl;
    cout<<"Parsing dij: directives:"<<endl;
    vector<string> Dij_files;
    for(auto & l : lines){
        if(l.rfind("dij:")==0){
            vector<string> tok=strtokens(l.substr(4));
            string name = tok[0].substr(0,tok[0].size()-4);
            Dij_files.push_back(tok[0]);
        }
    }

    if(idbflg) cout<<"After BuidDgoal size red2denseIdx = "<<myOpti.red2denseIdx.size()<<"  numVxl "<<myOpti.numVxl<<" myOpti.Von.size() "<<myOpti.Von.size()<<endl;


    if(foundDoseCTRL){
        myOpti.selectDoseCtrlVolume_bin(Dij_files);

        map.clear(); map.resize(N,0);
        ROI &roi = myOpti.OARs.back();
        for(int64 i=0;i<(int64)roi.Vi.size();i++){map[roi.Vi[i]]=roi.wei[i];}
        write_mhd("out/DoseCTRL.mhd",map,nn,x0,L,l,u,f);
        // exit(0);
    }
    
    if(foundDoseCTRL){
        myOpti.selectNT_bin(Dij_files);
        if(idbflg) {
          cout<<endl<<"After selectNT size red2denseIdx = "<<myOpti.red2denseIdx.size()<<endl;
        }
    }

    cout<<endl;
    cout<<"++++++++++        CHECK       ++++++++++"<<endl;
    cout<<"numVxl    = "<<myOpti.numVxl<<endl;    
    cout<<"numVxlPTV = "<<myOpti.numVxlPTV<<endl;    
    cout<<"numVxlOAR = "<<myOpti.numVxlOAR<<endl;    
    cout<<"++++++++++++++++++++++++++++++++++++++++"<<endl;
    cout<<endl;
    // exit(0);


    cout<<endl<<"Entering load Dij bin "<<endl;
    myOpti.loadDij_bin(Dij_files);
    if(idbflg) {
      cout<<"After loadDij_bin size red2denseIdx = "<<myOpti.red2denseIdx.size()<<"  numVxl "<<myOpti.numVxl<<" myOpti.Von.size() "<<myOpti.Von.size()<<endl;
      for(int ii = 0; ii<(int)myOpti.red2denseIdx.size();ii++){
    if(myOpti.red2denseIdx[ii]<0 ){
      cout<<"WARNING : for voxel" <<ii<< "  red2denseIdx[i] = "<< myOpti.red2denseIdx[ii]<<endl;
    }
      }
    }
    cout<<"----------------------------------------------------------------------"<<endl;
    cout<<"Parsing fluences: directive"<<endl;
    bool lfluences = false;
    string initFluences = "ones";
    for(auto & l : lines){
        if(l.rfind("fluences:")==0){
            if(lfluences) throw 50;
            lfluences = true;
            vector<string> tok=strtokens(l.substr(9));
	    initFluences = tok[0].substr(0,tok[0].size()-9);

	    if(initFluences=="flat") {
	      myOpti.pthreads_max_num = 1; // number of concurrent POSIX threads (1 = serial execution)
	      myOpti.f_flatOptimizer = true;
	    }

	}
    }
    myOpti.initParticleNumbers(initFluences);
    cout<<"----------------------------------------------------------------------"<<endl;
    if(initflg>0){
      (myOpti.f_out)->Write();
      (myOpti.f_out)->Close();
      
      exit(0);
    }
    cout <<"  Fluences initializated "<<endl;

    map.clear(); map.resize(N,0);
    for(int64 i=0;i<myOpti.numVxl;i++){map[myOpti.red2denseIdx[i]]=myOpti.Von[i];}

    string fname = "rois_"+outFlag+".mhd";
    write_mhd(fname,map,nn,x0,L,l,u,f);

    map.clear(); map.resize(N,0);
    for(int64 i=0;i<myOpti.numVxl;i++){map[myOpti.red2denseIdx[i]]=myOpti.DMF[i];}

    fname = "alldmf_"+outFlag+".mhd";
    write_mhd(fname,map,nn,x0,L,l,u,f);

    // exit(0);
    
    myOpti.run();

    myOpti.buildOptiDose(outFlag,1);

    myOpti.writePBs(outFlag);

    // myOpti.printPBs();
    }
    catch (int ie){
        cerr<<"error: ";
        switch(ie){
            case 1: // missing required parameter
            cerr<<"missing required parameter "<<lastParsedParameter<<endl;
            break;
            case -1: // missing value
            cerr<<"missing value for parameter "<<lastParsedParameter<<endl;
            break;
            case -2: // wrong syntax
            cerr<<"wrong syntax for parameter "<<lastParsedParameter<<endl;
            cerr<<"offending input: "<<lastParsedDefinition<<endl;
            break;
            case 10: // negative value not allowed
            cerr<<lastParsedParameter<<" must be a greater than zero"<<endl;
            cerr<<"offending input: "<<lastParsedDefinition<<endl;
            break;
            case 20: // wei out of range
            cerr<<lastParsedParameter<<" must be in the [0,100] range"<<endl;
            cerr<<"offending input: "<<lastParsedDefinition<<endl;
            break;
            case 50: // too many fluences: 
            cerr<<"only one fluences: directive is allowed "<<endl;
            break;
            default:
            cerr<<"ops...unknown error..."<<endl;
        }
        exit(1);
    }
}
