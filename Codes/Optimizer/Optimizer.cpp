//
//  Optimizer.cpp
//  Fred
//
//  Created by A. Schiavi on 09/04/13.
//  Copyright (c) 2013 University of Rome La Sapienza. All rights reserved.
//
//


#include "Optimizer.h"
#include "InputTools.h"
#include "mhd_io.h"
#include "TROOT.h"
#include <set>

#define MIN( A, B) (((A)<(B))?(A):(B))
#define MAX( A, B) (((A)>(B))?(A):(B))

//////////////////////////////////////////////////////////////////////
////////////////////  PTHREAD STRUCTURES and PROTOTYPES  /////////////
#include <pthread.h>
typedef struct threadInfo_s {
    int thread_num;
    int num_threads;
    Optimizer *optObj;
} threadInfo;
pthread_mutex_t buildD_mutex;

void *takeStep_thread(void *info) ;
void *buildD_thread(void *info);
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
float getLowerCutoffFromCumulativeDistribution(vector<float> &di,float thres)
{
    const bool debug =  false;
    // get min and max and sum
    float min,max;
    double sum=0;
    min=max=di[0];
    for(auto v : di){
        if(min>v) min=v;
        if(max<v) max=v;
        sum += v;
    }
    // log10 spaced binning
    float logmin = log10(min);
    float logmax = log10(max);
    int nsub = 50; // subdivisions per decade
    float hbin = 1.0/nsub; // bin spacing
    int nbin = (logmax-logmin)/hbin+1;
    vector<double> hist(nbin);
    // fill histogram
    for(auto v : di){
        int ibin = log10(v/min)*nsub;
        hist[ibin]+=v;
    }
    // normalize
    for(int i=0;i<hist.size();i++) hist[i]/=sum;
    // get cumulative
    for(int i=1;i<hist.size();i++) hist[i]+=hist[i-1];
    // find lower cutoff for given threshold
    int icut;
    for(int i=0;i<hist.size();i++) {
        icut=i;
        if(hist[i]>=thres) break;
    }
    float logcut = logmin+hbin*icut;
    float cutoff = pow(10.0f,logcut);
    double check = 0;
    for(auto v : di){
        if(v<=cutoff) check+=v;
    }
    if(debug) cout<<"cutoff cumul: "<<'\t'<<min<<' '<<max<<' '<<sum<<endl;
    if(debug) cout<<'\t'<<logmin<<' '<<logmax<<' '<<nbin<<' '<<icut<<' '<<logcut<<' '<<cutoff<<endl;
    if(debug) cout<<'\t'<<check<<' '<<check/sum<<endl<<endl;
    return cutoff;
}


void ROI::info(){
    cout<<"ROI: "<<ID<<endl;
    if(plannedDose>0) {
        cout<<"\ttype = PTV"<<endl;
        cout<<"\tplanned dose = "<<plannedDose<<endl;
    } else {
        cout<<"\ttype = OAR"<<endl;
    }
    cout<<"\tmax dose = "<<maxDose<<endl;
    cout<<"\tnum active voxels = "<<Vi.size()<<endl;

    double sum=0; 
    for(auto w : wei) sum+=w;
    cout<<"\tweight sum = "<<sum<<endl;

    cout<<"\tDMF = "<<DMF<<endl;

    cout<<endl;
}

void ROI::weightRescale(float weight)
{
    float ww = weight/wei.size();  // asch Apr 2023 volume-normalized weight
    for(size_t i=0;i<wei.size();i++) wei[i]*=ww; 
}


//////////////////////////////////////////////////////////////////////////////////////////

int Optimizer::setReductionNT(float64 f){if((f<=0)&&(f>1)) return -1; Reduction_NT= f; return 0;}
int Optimizer::setConv(float64 f){if(f<=0) return -1; conv= f; return 0;}
int Optimizer::setScaleFactor(float64  f){if(f<=0) return -1; fScaleOptimizationStep=f; return 0;}
int Optimizer::setIterMax(int Imax){if(Imax<=0) return -1; itermax=Imax; return 0;}


void Optimizer::loadROI_txt(string fpath){
    ifstream fin(fpath);
    if(!fin){cerr<<"Error: could not open "<<fpath<<endl; exit(-1);}
    if(fin.peek()=='#') fin.ignore(4096,'\n');
    
    ROI roi;
    fin>>roi.ID>>roi.plannedDose>>roi.maxDose;
    roi.Vi.reserve(100000);
    roi.wei.reserve(100000);
    int64 ivxl;
    float32 wei;
    while(1){
        if(fin.eof()) break;
        fin>>ivxl>>wei;
        if(!fin) break; // error
        // cout<<ivxl<<' '<<wei<<endl;
        roi.Vi.push_back(ivxl);
        roi.wei.push_back(wei);
    }

    if(roi.plannedDose>0) PTVs.push_back(roi);
    else OARs.push_back(roi);

}



void Optimizer::loadROI_mhd(ROI &roi,string fpath){
    vector<float32> fmap;
    int nn[3];
    float x0[3],L[3],l[3],u[3],f[3];
    string dtype;
    size_t pos;
    string dataFile;

    read_mhd_header(fpath,dtype,nn,x0,L,l,u,f,dataFile,pos);
    fmap.resize(1UL*nn[0]*nn[1]*nn[2]);
    if(dtype=="MET_SHORT"){
        vector<int16> map;
        read_mhd(fpath,map,nn,x0,L,l,u,f);
        for(size_t i=0;i<fmap.size();i++) fmap[i]=map[i];
    } else
    if(dtype=="MET_INT"){
        vector<int32> map;
        read_mhd(fpath,map,nn,x0,L,l,u,f);
        for(size_t i=0;i<fmap.size();i++) fmap[i]=map[i];
    } else
    if(dtype=="MET_FLOAT"){
        read_mhd(fpath,fmap,nn,x0,L,l,u,f);
    } else {
        cerr<<"Error: data type "<<dtype<<" for "<<fpath<<" not yet supported"<<endl;
        exit(1);
    } 

    if(numVxlAll<0) { // not set => first ROI
        numVxlAll = 1UL*nn[0]*nn[1]*nn[2];
        for(int i=0;i<3;i++) {gnn[i]=nn[i];}
    } else {   // check that geometry is aligned with other ROIs
        if(numVxlAll != 1UL*nn[0]*nn[1]*nn[2]) {
            cerr<<"Error: number of voxels are different"<<endl;
            cerr<<"Expected: "<<numVxlAll<<endl;
            cerr<<"Found for "<<fpath<<" : "<<1UL*nn[0]*nn[1]*nn[2]<<endl;
        }
    }

    // store only nonzero voxels for each roi
    size_t nonzero=0;
    for(size_t i=0;i<fmap.size();i++) if(fmap[i]!=0) nonzero++; // count voxels != 0
    roi.Vi.resize(nonzero); 
    roi.wei.resize(nonzero);
    nonzero=0;
    for(size_t i=0;i<fmap.size();i++) {
        if(fmap[i]!=0) {
            roi.Vi[nonzero]=i;
            roi.wei[nonzero]=fmap[i]; // map value is the local weight
            nonzero++;
        }
    }

    // for(size_t i=0;i<ROI.size();i++) ROInew[i] = ROI[i]>0? 0 : 1;
    // write_mhd("aaa.mhd",ROInew,nn,x0,L,l,u,f);    
}

#define VOXEL_OFF 0
#define VOXEL_PTV 1
#define VOXEL_OAR 2


void Optimizer::loadPTV_mhd(string ID,string fpath,float plannedDose,float weight, float DMF){
    ROI ptv;
    ptv.ID = ID;
    ptv.plannedDose = plannedDose;
    ptv.maxDose = plannedDose;
    ptv.DMF = DMF;
    loadROI_mhd(ptv,fpath);
    ptv.weightRescale(weight);
    PTVs.push_back(ptv);
    PTV_paths.push_back(fpath);
}

void Optimizer::loadOAR_mhd(string ID,string fpath,float maxDose,float weight, float DMF){
    ROI oar;
    oar.ID = ID;
    oar.plannedDose = 0;
    oar.maxDose = maxDose;
    oar.DMF = DMF;
    loadROI_mhd(oar,fpath);
    oar.weightRescale(weight);
    OARs.push_back(oar);
    OAR_paths.push_back(fpath);
}


void Optimizer::loadDij_txt(string fpath){

  ifstream fin(fpath);
    if(!fin){cerr<<"Error: could not open "<<fpath<<endl; exit(-1);}
    if(fin.peek()=='#') fin.ignore(4096,'\n');
    fin>>npb>>numVxl;
    Dij.clear(); Dij.resize(1UL*numVxl*npb,0); // allocate dense Dij matrix and set it to zero
    pbID.clear(); pbID.resize(npb,0);
    fieldID.clear(); fieldID.resize(npb,0);
    isPBActive.clear(); isPBActive.resize(npb,true);

    for (int j=0;j<npb;j++){
        fin>>pbID[j]>>fieldID[j];
        int nvxl;
        fin>>nvxl;
        for(int k=0;k<nvxl;k++){
            int64 ivxl; float64 dose;
            fin>>ivxl>>dose;
            Dij[ivxl+j*numVxl]=dose;
        }
    }

    if(f_Debug) {
      for(int i=0;i<numVxl;i++){
        for(int j=0;j<npb;j++){
	  cout<<Dij[i+j*numVxl]<<' ';
        }
        cout<<endl;
      }
    }
    D.resize(numVxl);
    initParticleNumbers("ones");
}



void Optimizer::selectNT_bin(vector<string> &fpaths){
  

    Dij_paths.clear();
    Dij_npb.clear();

    vector<FILE *> files;
    vector<int> Dij_npb;
    FILE *fin;

    for(auto & path : fpaths){
        cout<<"\t"<<"file = "<<path<<endl;
	fin = fopen(path.c_str(),"rb");
	if(!fin){
	  cerr<<"Error: could not open file: "<<path<<endl;
	  exit(1);
	}
        Dij_paths.push_back(path);
        files.push_back(fin);
	
	int nn[3];
	fread(nn,3,sizeof(int),fin);
	cout<<"nn = "<<nn[0]<<' '<<nn[1]<<' '<<nn[2]<<endl;
	
	for(int i=0;i<3;i++) {
	  if(gnn[i]!=nn[i]){
            cerr<<"Error: Dij grid dimensions do not match with ROI dims"<<endl;
            exit(1);
	  }
	}
	
	
	float hs[3];
	fread(hs,3,sizeof(float),fin);
	cout<<"hs = "<<hs[0]<<' '<<hs[1]<<' '<<hs[2]<<endl;
	
	float x0[3];
	fread(x0,3,sizeof(float),fin);
	cout<<"x0 = "<<x0[0]<<' '<<x0[1]<<' '<<x0[2]<<endl;
	
	fread(&npb,1,sizeof(int),fin);
	cout<<"npb = "<<npb<<endl;
        Dij_npb.push_back(npb);
        cout<<endl;
    }
	
    vector<int> Vi(numVxlAll);
    vector<float> di(numVxlAll);

    int ipb,ifield,ni;
    int numVxlOverThrNT=0, numVxlSelectNT=0, num_crazy_VXL=0;
    
    cout<<"AllVxl:: "<<numVxlAll<<" NumVxl before NoT "<<numVxl<<endl<<endl;
    for (int i=0;i<(int)files.size();i++){
        fin = files[i];

	for(int j=0;j<Dij_npb[i];j++){
	  fread(&ipb,1,sizeof(int),fin);
	  fread(&ifield,1,sizeof(int),fin);
	  fread(&ni,1,sizeof(int),fin);
	  Vi.resize(ni);
	  di.resize(ni);
	
	  fread(Vi.data(),ni,sizeof(int)  ,fin);
	  fread(di.data(),ni,sizeof(float),fin);

	  for(int i=0;i<ni;i++) {
	    //	cout<<" Vx for PB:: "<<i<<" "<<Vi[i]<<endl;
	    int ired = dense2redIdx[Vi[i]];

	    //check if ired && Von[ired]!=VOXEL_PTV && Von[ired]!=VOXEL_OAR
	    //if(ired != -1  && Von[ired]!=VOXEL_PTV && Von[ired]!=VOXEL_OAR && di[i]>thre_NT){
	      
	      
	    if( (ired == -1) && (di[i]>thre_NT) ){		
	      numVxlOverThrNT++;
	      if(numVxlOverThrNT%Reduction_NT==0){
		numVxlSelectNT++;
		dense2redIdx.at(Vi[i]) = numVxl;
		red2denseIdx.push_back(Vi[i]);		
		numVxl++; 		
		Dgoal.push_back(maxDose_NT); 
		weight.push_back(weight_NT);
		Von.push_back(VOXEL_OAR);
	      }
	    }
	    /*
	      the code should not enter in the next loop!!!
	    */	    
	    else if(ired != -1  && Von[ired]!=VOXEL_PTV && Von[ired]!=VOXEL_OAR){
	      cout<<"SelectNT: crazy voxel-> ired= "<<ired<<" Vi[i]= "<<Vi[i]<<" i= "<<i<<
		" Von= "<<Von[ired]<< endl;
	      num_crazy_VXL++;
	      if(di[i]>thre_NT){// crazy!!
		
		dense2redIdx.at(Vi[i]) = numVxl;
		red2denseIdx.push_back(Vi[i]);
	      
		numVxl++; 
		
		Dgoal.push_back(maxDose_NT); 
		weight.push_back(weight_NT);
		Von.push_back(VOXEL_OAR);
	      }
	    }
	  }
	}
    }
    cout<<"NT voxel Over Thr= "<<numVxlOverThrNT<<" NT voxel selected= "<<numVxlSelectNT<<
      " num crazy voxel= "<<num_crazy_VXL<<endl;
    cout<<"NumVxl after NoT= "<<numVxl<<"  numVxl OAR= "<<numVxl<<endl;
    cout<<"Check exiting selectNT_bin : red2denseIdx.size()= "<<red2denseIdx.size()<< endl;
}

void Optimizer::selectDoseCtrlVolume_bin(vector<string> &fpaths){
    Dij_paths.clear();
    Dij_npb.clear();

    vector<FILE *> files;
    vector<int> Dij_npb;
    FILE *fin;

    for(auto & path : fpaths){
        cout<<"\t"<<"file = "<<path<<endl;
        fin = fopen(path.c_str(),"rb");
        if(!fin){
          cerr<<"Error: could not open file: "<<path<<endl;
          exit(1);
        }
        Dij_paths.push_back(path);
        files.push_back(fin);
    
        int nn[3];
        fread(nn,3,sizeof(int),fin);
        cout<<"nn = "<<nn[0]<<' '<<nn[1]<<' '<<nn[2]<<endl;
    
        for(int i=0;i<3;i++) {
          if(gnn[i]!=nn[i]){
                cerr<<"Error: Dij grid dimensions do not match with ROI dims"<<endl;
                exit(1);
          }
        }
    
    
        float hs[3];
        fread(hs,3,sizeof(float),fin);
        cout<<"hs = "<<hs[0]<<' '<<hs[1]<<' '<<hs[2]<<endl;
        
        float x0[3];
        fread(x0,3,sizeof(float),fin);
        cout<<"x0 = "<<x0[0]<<' '<<x0[1]<<' '<<x0[2]<<endl;
        
        fread(&npb,1,sizeof(int),fin);
        cout<<"npb = "<<npb<<endl;
        Dij_npb.push_back(npb);
        cout<<endl;
    }
    
    vector<int> Vi;
    vector<float> di;
    vector<int> DOSECTRL(numVxlAll,0);

    int ipb,ifield,ni;
    
    cout<<"AllVxl:: "<<numVxlAll<<" NumVxl before DoseCTRL "<<numVxl<<endl<<endl;

    cout<<"Parsing Dij for Dose Control Volume: ";
    //    double t0=1.0*clock()/CLOCKS_PER_SEC;
    for (int i=0;i<(int)files.size();i++){
        fin = files[i];

        for(int j=0;j<Dij_npb[i];j++){
            fread(&ipb,1,sizeof(int),fin);
            fread(&ifield,1,sizeof(int),fin);
            fread(&ni,1,sizeof(int),fin);
            Vi.resize(ni);
            di.resize(ni);

            fread(Vi.data(),ni,sizeof(int)  ,fin);
            fread(di.data(),ni,sizeof(float),fin);

            float doseCutoff = getLowerCutoffFromCumulativeDistribution(di,thre_DCTRL);
            for(int i=0;i<Vi.size();i++){
                if(di[i]>doseCutoff) DOSECTRL[Vi[i]]=1; // flag voxel
            }
            cout<<'.';
            cout.flush();
        }
    }
    cout<<endl;
    // double t1=1.0*clock()/CLOCKS_PER_SEC;
    // cout<<"timing "<<t1-t0<<endl;
    // exit(0);
    int numVxlDoseCTRL=0;
    for(int i=0;i<DOSECTRL.size();i++) {
        if(DOSECTRL[i]>0) numVxlDoseCTRL++;
    }

    cout<<"num voxel in Dose Control Volume before reduction "<<numVxlDoseCTRL<<endl;
    cout<<"purging voxels already under control (PTV or OAR): ";
    for(auto i : red2denseIdx){
        if(DOSECTRL[i]>0) numVxlDoseCTRL--;
        DOSECTRL[i]=0;
    }
    cout<<numVxlDoseCTRL<<endl;
    cout<<"num voxel in PTVs: "<<numVxlPTV<<endl;
    cout<<"required voxel percentage wrt PTVs: "<<Reduction_DCTRL*100<<"%"<<endl;
    int numVxlToSelect = numVxlPTV*Reduction_DCTRL;
    numVxlToSelect = MAX(numVxlToSelect,0);
    numVxlToSelect = MIN(numVxlToSelect,numVxlDoseCTRL);
    // numVxlToSelect = numVxlDoseCTRL;
    cout<<"num voxel in Dose Control Volume: "<<numVxlToSelect<<endl;

    ROI oar;
    oar.ID = "DoseCTRL";
    oar.plannedDose = 0;
    oar.maxDose = maxDose_DCTRL;
    oar.DMF = DMF_DCTRL;
    // loadROI_mhd(oar,fpath);

    vector<int> vv;
    vv.reserve(numVxlDoseCTRL);
    for(int i=0;i<DOSECTRL.size();i++) {
        if(DOSECTRL[i]>0) vv.push_back(i);
    }

    oar.Vi.resize(numVxlToSelect,-1);
    // for(int i=0;i<numVxlToSelect;i++) oar.Vi[i]=vv[i];
    oar.wei.resize(numVxlToSelect,1);

    cout<<"random select voxels in original Dose Control Volume"<<endl;
    // random select voxels in original Dose Control Volume
    for(int i=0;i<numVxlToSelect;i++){
        float r = 1.0*rand()/RAND_MAX;
        int idx = int(r*vv.size());
        oar.Vi[i] = vv[idx];
        vv.erase(vv.begin()+idx);
    }

    oar.weightRescale(weight_DCTRL);
    OARs.push_back(oar);
    OAR_paths.push_back("out/DoseCTRL.mhd");

    oar.info();

    // add OAR to reduced matrices 
    for(int i=0;i<oar.Vi.size();i++){
        dense2redIdx[oar.Vi[i]] = numVxl++;
        numVxlOAR++;
        red2denseIdx.push_back(oar.Vi[i]);
        Dgoal.push_back(oar.maxDose); 
        weight.push_back(oar.wei[i]);
        Von.push_back(VOXEL_OAR);
    }
    cout<<"NumVxl after DoseCTRL= "<<numVxl<<"  numVxl OAR= "<<numVxlOAR<<endl;
}

void Optimizer::loadDij_bin(vector<string> &fpaths){

    Dij_paths.clear();
    Dij_npb.clear();

    vector<FILE *> files;
    vector<int> Dij_npb;
    FILE *fin;

    for(auto & path : fpaths){
        cout<<"\t"<<"file = "<<path<<endl;
        fin = fopen(path.c_str(),"rb");
        if(!fin){
            cerr<<"Error: could not open file: "<<path<<endl;
            exit(1);
        }
        Dij_paths.push_back(path);
        files.push_back(fin);

        int nn[3];
        fread(nn,3,sizeof(int),fin);
        

        for(int i=0;i<3;i++) {
            if(gnn[i]!=nn[i]){
                cerr<<"Error: Dij grid dimensions do not match with ROI dims"<<endl;
                cout<<'\t'<<"nn = "<<nn[0]<<' '<<nn[1]<<' '<<nn[2]<<endl;
                cout<<'\t'<<"nn = "<<gnn[0]<<' '<<gnn[1]<<' '<<gnn[2]<<endl;
                exit(1);
            }
        }
        float hs[3];
        fread(hs,3,sizeof(float),fin);
        // cout<<'\t'<<"hs = "<<hs[0]<<' '<<hs[1]<<' '<<hs[2]<<endl;

        float x0[3];
        fread(x0,3,sizeof(float),fin);
        // cout<<'\t'<<"x0 = "<<x0[0]<<' '<<x0[1]<<' '<<x0[2]<<endl;

        fread(&npb,1,sizeof(int),fin);
        cout<<'\t'<<"npb = "<<npb<<endl;
        Dij_npb.push_back(npb);
        cout<<endl;
    }

    vector<int> Vi(numVxlAll);
    vector<float> di(numVxlAll);

    int ipb,ifield,ni;
    
    npb=0; // reset
    for(auto & n : Dij_npb) npb+=n;
    cout<<"Total number of PB: "<<npb<<endl;

    // exit(0);

    // matrices are in column-major order as in Fred
    Dij.clear(); Dij.resize(1UL*numVxl*npb,0); // allocate reduced Dij matrix and set it to zero

    cout<<"We managed to allocate # PB: "<<npb<<endl;
    pbID.clear(); pbID.resize(npb,0);
    fieldID.clear(); fieldID.resize(npb,0);
    isPBActive.clear(); isPBActive.resize(npb,true); 


    cout<<"We are starting the loop "<<endl;
    int jcurr=0;
    for (int i=0;i<(int)files.size();i++){
        fin = files[i];
	//	cout<<" File:: "<<i<< " "<<npb<<" "<<Dij_npb[i]<<endl;
        for(int j=0;j<Dij_npb[i];j++){
            fread(&ipb,1,sizeof(int),fin);
            fread(&ifield,1,sizeof(int),fin);
            fread(&ni,1,sizeof(int),fin);
            Vi.resize(ni);
            di.resize(ni);
            pbID[jcurr]=ipb;
            fieldID[jcurr]=ifield;
            fread(Vi.data(),ni,sizeof(int)  ,fin);
            fread(di.data(),ni,sizeof(float),fin);
            for(int i=0;i<ni;i++) {
                int ired = dense2redIdx[Vi[i]];
                if(ired>=0){
		  //		  cout<<"  "<<ired<<" "<<ired+jcurr*numVxl<<" "<<npb*numVxl<<endl;
		  Dij[ired+jcurr*numVxl]=di[i]*DMF[Vi[i]]; // set voxel dose applying DMF
                }
            }      
            jcurr++; // update current global spot index
        }

    }    
    
    D.resize(numVxl);
    //    initParticleNumbers("ones");

    // exit(0);
  
    cout<<"Dij read from file(s): ";
    for(auto & path : Dij_paths) cout<<path<<' ';
    cout<<endl;
    printf("\tnum of reduced matrix elements: %lld\n",numVxl*npb);
    size_t Bytes = 1UL*numVxl*npb*sizeof(Dij[0]);
    size_t KB = Bytes/1024;
    size_t MB = KB/1024;
    size_t GB = MB/1024;
  
    if(GB>0)
        printf("\tmemory allocated: %.1f GB\n",MB/1024.);
    else if(MB>0)
        printf("\tmemory allocated: %.1f MB\n",KB/1024.);
    else
        printf("\tmemory allocated: %lu KB\n",KB);

    for (int i=0;i<(int)files.size();i++){
        fclose(files[i]);
    }

    // exit(0);
}

void  Optimizer::init_DMF()
{
    DMF.clear(); DMF.resize(numVxlAll,DMF_NT); // allocate dense map for DMF setting it to normal tissue DMF

    // update DMF for OAR(s)
    for(auto & roi : OARs) {
      for(int k=0;k<(int)roi.Vi.size();k++) {
            DMF[roi.Vi[k]]=roi.DMF;
        }
    }

    // PTV should be the last to be updated.    
    // update DMF for PTV(s)
    for(auto & roi : PTVs) {
      for(int k=0;k<(int)roi.Vi.size();k++) {
            DMF[roi.Vi[k]]=roi.DMF;
        }
    }

}



void Optimizer::buildDgoal(){

  numVxlPTV=numVxlOAR=0; 
  for(auto roi : PTVs) numVxlPTV+=roi.Vi.size();
  for(auto roi : OARs) numVxlOAR+=roi.Vi.size();

  numVxl = numVxlPTV+numVxlOAR;
  
  Dgoal.clear(); Dgoal.resize(numVxl,0);
  Von.clear(); Von.resize(numVxl,0); // voxel label 0=off, 1=on (PTV), 2=on (OAR)
  weight.clear(); weight.resize(numVxl,1);
  dense2redIdx.clear(); dense2redIdx.resize(numVxlAll,-1);
  red2denseIdx.clear(); red2denseIdx.resize(numVxl,-1);

  int32 i = 0;
  int32 iDou = 0;
  for(auto roi : PTVs) {
    for(int k=0;k<(int)roi.Vi.size();k++,i++) {
      dense2redIdx[roi.Vi[k]] = i;
      red2denseIdx[i] = roi.Vi[k];
      Dgoal[i]=roi.plannedDose; 
      weight[i]=roi.wei[k];
      Von[i]=VOXEL_PTV; 
    }
  }

  for(auto roi : OARs) {
    for(int k=0;k<roi.Vi.size();k++) {
      if (dense2redIdx[roi.Vi[k]] != -1)   { // already assigned to PTV
	//PTV has the priority.
	iDou++;
	numVxl--;
      } else {
	//before assignining...
	dense2redIdx[roi.Vi[k]] = i;
	red2denseIdx[i] = roi.Vi[k];
	Dgoal[i]=roi.maxDose;
	weight[i]=roi.wei[k];
	Von[i]=VOXEL_OAR;
	i++;
      }
    }
  }
  red2denseIdx.resize(numVxl);
  Von.resize(numVxl);
  Dgoal.resize(numVxl);
  weight.resize(numVxl);
  cout<<iDou<<" voxels are flagged both as PTV & OAR!!!  "<<endl;
  
  cout<<" This is interesting! Who enters DGoal? "<<i<<endl;
    
}
    
void Optimizer::initParticleNumbers(string startFile){
    N.clear(); N.resize(npb,1);

    if(startFile=="ones"){
        cout<<"\tsetting initial fluences to 1 [default]"<<endl;
        // do nothing => default
    }
    else if(startFile=="flat"){
        cout<<"\tsetting initial fluences to a random value flat for all pb of a given field"<<endl;
        srand((unsigned) time(NULL)); // this is used by the STL
        for(int ifi=0;ifi<(int)7;ifi++) {
	  double fi_flue = 1.+9.*rand()/RAND_MAX;
	  for(int i=0;i<(int)N.size();i++) {
	    if(fieldID[i] != ifi) continue;
	    N[i]= fi_flue;
	  }
	}
    }
    else if(startFile=="random"){
        cout<<"\tsetting initial fluences to a random value between 1. and 10."<<endl;
        srand((unsigned) time(NULL)); // this is used by the STL
        for(int i=0;i<(int)N.size();i++) N[i]= 1.+9.*rand()/RAND_MAX;
    }
    else { // load initial guess from file
        cout<<"\tsetting initial fluences using file "<<startFile<<endl;
        vector<string> lines = readLines(startFile);
        if (lines.empty()) exit(1);

        for(int i=0;i<(int)N.size();i++) N[i]= -1; // invalidate all
        int idx,fid,pbid;
        float fluence;


        for(auto & line : lines){
            istringstream iss(line);
            iss>>idx>>fid>>pbid>>fluence;
            if(!iss){
                cerr<<"error: cannot parse input line: "<<line<<endl;
                exit(1);
            }
            // cout<<idx<<' '<<fid<<' '<<pbid<<' '<<fluence<<endl;
            for(int j=0;j<(int)N.size();j++){
                if(fid==fieldID[j] and pbid==pbID[j]){
                    if(N[j]<0){
                        N[j]=fluence; 
                    } else {
                        cerr<<"error:  multiple definitions for FID="<<fieldID[j]<<" PBID="<<pbID[j]<<endl;
                        exit(1);
                    }
                    
                    break;
                }
            }
        }
        // final check that all PB have been set
        int notSet=0;
        for(int j=0;j<(int)N.size();j++){
            if(N[j]<0) {
                notSet++;
                cerr<<"error: FID="<<fieldID[j]<<" PBID="<<pbID[j]<<" not set"<<endl;
            }
        }
        if(notSet>0){
            cerr<<"error: not all PBs have been set => check input file "<<startFile<<endl;
            exit(1);
        }

    }
    Nnew=N;
    // exit(0);
}

void Optimizer::report(){

    double s0=0,s1=0,delta;
    double deltamin = 1e99, deltamax = -1e99;
    buildD();
    for (int i=0; i<numVxl; i++) {
        if(Von[i]!=VOXEL_PTV) continue; // skip this voxel

        if(Dgoal[i]>0.1){
            s0++;
            delta = fabs((D[i]-Dgoal[i])/Dgoal[i]);

	    ((TH1D*)gDirectory->Get("DoseDiscrepant"))->Fill(Dgoal[i]);
	    ((TH1D*)gDirectory->Get("DeltaDiscrepant"))->Fill(delta);
	    
            s1 += delta;
            deltamax = std::max(deltamax,D[i]/Dgoal[i]);
            deltamin = std::min(deltamin,D[i]/Dgoal[i]);
        }
    }
    if(s0<2) return; // no statistics!

    double deltaavg = s1/s0;
    cout<<"  PTV: "<<setw(3)<<" Deviation Avg: "<<setw(6)<<deltaavg*100.
    <<"% Min: "<<setw(6)<<(deltamin-1.)*100.
    <<"% Max: "<<setw(6)<<(deltamax-1.)*100.<<"%"<<endl;    

    f_out->Write();
    f_out->Close();
}


void Optimizer::buildD(){
    // old serial code >>>
    // for(int i=0;i<numVxl;i++){
    //     D[i]=0;
    //     for(int j=0;j<npb;j++){
    //         D[i]+=Dji[i*npb+j]*Nnew[j];
    //    }
    // }
    // return;
    // old serial code <<<

    for(int i=0;i<numVxl;i++) D[i]=0; // reset dose 

    vector<pthread_t> threads(pthreads_max_num);
    vector<threadInfo> tInfo(pthreads_max_num);
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    pthread_mutex_init(&buildD_mutex, NULL);

    for (int i=0; i<pthreads_max_num; i++) {
        tInfo[i].thread_num = i;
        tInfo[i].num_threads = pthreads_max_num;
        tInfo[i].optObj = this;
        pthread_create(&threads[i], &attr, buildD_thread, (void *) &tInfo[i]);
    }

    /* Wait for all threads to complete */ 
    for (int i=0; i<pthreads_max_num; i++) {
        pthread_join(threads[i], NULL);
    }
    
}

void Optimizer::MonitorDoseRatios(){
  double dose_in(0.), dose_out(0.);

  for(int j=0;j<npb;j++){
    //        isPBActive[j]=false;
    dose_in = dose_out = 0;
    for(int i=0;i<numVxl;i++){
      if(Dij[i+j*numVxl]>0){
	if(Von[i]==VOXEL_PTV){
	  //                isPBActive[j]=true;
	  //                break;
	  dose_in += Dij[i+j*numVxl]*Nnew[j];
	} else {
	  dose_out += Dij[i+j*numVxl]*Nnew[j];
	}
      }
    }
    ((TH1D*)gDirectory->Get("DoseRatioPBS"))->Fill(dose_out/(dose_in+dose_out));
    
  }

}

void Optimizer::selectActivePB(){

    for(int j=0;j<npb;j++){
        isPBActive[j]=false;
        for(int i=0;i<numVxl;i++){
            if(Von[i]==VOXEL_PTV and Dij[i+j*numVxl]>0){
                isPBActive[j]=true;
                break;
            }
        }
        if(not isPBActive[j]) N[j]=0; // suppress non active PBs by setting fluence to zero
    }
    
}


void Optimizer::printD(){
    cout<<"Dose:"<<endl;
    for(int i=0;i<numVxl;i++){
        cout<<"\t"<<D[i]<<endl;
    }
    cout<<endl;
}

void Optimizer::printDgoal(){
    cout<<"Planned dose :"<<endl;
    for(int i=0;i<numVxl;i++){
        cout<<"\t"<<Dgoal[i];
        if(Von[i]==VOXEL_PTV) cout<<" <- PTV";
        if(Von[i]==VOXEL_OAR) cout<<" <- OAR";
        cout<<endl;
    }
    cout<<endl;
}

void Optimizer::printPBs(){
    cout<<"Pencil beams :"<<endl;
    for(int j=0;j<npb;j++){
        cout<<"\t"<<j<<" ("<<fieldID[j]<<','<<pbID[j]<<") : "<<N[j];
        if(N[j]==0) {cout<<" <- suppressed";}
        else if(isPBActive[j]) cout<<" <- active";
        cout<<endl;
    }
    cout<<endl;
}

void Optimizer::writePBs(string outFlag){
    string fname = "optiPlan_" + outFlag + ".txt";
    ofstream fout(fname.c_str());
    fout<<"# isource fieldID pbID NumParticles"<<endl;
    for(int j=0;j<npb;j++){
        fout<<j<<'\t'<<fieldID[j]<<'\t'<<pbID[j]<<"\t"<<N[j];
        fout<<endl;
    }
    cout<<"Optimized plan written to "<<fname<<endl;
}

void Optimizer::logFluences(int iter){
    char fname[2048];
    snprintf(fname,sizeof(fname),"./out/flu%05d.txt",iter);
    ofstream fout(fname);
    fout<<"# isource fieldID pbID NumParticles"<<endl;
    for(int j=0;j<npb;j++){
        fout<<j<<'\t'<<fieldID[j]<<'\t'<<pbID[j]<<"\t"<<Nnew[j];
        fout<<endl;
    }
}

 float64 Optimizer::getCost(){
    float64 chi2=0,den,num;
    for (int32 i=0; i<numVxl; i++) {
        if(Von[i]==VOXEL_OFF) continue; // skip this voxel
        den = std::max(Dgoal[i], 1e-10);
        num = fabs(D[i]-Dgoal[i]);
        if(Von[i]==VOXEL_OAR and D[i]<Dgoal[i]) continue; // 
	//	if(Von[i]==VOXEL_OAR) 
	chi2 += weight[i]*(num*num)/(den*den); // MCTPS
	
//      chi2 += weight[i]*weight[i]*(num*num);  // Lomax
    }
    return chi2;
}

void Optimizer::takeStep(){
    // old serial code >>>
    // for(int32 j=0;j<npb;j++) {

    // float64 Up=0,Down=0;
    // for(int32 i=0;i<numVxl;i++) {
    
    //     if(Von[i]==VOXEL_OFF) continue; // skip this voxel
    //        if(D[i]>0 && Dgoal[i]>0){
    //             float64 square = Dij[i+j*numVxl]*Dij[i+j*numVxl];
    //             int heaviside = 1;
    //             if (Von[i]==VOXEL_OAR and (D[i]<Dgoal[i]))
    //                 heaviside = 0;
    //             Up   += weight[i]*square*Dgoal[i]/D[i]*heaviside; // MCTP
    //             Down += weight[i]*square*heaviside; // MCTP
    //             //    Up   += weight[i]*weight[i]*square*Dgoal[i]/D[i]*heaviside; // Lomax
    //             //    Down += weight[i]*weight[i]*square*heaviside; // Lomax
    //         }
    //     }
    //     if(Down>0){
    //         Nnew[j] *= (1.+fScaleOptimizationStep*(Up/Down-1.));
    //     //  Nnew[j] *= (Up/Down); // Lomax
    //     }
    // }
    // return;
    // old serial code <<<
    vector<pthread_t> threads(pthreads_max_num);
    vector<threadInfo> tInfo(pthreads_max_num);
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for (int i=0; i<pthreads_max_num; i++) {
        tInfo[i].thread_num = i;
        tInfo[i].num_threads = pthreads_max_num;
        tInfo[i].optObj = this;
        pthread_create(&threads[i], &attr, takeStep_thread, (void *) &tInfo[i]);
    }

    /* Wait for all threads to complete */ 
    for (int i=0; i<pthreads_max_num; i++) {
        pthread_join(threads[i], NULL);
    }
}


void Optimizer::run(){
    cout<<"Running the optimizer with the following parameters:"<<endl;
    cout<<'\t'<<"convergence criterium: "<<100*conv<<"%"<<endl;
    cout<<'\t'<<"max iterations: "<<itermax<<endl;
    cout<<'\t'<<"opt. step scaling factor: "<<fScaleOptimizationStep<<endl;
    cout<<endl;

    // buildDgoal();
    // printDgoal();
    selectActivePB();
    
    Nnew = N;

    printPBs();
    buildD();

    buildOptiDose();

    set<int> iterLog={1,2,5,10,20,50,100,200,500,1000,2000,3000,4000,5000};
    ofstream CostVsIter("./out/cost.txt");
    CostVsIter<<"# iter cost"<<endl;

    float64 cost,costPrev;
    int deactivated = 0;


    cout<<"##############   Iterations ##############"<<endl;

    if(f_Debug) {cout<<0<<": P ";  for (int32 i=0; i<numVxl; i++) {cout<<Dgoal[i]<<' ';} cout<<endl;}

    for (iter=0; iter<=itermax; iter++) {

        if(iter>1 and iter%100==0){
            float Nmaxnow = 0;
            for(int i=0;i<npb;i++) Nmaxnow = MAX(Nmaxnow,Nnew[i]);
            cout<<"max N now: "<<Nmaxnow<<endl;
            float thres = 1e-4;
            for(int i=0;i<npb;i++) {
                if(Nnew[i]>0 and Nnew[i]<thres*Nmaxnow){
                    Nnew[i]=0;
                    deactivated++;
                }
            }
            cout<<"total num of deactivated PB: "<<deactivated<<endl;
        }

        if (iter>0) {
            takeStep();
        }

        if(f_Debug)  {cout<<iter<<"  N ";for(int32 j=0;j<npb;j++) {cout<<Nnew[j]<<' ';} cout<<endl;}
        
        buildD();
        
        if(f_Debug) {cout<<iter<<"  D ";  for (int32 i=0; i<numVxl; i++) {cout<<D[i]<<' ';} cout<<endl;}

        costPrev = cost;
        cost = getCost();
        CostVsIter<<iter<<' '<<cost<<endl;

        if (iter<2) {
            printf("%-4d Cost: %12g \n",iter,cost);
        } else {
            float64 costChange=fabs((cost-costPrev)/std::max(1.E-20,costPrev));
            if(costPrev-cost<0.) cout<<"### WARNING: Divergent cost function! ###"<<endl;
            printf("%-4d Cost: %12g Change: %8.4f%%\n",iter,cost,100*costChange);
            if(costChange<conv) break;
        }
        if(iterLog.find(iter)!=iterLog.end()) {
            logFluences(iter);
            // char s[256];
            // sprintf(s,"iter%d",iter);
            // buildOptiDose(s,1);
        }

    }
    cout<<endl;

    cout<<"##########################################"<<endl;

    cout<<endl;
    report();
    cout<<endl;

    N=Nnew;
    printPBs();

}


void Optimizer::buildOptiDose(string out, int wr){

    vector<float> OptiDose(numVxlAll,0);
    vector<int> Vi(numVxlAll);
    vector<float> di(numVxlAll);

    vector<FILE *> files;
    FILE *fin;

    int nn[3];
    float hs[3];
    float x0[3];
    int nspot;

    int jcurr=0;
    for(auto & path : Dij_paths){
        fin = fopen(path.c_str(),"rb");
        if(!fin){
            cerr<<"Error: could not open file: "<<path<<endl;
            exit(1);
        }
        files.push_back(fin);

        fread(nn,3,sizeof(int),fin);
        fread(hs,3,sizeof(float),fin);
        fread(x0,3,sizeof(float),fin);
        fread(&nspot,1,sizeof(int),fin);

        int ipb,ifield,ni;
        for(int j=0;j<nspot;j++){
            fread(&ipb,1,sizeof(int),fin);
            fread(&ifield,1,sizeof(int),fin);
            fread(&ni,1,sizeof(int),fin);
            Vi.resize(ni);
            di.resize(ni);
            // cout<<"\t"<<ipb<<' '<<ifield<<" : "<<ni<<endl;
            fread(Vi.data(),ni,sizeof(int)  ,fin);
            fread(di.data(),ni,sizeof(float),fin);
            for(int i=0;i<ni;i++) {
                OptiDose[Vi[i]]+=di[i]*Nnew[jcurr]*DMF[Vi[i]];
            }
            jcurr++; 
        }
    }

    float l[3]={1,0,0},u[3]={0,1,0},f[3]={0,0,1};
    float L[3];
    for(int i=0;i<3;i++) L[i]=nn[i]*hs[i];
    for(int i=0;i<3;i++) {
        L[i]*=10;
        x0[i]*=10;
        hs[i]*=10;
    }

    if(wr) {
      string fname = "optiDose_" + out + ".mhd";
      write_mhd(fname,OptiDose,nn,x0,L,l,u,f);
      // cout<<"Opti dose = ";
      // for (int iz=0;iz<5;iz++)
      //     cout<<OptiDose[1+3*1+9*iz]<<' ';
      // cout<<endl;
      
      cout<<"Optimized dose map written to "<<fname<<endl;

    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void *takeStep_thread(void *info) 
{
    int tid = ((threadInfo *) info)->thread_num;
    int nthreads = ((threadInfo *) info)->num_threads;
    Optimizer *opt = ((threadInfo *) info)->optObj;

    if(!opt->f_flatOptimizer) {
      //Standard code
      for(int32 j=0;j<opt->npb;j++) {
        if(j%nthreads!=tid) continue; 
	
        float32* column_j = opt->Dij.data()+j*opt->numVxl;
	
        float64 Up=0,Down=0;
        for(int32 i=0;i<opt->numVxl;i++) {
	  
            if(opt->Von[i]==VOXEL_OFF) continue; // skip this voxel
	    
            if(opt->D[i]>0 && opt->Dgoal[i]>0){
	      float32 square = column_j[i]*column_j[i];
	      int heaviside = 1;
	      if (opt->Von[i]==VOXEL_OAR && (opt->D[i]<opt->Dgoal[i]))
                    heaviside = 0;
	      Up   += opt->weight[i]*square*opt->Dgoal[i]/opt->D[i]*heaviside; // MCTP
	      Down += opt->weight[i]*square*heaviside; // MCTP
	      //    Up   += weight[i]*weight[i]*square*Dgoal[i]/D[i]*heaviside; // Lomax
	      //    Down += weight[i]*weight[i]*square*heaviside; // Lomax
            }
        }
        if(Down>0){
	  opt->Nnew[j] *= (1.+opt->fScaleOptimizationStep*(Up/Down-1.));
	  //  Nnew[j] *= (Up/Down); // Lomax
        }
      }
    } else {
      //FLAT beam optimization
      int hmfields = (int)7;
      
      for(int32 ifi=1;ifi<hmfields+1;ifi++) {
	
	//Consider only one field at a time
	float64 Up=0,Down=0;
	for(int32 j=0;j<opt->npb;j++) {
	  
	  if(j%nthreads!=tid) continue; 
	  
	  //The associated field
	  int32 fID = opt->fieldID.at(j);
	  if(fID != ifi) continue;
	  
	  float32* column_j = opt->Dij.data()+j*opt->numVxl;
	  
	  for(int32 i=0;i<opt->numVxl;i++) {
	    
	    if(opt->Von[i]==VOXEL_OFF) continue; // skip this voxel
	    
	    if(opt->D[i]>0 && opt->Dgoal[i]>0){
	      float32 square = column_j[i]*column_j[i];
	      int heaviside = 1;
	      if (opt->Von[i]==VOXEL_OAR && (opt->D[i]<opt->Dgoal[i]))
		heaviside = 0;
	      Up   += opt->weight[i]*square*opt->Dgoal[i]/opt->D[i]*heaviside; // MCTP
	      Down += opt->weight[i]*square*heaviside; // MCTP
	      //    Up   += weight[i]*weight[i]*square*Dgoal[i]/D[i]*heaviside; // Lomax
	      //    Down += weight[i]*weight[i]*square*heaviside; // Lomax
	    }
	  }
	}
	for(int32 j=0;j<opt->npb;j++) {
	  
	  if(j%nthreads!=tid) continue; 
	  
	  //The associated field
	  int32 fID = opt->fieldID.at(j);
	  if(fID != ifi) continue;
	  
	  if(Down>0){
	    opt->Nnew[j] *= (1.+opt->fScaleOptimizationStep*(Up/Down-1.));
	    //  Nnew[j] *= (Up/Down); // Lomax
	  }
	}
	
      }
    }


    
    pthread_exit(NULL);
}

void *buildD_thread(void *info) 
{
    int tid = ((threadInfo *) info)->thread_num;
    int nthreads = ((threadInfo *) info)->num_threads;
    Optimizer *opt = ((threadInfo *) info)->optObj;

    // asch Sep 2020: memory leak when using pthreads and stl vectors!!! ==>
    //vector<float64> myD(opt->numVxl,0); // allocate and reset local dose buffer
    // asch Sep 2020: memory leak when using pthreads and stl vectors!!! <==


    float64 *myD = new float64[opt->numVxl]; // using old-fashioned mem allocation !!!
    for(int32 i=0;i<opt->numVxl;i++) myD[i]=0; // reset local buffer

    // for(int32 i=0;i<opt->numVxl;i++) {
    //     if(i%nthreads!=tid) continue; 

    //     float64* column_i = opt->Dji.data()+i*opt->npb;
        
    //     for(int32 j=0;j<opt->npb;j++) {
    //         opt->D[i]+=column_i[j]*opt->Nnew[j];
    //     }
    // }

    for(int32 j=0;j<opt->npb;j++) {
        if(j%nthreads!=tid) continue; 

        float32* column_j = opt->Dij.data()+j*opt->numVxl;

        float32  Nj = opt->Nnew[j];

        for(int32 i=0;i<opt->numVxl;i++) {
            myD[i] += column_j[i]*Nj;
        }
    }

    pthread_mutex_lock (&buildD_mutex);
    for(int32 i=0;i<opt->numVxl;i++) {
        opt->D[i] += myD[i]; // sum reduce thread pool on global D array
    }
    pthread_mutex_unlock (&buildD_mutex);

    delete [] myD; // release local buffer

    pthread_exit(NULL);
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

void Optimizer::bookHistos(string out){

  char tmpOut[200];
  snprintf(tmpOut,sizeof(tmpOut),"Debug_%s.root",out.data());
  f_out = new TFile(tmpOut,"RECREATE");
  
  h = new TH1D("DoseDiscrepant","Dose discrepant",100,0,1.);
  h = new TH1D("DeltaDiscrepant","Delta for the  discrepant",200,-1,1.);

  h = new TH1D("DoseRatioPBS","Ratio PTV out",100,0,1.);
  h = new TH1D("DoseRatioPBS_1","Ratio PTV out",100,0,1.);
  h = new TH1D("DoseRatioPBS_2","Ratio PTV out",100,0,1.);
  h = new TH1D("PixDoseInside","Ratio PTV out",100,0,5.);
  h = new TH1D("PixDoseOutside","Ratio PTV out",100,0,5.);

  h2 = new TH2D("2DDistsInside","",1700,0,1700,100,0,1.);
  h2 = new TH2D("2DDistOutside","",1700,0,1700,100,0,1.);

  /*
    vp: check histos for the PB dose in ptv and oar
  */
  h = new TH1D("hLogVoxelDose","Voxel dose after filtering",1000,-20,-5);
  h = new TH1D("hLogVoxelDoseNT","Voxel dose  NT",1000,-20,-5);
  h = new TH1D("hLogVoxelDosePTV","Voxel dose in PTV",1000,-20,-5);
  h = new TH1D("hLogVoxelDoseOAR","Voxel dose in OAR",1000,-20,-5);
  h2 = new TH2D("hLogVoxelDoseVSpb","dose voxe vs PB ",150,-20,-3,2001,-0.5,2000.5);
  /*
    h = new TH1D("hLogPbDosePTV","PB dose inside PTV",200,-10,-3.);
    h = new TH1D("hLogPbDoseOAR","PB dose inside OAR",200,-10,-3.);
  */
  
}

