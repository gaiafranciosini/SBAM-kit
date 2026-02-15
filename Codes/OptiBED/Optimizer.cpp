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
    for(int i=0;i<(int)hist.size();i++) hist[i]/=sum;
    // get cumulative
    for(int i=1;i<(int)hist.size();i++) hist[i]+=hist[i-1];
    // find lower cutoff for given threshold
    int icut(0);
    for(int i=0;i<(int)hist.size();i++) {
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
    cout<<"\tFMFmin = "<<FMFmin<<endl;
    cout<<"\tAlpha/Beta = " << alphabeta << endl;
    cout<<"\tDthr="<<Dthr<<endl;

    cout<<endl;
}

void ROI::weightRescale(float weight){
    float ww = weight/wei.size();  // asch Apr 2023 volume-normalized weight
    for(size_t i=0;i<wei.size();i++) wei[i]*=ww;
    //wei.size() = numero di voxel attivi dell'oar/ptv considerato 
}


//////////////////////////////////////////////////////////////////////////////////////////
int Optimizer::numFractions;

int Optimizer::setReductionNT(float64 f){if((f<=0)&&(f>1)) return -1; Reduction_NT= f; return 0;}
int Optimizer::setConv(float64 f){if(f<=0) return -1; conv= f; return 0;}
int Optimizer::setScaleFactor(float64  f){if(f<=0) return -1; fScaleOptimizationStep=f; return 0;}
int Optimizer::setIterMax(int Imax){if(Imax<=0) return -1; itermax=Imax; return 0;}

void Optimizer::setNumFractions(int value) {
    numFractions = value;
}
int Optimizer::getNumFractions() {
  return numFractions;
}


void Optimizer::loadROI_txt(string fpath){ //non usato, diamo input con mhd
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
  cout << "1) loadROI_mhd called =======================================================================" << endl;
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

    //controlla che tutte le ROI appartengano ad una grid con le stesse dimensioni
    if(numVxlAll<0) { // not set => first ROI
        numVxlAll = 1UL*nn[0]*nn[1]*nn[2];
        for(int i=0;i<3;i++) {gnn[i]=nn[i];}
    } else {   // check that geometry is aligned with other ROIs
      if(numVxlAll != (int64)1UL*nn[0]*nn[1]*nn[2]) {
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
    
    //alphabeta.resize(4, vector<float>(nonzero,1000000.0)); //qui 26 mar

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


void Optimizer::loadPTV_mhd(string ID,string fpath,float plannedDose,float weight, float FMFmin, float alpha_beta, float Dthr){
   cout << "2) loadPTV_mhd called =======================================================================" << endl;
    ROI ptv;
    ptv.ID = ID;
    ptv.plannedDose = plannedDose;
    ptv.maxDose = plannedDose;
    ptv.FMFmin = FMFmin;
    ptv.alphabeta = alpha_beta;
    ptv.Dthr=Dthr;
    loadROI_mhd(ptv,fpath);
    ptv.weightRescale(weight);
    PTVs.push_back(ptv);
    PTV_paths.push_back(fpath);

}

void Optimizer::loadOAR_mhd(string ID,string fpath,float maxDose,float weight, float FMFmin, float alpha_beta,float Dthr){
    cout << "3) loadOAR_mhd called =======================================================================" << endl;
    ROI oar;
    oar.ID = ID;
    oar.plannedDose = 0;
    oar.maxDose = maxDose;
    oar.FMFmin = FMFmin;
    oar.alphabeta = alpha_beta;
    oar.Dthr=Dthr;
    loadROI_mhd(oar,fpath);
    oar.weightRescale(weight);
    OARs.push_back(oar);
    OAR_paths.push_back(fpath);

}


void Optimizer::loadDij_txt(string fpath){
  cout << "4) loadDij_txt called =======================================================================" << endl;
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

    int NF = getNumFractions();
    
    D.resize(NF); //dose D
    for(auto &vec : D){
      vec.resize(numVxl);
    }
    initParticleNumbers("ones");
}



void Optimizer::selectNT_bin(vector<string> &fpaths){
    cout << "5) selectNT_bin called =======================================================================" << endl;

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
    cout<<"NT voxel Over Thr= "<<numVxlOverThrNT<<" NT voxel selected= "<<numVxlSelectNT<< " num crazy voxel= "<<num_crazy_VXL<<endl;
    cout<<"NumVxl after NoT= "<<numVxl<<"  numVxl OAR= "<<numVxl<<endl;
    cout<<"Check exiting selectNT_bin : red2denseIdx.size()= "<<red2denseIdx.size()<< endl;
}

void Optimizer::selectDoseCtrlVolume_bin(vector<string> &fpaths){
  cout << "6) selectDoseCtrlVolume_bin called =======================================================================" << endl;
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
            for(int i=0;i<(int)Vi.size();i++){
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
    for(int i=0;i<(int)DOSECTRL.size();i++) {
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
    oar.FMFmin = FMFmin_DCTRL;
    oar.alphabeta = alphabeta_DCTRL;
    oar.Dthr=Dthr_DCTRL;
    //qui 9 aprile
    // loadROI_mhd(oar,fpath);

    vector<int> vv;
    vv.reserve(numVxlDoseCTRL);
    for(int i=0;i<(int)DOSECTRL.size();i++) {
        if(DOSECTRL[i]>0) vv.push_back(i);
    }

    oar.Vi.resize(numVxlToSelect,-1);
    // for(int i=0;i<numVxlToSelect;i++) oar.Vi[i]=vv[i];
    oar.wei.resize(numVxlToSelect,1);
    /*
    oar.alphabeta.resize(4);
    for(int iField = 0; iField < 4; iField++){
      oar.alphabeta.at(iField) = 1000000.0;
    } //qui 26 mar
    */

    cout<<"random select voxels in original Dose Control Volume"<<endl;
    // random select voxels in original Dose Control Volume
    //unsigned int seed = 181527;
    //unsigned int seed = 181528;
    //srand(seed);
    for(int i=0;i<numVxlToSelect;i++){
        float r = 1.0*rand()/RAND_MAX;
	//cout << "random selection of voxels: r = " << r << endl;
        int idx = int(r*vv.size());
        oar.Vi[i] = vv[idx];
        vv.erase(vv.begin()+idx);
    }

    oar.weightRescale(weight_DCTRL);
    OARs.push_back(oar);
    OAR_paths.push_back("out/DoseCTRL.mhd");

    oar.info();

    // add OAR to reduced matrices 
    for(int i=0;i<(int)oar.Vi.size();i++){
        dense2redIdx[oar.Vi[i]] = numVxl++;
	//cout << "DCTRL dense2redIdx[oar.Vi[i]] = " << dense2redIdx[oar.Vi[i]] << endl;
        numVxlOAR++;
        red2denseIdx.push_back(oar.Vi[i]);
	//cout << "DCTRL red2denseIdx[i] = " << oar.Vi[i] << endl;
        Dgoal.push_back(oar.maxDose);
	alphabeta.push_back(alphabeta_DCTRL);
	//Dthr.push_back(Dthr_DCTRL);
        weight.push_back(oar.wei[i]);
        Von.push_back(VOXEL_OAR);
    }
    cout<<"NumVxl after DoseCTRL= "<<numVxl<<"  numVxl OAR= "<<numVxlOAR<<endl;
}

void Optimizer::loadDij_bin(vector<string> &fpaths){
    cout << "7) loadDij_bin called =======================================================================" << endl;
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
    //vector<float> FMF(numVxlAll);//MODIFICA
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
		  //if(di[i]<Dthr[ired]){//MODIFICA
		  //  FMF[Vi[i]]=1;
		  //} else {
		  //  FMF[Vi[i]]=(1-DMF[Vi[i]])*(Dthr[ired]/di[i])+DMF[Vi[i]];
		    // }
		  //Dij[ired+jcurr*numVxl]=di[i]*FMF[Vi[i]]; // set voxel dose applying DMF
		  Dij[ired+jcurr*numVxl]=di[i]; 
                }
            }      
            jcurr++; // update current global spot index
        }

    }    
    
    int NF = getNumFractions();
    
    D.resize(NF); //dose D
    for(auto &vec : D){
      vec.resize(numVxl);
    }
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
      cout << "8) init_DMF called =======================================================================" << endl;
    FMFmin.clear(); FMFmin.resize(numVxlAll,FMFmin_DCTRL); // allocate dense map for DMF setting it to normal tissue DMF
    Dthr.clear(); Dthr.resize(numVxlAll,Dthr_DCTRL);

    
    // update FMFmin Dthr for OAR(s)
    for(auto & roi : OARs) {
      for(int k=0;k<(int)roi.Vi.size();k++) {
	FMFmin[roi.Vi[k]]=roi.FMFmin;
	Dthr[roi.Vi[k]]=roi.Dthr;
        }
    }

    // PTV should be the last to be updated.    
    // update DMF for PTV(s)
    for(auto & roi : PTVs) {
      for(int k=0;k<(int)roi.Vi.size();k++) {
            FMFmin[roi.Vi[k]]=roi.FMFmin;
	    Dthr[roi.Vi[k]]=roi.Dthr;
        }
    } 
}

double Optimizer::compute_FMF(double d, int idense){
  double FMF = 1;
  
  if(d > Dthr[idense]){
    FMF = (1 - FMFmin[idense]) * (Dthr[idense]/d) + FMFmin[idense];
    /*if (vxl == 397483) {
      cout << "vxl " << vxl << " (red), " << i << " nella dense" << endl;
      cout << "FMFmin = " << FMFmin[i] << endl;
      cout << "Dthr = " << Dthr[i] << endl;
      cout << "d = " << d << endl;
      cout << "FMF = " << FMF << endl;
      }*/
    //if() cout << "inside compute_FMF with d > Dthr for vxl " << vxl << "(red) " << i << " (dense)" << endl;
  }
    return FMF;
}


void Optimizer::buildDgoal(){
    cout << "9) buildDgoal called =======================================================================" << endl;
    //int num_fractions;
   
  numVxlPTV=numVxlOAR=0; 
  for(auto roi : PTVs) numVxlPTV+=roi.Vi.size();
  for(auto roi : OARs) numVxlOAR+=roi.Vi.size();
  
  numVxl = numVxlPTV+numVxlOAR;
  cout << "inizio di builDgoal " << numVxl << endl << endl;
  
  Dgoal.clear(); Dgoal.resize(numVxl,0);
  alphabeta.clear();
  alphabeta.resize(numVxl, 0);
  //Dthr.clear();
  //Dthr.resize(numVxl,0);
  Von.clear(); Von.resize(numVxl,0); // voxel label 0=off, 1=on (PTV), 2=on (OAR)
  weight.clear(); weight.resize(numVxl,1);
  dense2redIdx.clear(); dense2redIdx.resize(numVxlAll,-1);
  red2denseIdx.clear(); red2denseIdx.resize(numVxl,-1);
 
  int32 i = 0;
  int32 iDou = 0;
  int32 iDouPTV = 0;
  for(auto roi : PTVs) {
    for(int k=0;k<(int)roi.Vi.size();k++) {
      if (dense2redIdx[roi.Vi[k]] != -1)   { // already assigned
      	iDouPTV++;
	iDou++;
     	numVxl--;
      } else {
	dense2redIdx[roi.Vi[k]] = i;
	red2denseIdx[i] = roi.Vi[k];
	//cout << "roi: " << roi.ID << endl;
	//cout << "dense2redIdx[roi.Vi[k]] = " << dense2redIdx[roi.Vi[k]] << endl;
	//cout << "red2denseIdx[i] = " << red2denseIdx[i] << endl;
	Dgoal[i]=roi.plannedDose;
	alphabeta[i]=roi.alphabeta;
	weight[i]=roi.wei[k];
	Von[i]=VOXEL_PTV;
	i++; //voglio scorrere nella matrice ridotta solo quando ttrovo un vxl non ancora assegnato
      }
    }
  }


  for(auto roi : OARs) {
    for(int k=0;k<(int)roi.Vi.size();k++) {
      if (dense2redIdx[roi.Vi[k]] != -1)   { // already assigned to PTV
	//PTV has the priority.
	iDou++;
	numVxl--;
      } else {
	//before assignining...
	dense2redIdx[roi.Vi[k]] = i;
	red2denseIdx[i] = roi.Vi[k];
	Dgoal[i]=roi.maxDose;
	alphabeta[i]=roi.alphabeta;
	//cout << "roi: " << roi.ID << endl;
	//cout << "dense2redIdx[roi.Vi[k]] = " << dense2redIdx[roi.Vi[k]] << endl;
	//cout << "red2denseIdx[i] = " << red2denseIdx[i] << endl;
	weight[i]=roi.wei[k];
	Von[i]=VOXEL_OAR;
	i++;
      }
    }
  }
  
  /*
  for(size_t iField = 0; iField < alphabeta.size(); iField++){
    for(int i_vxl = 0; i_vxl < numVxl; i_vxl++){
      if(alphabeta.at(iField)[i_vxl] == 0.0){
	cout << "alphabeta at field " << iField << " in voxel " << i_vxl << " = " << alphabeta.at(iField)[i_vxl] << endl;
      }
    }
    }*/


  red2denseIdx.resize(numVxl);
  Von.resize(numVxl);
  Dgoal.resize(numVxl);
  alphabeta.resize(numVxl);
  //Dthr.resize(numVxl);
  weight.resize(numVxl);
  cout<<iDouPTV<<" voxels in common among the PTVs"<<endl;
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
	    //delta = fabs((D[i]-Dgoal[i])/Dgoal[i]); //dose D
	    
	    // Sum the doses over all fractions
            double Dsum = 0.0;
            for (size_t fid = 0; fid < D.size(); ++fid) {
	      Dsum += D.at(fid)[i];
            }

            delta = fabs((Dsum - Dgoal[i]) / Dgoal[i]);

            

	    ((TH1D*)gDirectory->Get("DoseDiscrepant"))->Fill(Dgoal[i]);
	    ((TH1D*)gDirectory->Get("DeltaDiscrepant"))->Fill(delta);
	    
            s1 += delta;
            deltamax = std::max(deltamax,Dsum/Dgoal[i]); //dose D
            deltamin = std::min(deltamin,Dsum/Dgoal[i]); //dose D
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

    //for(int i=0;i<numVxl;i++) D[i]=0; // reset dose //dose D

  int NF = getNumFractions();

  vector<vector<float64>> myBED;

  myBED.resize(NF);
  for(auto &vec : myBED){
    vec.resize(numVxl);
  }

  for (int f = 0; f < NF; ++f) {
    for (int i = 0; i < numVxl; ++i) {
      D.at(f)[i] = 0.0;
      myBED.at(f)[i] = 0.0;
    }
  }

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

    for(int i = 0; i < numVxl; i++){
      for(int f = 0; f < NF; f++){
	int idense = red2denseIdx[i];
	float64 d = D.at(f)[i];
	float FMF = compute_FMF(d,idense);
	myBED.at(f)[i] = d * FMF * (1.0 + (d * FMF)/alphabeta[i]);
      }
    }
    /*
    int dense12 = red2denseIdx[12];
    cout << "vxl 12 (red) / " << dense12 << " (dense)" << endl;
    cout << "d_0 = " << D.at(0)[12] << "d_1 = " << D.at(1)[12] << "d_2 = " << D.at(2)[12] << "d_3 = " << D.at(3)[12] << endl;
    cout << "alphabeta[12] = " << alphabeta[12] << endl;
    float myBED_0 = myBED.at(0)[12] + myBED.at(1)[12] + myBED.at(2)[12] + myBED.at(3)[12];
    cout << "BED[12] = " << myBED_0 << endl;
    cout << "FMF = " << compute_FMF(D.at(2)[12], dense12) << endl;
    */



     for(int i = 0; i < numVxl; i++){
      for(int f = 0; f < NF; f++){
	D.at(f)[i] = myBED.at(f)[i];
      }
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
    for (size_t f = 0; f < D.size(); ++f) { 
      for (int i = 0; i < numVxl; ++i) { 
	cout<<"\t"<<D.at(f)[i]<<endl; //dose D
      }
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
    int NF = getNumFractions();
    //vector<float> chi2fx;
    //chi2fx.resize(NF,0.0);
    float Dsum = 0.0;
    for (int32 i=0; i<numVxl; i++) {
      if(Von[i]==VOXEL_OFF) continue; // skip this voxel
      for(int iFx = 0; iFx < NF; iFx++){
	Dsum += D.at(iFx)[i];
      }
      den = std::max(Dgoal[i], 1e-10);
      num = fabs(Dsum-Dgoal[i]); //dose D
      if(Von[i]==VOXEL_OAR and Dsum<Dgoal[i]) continue; //dose D
      //	if(Von[i]==VOXEL_OAR) 
      chi2 += weight[i]*(num*num)/(den*den); // MCTPS
      //chi2 += weight[i]*weight[i]*(num*num);  // Lomax
      Dsum = 0.0;
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
    // cout << "sono dopo printPBs e prima di buildD =================================" << endl;
    buildD();
    //cout << "sono dopo buildD e prima di buildOptiDose =================================" << endl;

    buildOptiDose();

    set<int> iterLog={1,2,5,10,20,50,100,200,500,1000,2000,3000,4000,5000};
    ofstream CostVsIter("./out/cost.txt");
    CostVsIter<<"# iter cost"<<endl;

    float64 cost(1.e9),costPrev;
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
        
        if(f_Debug) {cout<<iter<<"  D ";  for (int32 i=0; i<numVxl; i++) {cout<<D.at(1)[i]<<' ';} cout<<endl;} //dose D

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

  int NF = getNumFractions();

  vector<vector<float>> OptiDose2;
  OptiDose2.resize(NF);
  for(auto &vec : OptiDose2) {
    vec.resize(numVxlAll,0.0);
    }
  vector<float> OptiDose(numVxlAll,0);
  vector<float> OptiEQD(numVxlAll,0);
  vector<float> OptiBED(numVxlAll,0);
  vector<float> OptiDoseFLASH(numVxlAll,0);
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
            fread(Vi.data(),ni,sizeof(int)  ,fin);
            fread(di.data(),ni,sizeof(float),fin);
	    for (int i = 0; i < ni; i++) {
	      // int ison = dense2redIdx.at(i); //MODIFICA 9 giugno
	      /*
	      if(di[i]<Dthr[ison]){//MODIFICA
		    FMF[Vi[i]]=1;
		  } else {
		    FMF[Vi[i]]=(1-DMF[Vi[i]])*(Dthr[ison]/di[i])+DMF[Vi[i]];
		    } */
	      // OptiDose2.at(ifield-1)[Vi[i]]+=di[i]*Nnew[jcurr]*FMF[Vi[i]]; MODIFICA
	      OptiDose2.at(ifield-1)[Vi[i]]+=di[i]*Nnew[jcurr];
	      //OptiDose[Vi[i]]+=di[i]*Nnew[jcurr]*FMF[Vi[i]];
	      OptiDose[Vi[i]]+=di[i]*Nnew[jcurr];
	    }
            jcurr++;
	}
    }
    

    // cout << "OptiDose2[19864911]" << OptiDose2.at(0)[19864911] << " " << OptiDose2.at(1)[19864911] << " " << OptiDose2.at(2)[19864911] << " " << OptiDose2.at(3)[19864911] << endl;
    
  
    float l[3]={1,0,0},u[3]={0,1,0},f[3]={0,0,1};
    float L[3];
    for(int i=0;i<3;i++) L[i]=nn[i]*hs[i];
    for(int i=0;i<3;i++) {
        L[i]*=10;
        x0[i]*=10;
        hs[i]*=10;
    }

    /*
    for(int i=0;i<all_Vi.size();i++) {
      int ison = dense2redIdx.at(all_Vi[i]);
      for(int iFx = 0; iFx < NF; iFx++){
	double diFx = OptiDose2.at(iFx)[all_Vi[i]];
       	double ab = alphabeta_DCTRL; //conversione in EQD dei voxel non appartenenti alle roi inserite dal .inp
	//ison == -1 per i vxl che non appartengono a ptv, oar messi nel .inp e doseCTRLvol
	if(ison!=-1) ab = alphabeta[ison];
	OptiEQD[all_Vi[i]] += diFx * ((diFx + ab) / (15.0 + ab));
      }
    }
    */
    
    for(int i = 0; i<numVxlAll; i++){
      for(int iFx = 0; iFx < NF; iFx++){
	//int iRed = dense2redIdx[i];
	double diFx = OptiDose2.at(iFx)[i];
	double FMF = compute_FMF(diFx,i);
	OptiDoseFLASH[i] += diFx * FMF;
      }
      }
    
    for(int i=0;i<numVxlAll;i++) {
      //cout << i << " --- " <<  dense2redIdx[i] << endl;
      //int iRed = dense2redIdx[i];
      int ison = dense2redIdx.at(i);
      for(int iFx = 0; iFx < NF; iFx++){
	double diFx = OptiDose2.at(iFx)[i];
	double FMF = compute_FMF(diFx,i);
	double ab = alphabeta_DCTRL;
	if(ison!=-1) ab = alphabeta[ison];
	//OptiEQD[i]+=diFx  * ((diFx + ab) / (2 + ab)); //non flash
	//OptiBED[i]+=diFx  * (1.0 + (diFx/ab)); //non flash
	OptiEQD[i] += diFx * FMF * (((diFx * FMF) + ab) / (2 + ab)); //flash
	OptiBED[i] += diFx * FMF * (1.0 + ((diFx * FMF)/ab)); //flash
      }
    }
    
    // cout << "OptiEQD2[19864911] = " << OptiEQD[19864911] << endl;

    if(wr) {
      string fname = "optiDose_" + out + ".mhd";
      // cout<<"Opti dose = ";
      // for (int iz=0;iz<5;iz++)
      //     cout<<OptiDose[1+3*1+9*iz]<<' ';
      // cout<<endl;
      cout<<"Optimized dose map written to "<<fname<<endl;
    }
    
    if(wr) {
      string fname = "optiDoseFLASH_" + out + ".mhd";
      write_mhd(fname,OptiDoseFLASH,nn,x0,L,l,u,f);
      // cout<<"Opti dose = ";
      // for (int iz=0;iz<5;iz++)
      //     cout<<OptiDose[1+3*1+9*iz]<<' ';
      // cout<<endl;
      cout<<"Optimized dose map written to "<<fname<<endl;
      }

 
    if(wr) {
      string fname = "optiBED_" + out + ".mhd";
      write_mhd(fname,OptiBED,nn,x0,L,l,u,f);
      // cout<<"Opti dose = ";
      // for (int iz=0;iz<5;iz++)
      //     cout<<OptiDose[1+3*1+9*iz]<<' ';
      // cout<<endl;
      cout<<"Optimized BED map written to "<<fname<<endl;
    }

    if(wr) {
      string fname = "optiEQD2_" + out + ".mhd";
      write_mhd(fname,OptiEQD,nn,x0,L,l,u,f);
      // cout<<"Opti dose = ";
      // for (int iz=0;iz<5;iz++)
      //     cout<<OptiDose[1+3*1+9*iz]<<' ';
      // cout<<endl;
      cout<<"Optimized EQD2 map written to "<<fname<<endl;
    }
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

void *takeStep_thread(void *info) 
{
    int tid = ((threadInfo *) info)->thread_num;
    int nthreads = ((threadInfo *) info)->num_threads;
    Optimizer *opt = ((threadInfo *) info)->optObj;
    int NF = opt->numFractions;
    float Dsum = 0.0;

    // if(!opt->f_flatOptimizer) {
    //cout << "SONO QUI ==================================================================" << endl;
      //Standard code
      for(int32 j=0;j<opt->npb;j++) {
        if(j%nthreads!=tid) continue;
	
        float32* column_j = opt->Dij.data()+j*opt->numVxl;
	//int32 fID = opt->fieldID.at(j);
	
        float64 Up=0,Down=0;
        for(int32 i=0;i<opt->numVxl;i++) {
	  for(int iFx = 0; iFx < NF; iFx++){
	    Dsum += opt->D.at(iFx)[i];    
	  }
	  
	  //float32 phys_D = column_j[i];
	  //float32 bio_D = phys_D * (1.0f + phys_D / opt->alphabeta[i]);
	  if(opt->Von[i]==VOXEL_OFF) continue; // skip this voxel
	    
	  if(Dsum>0 && opt->Dgoal[i]>0){ //dose D
	      //float32 square = bio_D * bio_D;
	    float32 square = column_j[i]*column_j[i];
	    int heaviside = 1;
	    
	    if (opt->Von[i]==VOXEL_OAR && (Dsum<opt->Dgoal[i])){ //dose D
		heaviside = 0;
	    }
	    Up   += opt->weight[i]*square*(opt->Dgoal[i])/Dsum*heaviside; // MCTP //dose D
	    Down += opt->weight[i]*square*heaviside; // MCTP
	    //    Up   += weight[i]*weight[i]*square*Dgoal[i]/D[i]*heaviside; // Lomax
	    //    Down += weight[i]*weight[i]*square*heaviside; // Lomax
	  }
	  Dsum = 0.0;
        }
	//cout << "opt step" << opt->fScaleOptimizationStep << endl;
        if(Down>0){
	  opt->Nnew[j] *= (1.+opt->fScaleOptimizationStep*(Up/Down-1.));
	  //  Nnew[j] *= (Up/Down); // Lomax
        }
      }
      /*} else {
      //FLAT beam optimization
      int hmfields = (int)7;
      cout << "SONO QUI ==================================================================" << endl;
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
	    float32 phys_D = column_j[i];
	    float32 bio_D = phys_D * (1.0f + phys_D / opt->alphabeta[i]);
	    if(opt->Von[i]==VOXEL_OFF) continue; // skip this voxel
	    
	    if(opt->D.at(fID-1)[i]>0 && opt->Dgoal[i]/NF>0){ //dose D
	      float32 square = bio_D*bio_D;
	      int heaviside = 1;
	      if (opt->Von[i]==VOXEL_OAR && (opt->D.at(fID-1)[i]<opt->Dgoal[i]/NF)) //dose D
		heaviside = 0;
	      Up   += opt->weight[i]*square*(opt->Dgoal[i]/NF)/opt->D.at(fID-1)[i]*heaviside; // MCTP //dose D
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
      */

    
    pthread_exit(NULL);
}


void *buildD_thread(void *info) 
{
    int tid = ((threadInfo *) info)->thread_num;
    int nthreads = ((threadInfo *) info)->num_threads;
    Optimizer *opt = ((threadInfo *) info)->optObj;

    int NF = opt->numFractions;
    int numVxl = opt->numVxl;

    float64 **myDose = new float64*[NF];
    //float64 **myBED = new float64*[NF];
    for (int f = 0; f < NF; ++f) {
        myDose[f] = new float64[numVxl];
	//myBED[f] = new float64[numVxl];
        for (int i = 0; i < numVxl; ++i) {
            myDose[f][i] = 0.0;
	    //myBED[f][i] = 0.0;
	}
    }

    // Loop su tutti i pencil beam
    for (int32 j = 0; j < opt->npb; ++j) {
      if (j % nthreads != tid)	continue; //skippiamo i pb non assegnati a questo thread
      
      //cout << "PB j = " << j << " -- tid = " << tid << endl;
        int fieldID = opt->fieldID.at(j); // fieldID ∈ [1, NF]
        int f = fieldID - 1;              // convertiamo in [0, NF-1]
	
        float32* column_j = opt->Dij.data() + j * numVxl;
        float32 Nj = opt->Nnew[j];

        for (int i = 0; i < numVxl; ++i) {
            myDose[f][i] += column_j[i] * Nj; // somma la dose fisica per field
        }
    }
    
    /* for (int f = 0; f < NF; ++f) {
        for (int i = 0; i < numVxl; ++i) {
            float64 d = myDose[f][i];
            float64 ab = opt->alphabeta[i];
            if (ab > 0.0)
                myBED[f][i] = d * (1.0 + d / ab);
            else {
                myBED[f][i] = 0.0;
                cout << "Warning: alphabeta = 0 nel voxel " << i << endl;
            }
        }
	}*/

    // Scrittura thread-safe su opt->D[f][i]
    pthread_mutex_lock(&buildD_mutex);
    for (int f = 0; f < NF; ++f) {
        for (int i = 0; i < numVxl; ++i) {
	  //opt->D.at(f)[i] += myBED[f][i];
	  opt->D.at(f)[i] += myDose[f][i];
        }
    }
    pthread_mutex_unlock(&buildD_mutex);

    for (int f = 0; f < NF; ++f) {
        delete[] myDose[f];
        //delete[] myBED[f];
    }
    delete[] myDose;
    //delete[] myBED;

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

