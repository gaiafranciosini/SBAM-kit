

//
//  Optimizer.h
//  Fred
//
//  Created by Agnul on 09/04/13.
//  and modified by S. Pioli
//
//  Copyright (c) 2013 University of Rome La Sapienza. All rights reserved.
//

#ifndef __Fred__Optimizer__
#define __Fred__Optimizer__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

using namespace std;

// type definitions
typedef int         int32;		// 32 bit signed
typedef long long   int64;		// 64 bit signed
typedef float       float32;    // 32 bit floating point (single precision)
typedef double      float64;    // 64 bit floating point (double precision)

class ROI {
public:
	string ID;
	vector<int64> Vi; // voxel indexes
	vector<float64> wei; // voxel weights
	float64 plannedDose,maxDose;
	float DMF; 

	void info();
	void weightRescale(float weight);

};

class Optimizer{
	public:
	int32 npb; // num of pencilbeams or spots
	int32 nroi;  // num of ROIs
	int64 numVxl; // num of voxels used for optimization => reduced grid
	int64 numVxlAll; // num of voxels of original grid => dense grid;
	
	int64 numVxlPTV; // num of voxels beloging to PTV
	int64 numVxlOAR; // num of voxels beloging to OAR
	int64 numVxlDoseCTRL; // num of voxels beloging to Dose Control Volume


	float DMF_NT = 1.0;     // DMF for Normal Tissue
	float maxDose_NT = 1.0; // maxDose for Normal Tissue
	float weight_NT = 10;   // weights for Normal Tissue
	float thre_NT = 5;   // thresholds for NT 
	int Reduction_NT = 1;   // reduction factor for NT 

	float DMF_DCTRL = 1.0;     // DMF for Dose Control Volume
	float maxDose_DCTRL = 1.0; // maxDose for Dose Control Volume
	float weight_DCTRL = 10;   // weight for Dose Control Volume
	float thre_DCTRL = 10;   // threshold for Dose Control Volume (percentage of cumulative dose distrib)
	float Reduction_DCTRL = 1.0;   // reduction factor for Dose Control Volume



	vector<float> DMF; // dense DMF map

	vector<ROI> PTVs;
	vector<ROI> OARs;

	vector<int32> pbID;
	vector<int32> fieldID;

	vector<int32> dense2redIdx; // map of dense to reduced voxel index
	vector<int32> red2denseIdx; // map of reduced to dense voxel index

	int32 gnn[3]; // global
	// float32 gx0[3],ghs[3];

	vector<string> PTV_paths;
	vector<string> OAR_paths;
	vector<string> Dij_paths;
	vector<int> Dij_npb;


	vector<float64> Dgoal,D;
	vector<float32> weight;
	vector<int32> Von;
	vector<float32> Dij; // column-major!!!
	vector<float32> N,Nnew;
	vector<bool> isPBActive;

	int32 iter,itermax; // current and maximum iteration number

	float64 fScaleOptimizationStep;
	float64 conv;
	bool f_Debug;
        bool f_flatOptimizer;
	TH1D *h;
	TH2D *h2;
	TFile *f_out;

	int pthreads_max_num;

	Optimizer() {
		conv = 0.0005;
		fScaleOptimizationStep = 1.0;
		itermax = 10000;
		numVxl=npb=-1;// not set!!!
		numVxlAll=-1;
		numVxlPTV=numVxlOAR=-1;
		f_flatOptimizer = false;
		f_Debug = false;
		pthreads_max_num=6;
	};

	~Optimizer(){};

	void run();
	inline void EnableDebug() {f_Debug = true;}
	int setConv(float64 f);
        int setReductionNT(float64 f);
	int setScaleFactor(float64 f);
	int setIterMax(int Imax);
	float64 getPBPrimaryNum(int jpb);

	// private:

	void loadDij_txt(string fpath);
	void loadDij_bin(vector<string> &fpaths);
  //	void filterDij_bin(string fpath, int call);
	void selectNT_bin(vector<string> &fpaths);
	void selectDoseCtrlVolume_bin(vector<string> &fpaths);
        void loadROI_txt(string fpath);	
	void loadROI_bin(string fpath);

	void init_DMF();

	void loadROI_mhd(ROI &roi,string fpath);
	void loadPTV_mhd(string ID,string fpath,float plannedDose,float weight, float DMF);
	void loadOAR_mhd(string ID,string fpath,float maxDose,float weight, float DMF);
	
	void buildDgoal();
	void selectActivePB();
	void MonitorDoseRatios();

	void initParticleNumbers(string start);
	void makeDji();

	void buildD();
	void takeStep();
	float64 getCost();

	void report();
	void bookHistos(string out="def");
	void saveOptimizedPB();
	void printD();
	void printDgoal();
	void printPBs();

    void buildOptiDose(string out="def", int wr = 0);
	void writePBs(string out="def");

	void logFluences(int iter);

};

#endif /* defined(__Fred__Optimizer__) */

