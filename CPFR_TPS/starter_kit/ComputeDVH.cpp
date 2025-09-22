#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <set>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <ostream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <bits/stdc++.h>
#include <ctime>
#include "mhd_io.h"

using namespace std;

std::vector<std::string> args;

double closest(std::vector<double> const& vec, double value) {
    auto const it = std::lower_bound(vec.begin(), vec.end(), value);
    if (it == vec.begin()) { return -1; }
    else
      return *(it);
    }
   
int findIndex(const vector<double> &arr, double item) {

    for (auto i = 0; i < arr.size(); ++i) {
        if (arr[i] == item)
            return i;
    }
     return -1;
}

int main (int argc, char *argv[]){
  std::vector<std::string> filename;
  string dosefile;
  string roitype;
  string fileLabel;
  string dir;
  float Dgoal = 0.0f;
  int nn[3];
  float x0[3];
  float L[3];
  float hs[3];
  float left[3];
  float up[3];
  float front[3];
  int binningDVH = 100;
  bool plot = false;
  string allROINames;
    

    if (argc == 2 && std::string(argv[1]) == "-h") {
        std::cout << "Usage: " << argv[0] << " -Dgoal <value> -roi <filename1> <filename2> ... -dose <filename> -type <int> or <float> -fileLabel <value> -dir <outputdirectory> -plot or not if you want to display the DVH" << std::endl;
        return 0; 
    }

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " -Dgoal <value> -roi <filename1> <filename2> ... -dose <filename> -type <int> or <float> -fileLabel <value> -dir <outputdirectory> -plot or not if you want to display the DVH" << std::endl;
        return 1;
    }

    for (int i = 1; i < argc - 1; ++i) {

      if (std::string(argv[i]) == "-Dgoal") {
	if (i < argc - 1) {
	  Dgoal = std::stof(argv[i + 1]);
	  ++i;
	} else {
	  std::cerr << "Missing value for -Dgoal option" << std::endl;
	  return 1;
	}
      }else if(std::string(argv[i]) == "-plot") {
	plot = true;
      }else if (std::string(argv[i]) == "-roi") {
	while (i + 1 < argc && std::string(argv[i + 1]) != "-Dgoal" && std::string(argv[i + 1]) != "-dose") {
	  filename.push_back(argv[++i]);
	}
      } else if (std::string(argv[i]) == "-dose") {
	if (i < argc - 1) {
	  dosefile = std::string(argv[i + 1]);
	  ++i;
	} else {
	  std::cerr << "Missing value for -dose option" << std::endl;
	  return 1;
	}
      } else if (std::string(argv[i]) == "-type") {
	if (i < argc - 1) {
	  roitype = std::string(argv[i + 1]);
	  ++i;
	} else {
	  std::cerr << "Missing value for -type option" << std::endl;
	  return 1;
	}
      } else if(string(argv[i]) == "-dir"){
	if(i<argc-1){
	  dir = string(argv[i+1]);
	  ++i;
	}else{
	  std::cerr << "Missing value for -dir option" << std::endl;
	}
      } else if (std::string(argv[i]) == "-fileLabel") {
	if (i < argc - 1) {
	  fileLabel = std::string(argv[i + 1]);
	  ++i;
	} else {
	  std::cerr << "Missing value for -fileLabel option" << std::endl;
	  return 1;
	}
      }else {
	std::cerr << "Unknown option: " << argv[i] << std::endl;
	return 1;
      }
    }
    
    if (Dgoal == 0.0f || filename.empty() || dosefile.empty() || roitype.empty() || fileLabel.empty()) {
      std::cerr << "Missing required options or the input line has wrong options" << std::endl;
      return 1;
    }
    
    vector <float32> dose;
    cout<<"Dose file: "<<dosefile<<endl;
    
    read_mhd(dosefile, dose, nn, x0, L, left, up, front);
    
    for(int i = 0; i < 3; i++){
      L[i]/=10.;
      x0[i]/=10.;
    }
    
    // -------- loop over the ROIs and DVHs calculation ------------//
    
    for(int i = 0 ; i<filename.size();i++){

      vector <int16> roi_int;
      vector <float32> roi_float;
      
      int nVoxelROI = 0;
      float voxelVolume = 0.0f;
      float ROIVolume = 0.0f;
      
      
      cout<<"roitype: "<<roitype<<endl;
      
      if(roitype == "int"){
	
	read_mhd(filename.at(i), roi_int, nn, x0, L, left, up, front);
	size_t lastSlash = filename.at(i).find_last_of("/");
	size_t lastDot = filename.at(i).find_last_of(".");
	std::string ROI = filename.at(i).substr(lastSlash + 1, lastDot - lastSlash - 1);
	allROINames += ROI;
	if (i != filename.size() - 1) {
	  allROINames += " ";
	}
	//std::string ROI = filename.at(i).substr(lastSlash + 1);
	std::vector<double> D(roi_int.size(),0);
	for(int i = 0; i < 3; i++){
	  L[i]/=10.;
	  x0[i]/=10.;
	  hs[i]=L[i]/nn[i];
	}

     voxelVolume = hs[0]*hs[1]*hs[2];

     cout<<"VoxelVolume: "<<voxelVolume<<endl;
     
      for(int n = 0; n < dose.size();n++){
	if(roi_int[n]>0){
	  ROIVolume += roi_int[n]*voxelVolume;
	  nVoxelROI++;
	  D[n]=(dose[n]*100);
	}
       }

     cout<<"ROI: "<<ROI<<" Non-zero voxels: "<<nVoxelROI<< " VolumeRoi: "<<ROIVolume<<endl;
     
     float maxDose = *std::max_element(D.begin(), D.end());
     cout<<"MaxDose: "<<maxDose<<endl; 
     float binSize = Dgoal/binningDVH;
     int nbin = std::max(binningDVH, static_cast<int>(maxDose / binSize) + 1);
     cout<<maxDose<<" "<<binSize<< " "<<nbin<<" "<<Dgoal<<endl;
     std::vector<float> bins(nbin, 0.0);
     std::vector<float> hist(nbin, 0.0);

     // Calcolo di bins
     for (int i = 0; i < nbin; i++) {
       bins[i] = (i + 0.5) * binSize;
    }
     
    // Calcolo di hist
     for (size_t i = 0; i < D.size(); ++i) {
       int imax = static_cast<int>(D[i] / binSize);
       for (int j = 0; j <= imax; ++j) {
	 hist[j] += ((1*roi_int[i])*voxelVolume);
       }
     }
     //string ROIDVH = ROI + fileLabel + ".txt";
     string ROIDVH = dir+ "/" +ROI + fileLabel + ".txt";
     ofstream fout (ROIDVH);
     for (int i = 0; i<=nbin -1; i++){
       fout << bins[i]<<" "<< hist[i]/ROIVolume*100.<< endl;
     }

     
   }else{

	
	read_mhd(filename.at(i), roi_float, nn, x0, L, left, up, front);
	/*	for(int n = 0; n < roi_float.size();n++){
	  if(roi_float[n]>0)cout<<roi_float[n]<<endl;
	  }*/
	size_t lastSlash = filename.at(i).find_last_of("/");
	size_t lastDot = filename.at(i).find_last_of(".");
	std::string ROI = filename.at(i).substr(lastSlash + 1, lastDot - lastSlash - 1);
	allROINames += ROI;
	if (i != filename.size() - 1) {
	  allROINames += " ";
	}
	//std::string ROI = filename.at(i).substr(lastSlash + 1);
	std::vector<double> D(roi_float.size(),0);
	cout<<roi_float.size()<<endl;
	for(int i = 0; i < 3; i++){
	  L[i]/=10.;
	  x0[i]/=10.;
	  hs[i]=L[i]/nn[i];
	}

	voxelVolume = hs[0]*hs[1]*hs[2];

	cout<<"VoxelVolume: "<<voxelVolume<<endl;
     
      for(int n = 0; n < dose.size();n++){
	if(roi_float[n]>0){
	  ROIVolume += roi_float[n]*voxelVolume;
	  nVoxelROI++;
	  D[n]=(dose[n]*100);
	}
       }

     cout<<"ROI: "<<ROI<<" Non-zero voxels: "<<nVoxelROI<< " VolumeRoi: "<<ROIVolume<<endl;

	float maxDose = *std::max_element(D.begin(), D.end());
	float binSize = Dgoal/binningDVH;
	int nbin = std::max(binningDVH, static_cast<int>(maxDose / binSize) + 1);

	std::vector<float> bins(nbin, 0.0);
	std::vector<float> hist(nbin, 0.0);


	for (int i = 0; i < nbin; i++) {
	  bins[i] = (i + 0.5) * binSize;
	}
     

	for (size_t i = 0; i < D.size(); ++i) {
	  int imax = static_cast<int>(D[i] / binSize);
	  for (int j = 0; j <= imax; ++j) {
	    hist[j] += ((1*roi_float[i])*voxelVolume);
	  }
	}
	string ROIDVH = dir+ "/" +ROI + fileLabel + ".txt";
	cout<<ROIDVH<<endl;
	ofstream fout (ROIDVH);
	for (int i = 0; i<=nbin -1; i++){
	  fout << bins[i]<<" "<< hist[i]/ROIVolume*100.<< endl;
	}

       
   }

     
 }

 cout<<allROINames<<endl;
 if(plot){
  std::string command = "python3 plotDVH.py -label1 " + fileLabel + " -dir1 " + dir + " -roi " + allROINames + " -Dgoal " + std::to_string(Dgoal);
int result = system(command.c_str());
 }
 
}
      
  
  
  
  

