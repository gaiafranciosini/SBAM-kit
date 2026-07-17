#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <sys/stat.h>
#include <unistd.h>
#include "mhd_io.h"

using namespace std;

int main(int argc, char *argv[]) {
    vector<string> filename;
    string dosefile = "";
    string roitype = "";
    string fileLabel = "";
    string dir = ".";
    int npulse = 1;
    float prf = 1.0f; 
    
    int nn[3];
    float x0[3], L[3], hs[3];
    float left[3], up[3], front[3];
    int binningDVH = 400; 
    float maxRangeX = 200.0f; 
    bool plot = false;
    string allROINames = "";

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-plot") {
            plot = true;
        } else if (arg == "-roi") {
            while (i + 1 < argc && argv[i+1][0] != '-') {
                filename.push_back(argv[++i]);
            }
        } else if (arg == "-dose" && i + 1 < argc) {
            dosefile = argv[++i];
        } else if (arg == "-type" && i + 1 < argc) {
            roitype = argv[++i];
        } else if (arg == "-dir" && i + 1 < argc) {
            dir = argv[++i];
        } else if (arg == "-fileLabel" && i + 1 < argc) {
            fileLabel = argv[++i];
        } else if (arg == "-npulse" && i + 1 < argc) {
            npulse = stoi(argv[++i]);
        } else if (arg == "-prf" && i + 1 < argc) {
            prf = stof(argv[++i]);
        }
    }

    if (filename.empty() || dosefile.empty() || roitype.empty() || fileLabel.empty()) {
        cerr << "Error: Missing arguments." << endl;
        return 1;
    }

    float T = (prf > 0) ? (static_cast<float>(npulse) / prf) : 1.0f;
    
    vector<float32> dose;
    read_mhd(dosefile, dose, nn, x0, L, left, up, front);

    for (size_t i = 0; i < dose.size(); ++i) {
        dose[i] /= T;
    }

    // Costruzione nome file: nome_originale + _ + fileLabel
    size_t lastDotDose = dosefile.find_last_of(".");
    string baseName = (lastDotDose == string::npos) ? dosefile : dosefile.substr(0, lastDotDose);
    size_t lastSlashDose = baseName.find_last_of("/");
    string simpleBaseName = (lastSlashDose == string::npos) ? baseName : baseName.substr(lastSlashDose + 1);

    string doseRateOutFile = dir + "/" + simpleBaseName + "_" + fileLabel;
    
    // Configurazione parametri write_mhd come richiesto
    string dtype = "MET_FLOAT"; 
    write_mhd(doseRateOutFile, dtype, (char*)dose.data(), nn, x0, L, left, up, front, true);
    
    // Passaggio a cm per calcolo DVH
    for (int i = 0; i < 3; i++) {
        L[i] /= 10.0;
        x0[i] /= 10.0;
    }

    for (size_t i = 0; i < filename.size(); i++) {
        vector<int16> roi_int;
        vector<float32> roi_float;
        float voxelVolume = 0.0f;
        float ROIVolume = 0.0f;

        size_t start = filename[i].find_last_of("/") + 1;
        size_t end = filename[i].find_last_of(".");
        string ROI_name = filename[i].substr(start, end - start);
        
        allROINames += ROI_name + (i == filename.size() - 1 ? "" : " ");
        vector<double> D_masked(dose.size(), 0.0);

        if (roitype == "int") {
            read_mhd(filename[i], roi_int, nn, x0, L, left, up, front);
            for (int k = 0; k < 3; k++) hs[k] = L[k] / nn[k];
            voxelVolume = hs[0] * hs[1] * hs[2];
            for (size_t n = 0; n < dose.size(); n++) {
                if (roi_int[n] > 0) {
                    ROIVolume += roi_int[n] * voxelVolume;
                    D_masked[n] = dose[n];
                }
            }
        } else {
            read_mhd(filename[i], roi_float, nn, x0, L, left, up, front);
            for (int k = 0; k < 3; k++) hs[k] = L[k] / nn[k];
            voxelVolume = hs[0] * hs[1] * hs[2];
            for (size_t n = 0; n < dose.size(); n++) {
                if (roi_float[n] > 0) {
                    ROIVolume += roi_float[n] * voxelVolume;
                    D_masked[n] = dose[n];
                }
            }
        }

        int nbin = binningDVH; 
        float binSize = maxRangeX / (float)nbin; 

        vector<float> bins(nbin, 0.0);
        vector<float> hist(nbin, 0.0);

        for (int b = 0; b < nbin; b++) {
            bins[b] = (b + 0.5) * binSize;
        }

        for (size_t n = 0; n < D_masked.size(); n++) {
            if (D_masked[n] > 0) {
                float currentDose = static_cast<float>(D_masked[n]);
                float weight = (roitype == "int") ? static_cast<float>(roi_int[n]) : roi_float[n];
                
                int b_idx;
                if (currentDose >= maxRangeX) {
                    b_idx = nbin - 1;
                } else {
                    b_idx = static_cast<int>(currentDose / binSize);
                }

                for (int j = 0; j <= b_idx && j < nbin; j++) {
                    hist[j] += weight * voxelVolume;
                }
            }
        }

        string outPath = dir + "/" + ROI_name + fileLabel + ".txt";
        ofstream fout(outPath);
        for (int b = 0; b < nbin; b++) {
            if (ROIVolume > 0) {
                fout << bins[b] << " " << (hist[b] / ROIVolume) * 100.0 << endl;
            }
        }
        fout.close();
    }

    if (plot) {
        string cmd = "python3 plotDRVH.py -label1 " + fileLabel + " -dir1 " + dir + " -roi " + allROINames;
        system(cmd.c_str());
    }

    return 0;
}
