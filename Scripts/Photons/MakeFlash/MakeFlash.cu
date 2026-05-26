#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include "mhd_io.h"

using namespace std;

// Kernel CUDA
__global__ void computeFlashKernel(int Nvxl, int nFields, float **d_doses, float **d_rates, 
                                   int16_t *ptv, float *doseFLASH, float DT, float FMFmin, 
                                   float Dr_T, int Nfraction) 
{
    int ivxl = blockIdx.x * blockDim.x + threadIdx.x;
    if (ivxl >= Nvxl) return;

    float Dose_flash_per_fraction = 0.0f;
    float Dose_non_flash_per_fraction = 0.0f;

    // Logica PTV: ptv < 1 indica tessuto sano (OAR)
    if (ptv[ivxl] < 6) { 
        for (int i = 0; i < nFields; ++i) {
            float dose_frac = d_doses[i][ivxl]; // Già per frazione
            
            if (d_rates[i][ivxl] > Dr_T) 
                Dose_flash_per_fraction += dose_frac;
            else 
                Dose_non_flash_per_fraction += dose_frac;
        }

        if (Dose_flash_per_fraction > DT) {
            float FMF = (1.0f - FMFmin) * (DT / Dose_flash_per_fraction) + FMFmin;
            Dose_flash_per_fraction *= FMF;
        }
    } 
    else { 
        for (int i = 0; i < nFields; ++i) {
            Dose_non_flash_per_fraction += d_doses[i][ivxl];
        }
    }

    // Risultato finale: (Dose Effettiva) * Numero frazioni
    doseFLASH[ivxl] = (Dose_non_flash_per_fraction + Dose_flash_per_fraction) * (float)Nfraction;
}

int main(int argc, char *argv[]) {
    vector<string> doseFiles;
    vector<float> deliveryTimes;
    string ptvFile, outFile = "output_FLASH.mhd";
    float DT = 5.0f, FMFmin = 0.6f, Dr_T = 10.0f;
    int Nfraction = 5; 

    // Parsing Argomenti (rimane invariato)
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-d" && i + 1 < argc) {
            while (i + 1 < argc && argv[i + 1][0] != '-') doseFiles.push_back(argv[++i]);
        } else if (arg == "-t" && i + 1 < argc) {
            while (i + 1 < argc && argv[i + 1][0] != '-') deliveryTimes.push_back(stof(argv[++i]));
        } else if (arg == "-ptv" && i + 1 < argc) ptvFile = argv[++i];
        else if (arg == "-n" && i + 1 < argc) Nfraction = stoi(argv[++i]);
        else if (arg == "-dr" && i + 1 < argc) Dr_T = stof(argv[++i]);
        else if (arg == "-dt" && i + 1 < argc) DT = stof(argv[++i]);
        else if (arg == "-fmin" && i + 1 < argc) FMFmin = stof(argv[++i]);
        else if (arg == "-out" && i + 1 < argc) outFile = argv[++i];
    }

    if (doseFiles.size() != deliveryTimes.size() || ptvFile.empty() || doseFiles.empty()) {
        cerr << "Errore: parametri mancanti o incoerenti." << endl;
        return 1;
    }

    int nFields = doseFiles.size();
    int nn[3]; float x0[3], L[3], hs[3], left[3], up[3], front[3];

    vector<int16_t> h_ptv;
    read_mhd(ptvFile, h_ptv, nn, x0, L, left, up, front);
    int Nvxl = nn[0] * nn[1] * nn[2];

    int16_t *d_ptv; float *d_out;
    cudaMallocManaged(&d_ptv, Nvxl * sizeof(int16_t));
    cudaMallocManaged(&d_out, Nvxl * sizeof(float));
    memcpy(d_ptv, h_ptv.data(), Nvxl * sizeof(int16_t));

    float **h_dosePtrs = (float**)malloc(nFields * sizeof(float*));
    float **h_ratePtrs = (float**)malloc(nFields * sizeof(float*));

    for (int i = 0; i < nFields; ++i) {
        vector<float> tmpD(Nvxl);
        read_mhd(doseFiles[i], tmpD, nn, x0, L, left, up, front);
        
        cudaMallocManaged(&h_dosePtrs[i], Nvxl * sizeof(float));
        cudaMallocManaged(&h_ratePtrs[i], Nvxl * sizeof(float));

        float time_i = deliveryTimes[i];
        for (int j = 0; j < Nvxl; ++j) {
            h_dosePtrs[i][j] = tmpD[j];
            h_ratePtrs[i][j] = tmpD[j] / time_i;
        }
    }

    // --- CORREZIONE ALLOCAZIONE DOPPI PUNTATORI ---
    float **d_dosePtrs_dev, **d_ratePtrs_dev;
    cudaMalloc(&d_dosePtrs_dev, nFields * sizeof(float*));
    cudaMalloc(&d_ratePtrs_dev, nFields * sizeof(float*)); // Era int*, ora float*

    cudaMemcpy(d_dosePtrs_dev, h_dosePtrs, nFields * sizeof(float*), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ratePtrs_dev, h_ratePtrs, nFields * sizeof(float*), cudaMemcpyHostToDevice);

    int NT = 256;
    int NB = (Nvxl + NT - 1) / NT;
    
    computeFlashKernel<<<NB, NT>>>(Nvxl, nFields, d_dosePtrs_dev, d_ratePtrs_dev, d_ptv, d_out, DT, FMFmin, Dr_T, Nfraction);
    
    cudaError_t err = cudaDeviceSynchronize();
    if(err != cudaSuccess) {
        printf("Errore CUDA: %s\n", cudaGetErrorString(err));
    }

    write_mhd(outFile, d_out, nn, x0, L, left, up, front);

    // Cleanup
    for(int i=0; i<nFields; i++) {
        cudaFree(h_dosePtrs[i]);
        cudaFree(h_ratePtrs[i]);
    }
    free(h_dosePtrs); free(h_ratePtrs);
    cudaFree(d_dosePtrs_dev); cudaFree(d_ratePtrs_dev); 
    cudaFree(d_ptv); cudaFree(d_out);

    cout << "Processo terminato con successo. Output: " << outFile << endl;
    return 0;
}
