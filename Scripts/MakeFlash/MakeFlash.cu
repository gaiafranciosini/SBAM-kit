#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include "mhd_io.h"

using namespace std;

// Kernel per gestire N campi e parametri variabili
__global__ void computeFlashKernel(int Nvxl, int nFields, float32 **d_doses, float32 **d_rates, 
                                   int16 *ptv, float32 *doseFLASH, float DT, float FMFmin, float Dr_T, int Nfraction) 
{
    int ivxl = blockIdx.x * blockDim.x + threadIdx.x;
    if (ivxl >= Nvxl) return;

    float Dose_flash_per_fraction = 0.0f;
    float Dose_non_flash_per_fraction = 0.0f;

    if (ptv[ivxl] < 1) { 
        for (int i = 0; i < nFields; ++i) {
            if (d_rates[i][ivxl] > Dr_T) 
                Dose_flash_per_fraction += d_doses[i][ivxl];
            else 
                Dose_non_flash_per_fraction += d_doses[i][ivxl];
        }

        // Calcolo FMF basato sulla dose FLASH della frazione
        if (Dose_flash_per_fraction > DT) {
            float FMF = (1.0f - FMFmin) * (DT / Dose_flash_per_fraction) + FMFmin;
            Dose_flash_per_fraction *= FMF;
        }
    } else { // Dentro il PTV 
        for (int i = 0; i < nFields; ++i) {
            Dose_non_flash_per_fraction += d_doses[i][ivxl];
        }
    }

    // Risultato finale scalato sulle frazioni totali
    doseFLASH[ivxl] = (Dose_non_flash_per_fraction + Dose_flash_per_fraction) * Nfraction;
}

int main(int argc, char *argv[]) {
    vector<string> doseFiles, rateFiles;
    string ptvFile, outFile = "output_FLASH.mhd";
    float DT = 6.0f, FMFmin = 0.6f, Dr_T = 10.0f;
    int Nfraction = 1;

    // Parser degli argomenti da linea di comando
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-d" && i + 1 < argc) {
            while (i + 1 < argc && argv[i + 1][0] != '-') doseFiles.push_back(argv[++i]);
        } else if (arg == "-r" && i + 1 < argc) {
            while (i + 1 < argc && argv[i + 1][0] != '-') rateFiles.push_back(argv[++i]);
        } else if (arg == "-ptv" && i + 1 < argc) ptvFile = argv[++i];
        else if (arg == "-n" && i + 1 < argc) Nfraction = stoi(argv[++i]);
        else if (arg == "-dr" && i + 1 < argc) Dr_T = stof(argv[++i]);
        else if (arg == "-dt" && i + 1 < argc) DT = stof(argv[++i]);      // Nuova soglia dose FMF
        else if (arg == "-fmin" && i + 1 < argc) FMFmin = stof(argv[++i]); // Nuova FMF min
        else if (arg == "-out" && i + 1 < argc) outFile = argv[++i];
    }

    if (doseFiles.size() != rateFiles.size() || ptvFile.empty() || doseFiles.empty()) {
        cerr << "Utilizzo: ./ComputeFlash.x -d [dosi...] -r [rate...] -ptv [file] -n [frazioni] -dr [soglia_rate] -dt [soglia_dose] -fmin [fmf_min]" << endl;
        return 1;
    }

    int nFields = doseFiles.size();
    int nn[3]; float x0[3], L[3], hs[3], left[3], up[3], front[3];

    // Caricamento PTV per geometria
    vector<int16> h_ptv;
    read_mhd(ptvFile, h_ptv, nn, x0, L, left, up, front);
    int Nvxl = nn[0] * nn[1] * nn[2];

    // Allocazione Output e PTV in memoria gestita
    int16 *d_ptv; float32 *d_out;
    cudaMallocManaged(&d_ptv, Nvxl * sizeof(int16));
    cudaMallocManaged(&d_out, Nvxl * sizeof(float32));
    memcpy(d_ptv, h_ptv.data(), Nvxl * sizeof(int16));

    // Liste di puntatori per gestire campi multipli
    float32 **h_dosePtrs = (float32**)malloc(nFields * sizeof(float32*));
    float32 **h_ratePtrs = (float32**)malloc(nFields * sizeof(float32*));

    cout << "Caricamento di " << nFields << " campi..." << endl;
    for (int i = 0; i < nFields; ++i) {
        vector<float32> tmpD, tmpR;
        read_mhd(doseFiles[i], tmpD, nn, x0, L, left, up, front);
        read_mhd(rateFiles[i], tmpR, nn, x0, L, left, up, front);
        
        cudaMallocManaged(&h_dosePtrs[i], Nvxl * sizeof(float32));
        cudaMallocManaged(&h_ratePtrs[i], Nvxl * sizeof(float32));
        
        memcpy(h_dosePtrs[i], tmpD.data(), Nvxl * sizeof(float32));
        memcpy(h_ratePtrs[i], tmpR.data(), Nvxl * sizeof(float32));
    }

    // Copia i riferimenti dei puntatori sulla GPU
    float32 **d_dosePtrs, **d_ratePtrs;
    cudaMalloc(&d_dosePtrs, nFields * sizeof(float32*));
    cudaMalloc(&d_ratePtrs, nFields * sizeof(float32*));
    cudaMemcpy(d_dosePtrs, h_dosePtrs, nFields * sizeof(float32*), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ratePtrs, h_ratePtrs, nFields * sizeof(float32*), cudaMemcpyHostToDevice);

    // Esecuzione Kernel
    int NT = 256;
    int NB = (Nvxl + NT - 1) / NT;
    cout << "Esecuzione calcolo FLASH (Dr_T=" << Dr_T << ", DT=" << DT << ", FMFmin=" << FMFmin << ")..." << endl;
    computeFlashKernel<<<NB, NT>>>(Nvxl, nFields, d_dosePtrs, d_ratePtrs, d_ptv, d_out, DT, FMFmin, Dr_T, Nfraction);
    cudaDeviceSynchronize();

    // Salvataggio
    vector<float32> h_out(Nvxl);
    memcpy(h_out.data(), d_out, Nvxl * sizeof(float32));
    //    write_mhd(outFile, h_out, nn, x0, L, left, up, front);
    write_mhd(outFile, h_out.data(), nn, x0, L, left, up, front);
    // Free
    for(int i=0; i<nFields; i++) { cudaFree(h_dosePtrs[i]); cudaFree(h_ratePtrs[i]); }
    free(h_dosePtrs); free(h_ratePtrs);
    cudaFree(d_dosePtrs); cudaFree(d_ratePtrs); cudaFree(d_ptv); cudaFree(d_out);

    cout << "Fine. File salvato in: " << outFile << endl;
    return 0;
}
