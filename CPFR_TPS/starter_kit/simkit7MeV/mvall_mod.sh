#!/bin/bash

for run in {1..5}
do
    ./bnn2mhd run_${run}R001_fort.23 dose_tot_run_${run}.mhd -Gy
done
mhd_combine.py -avg dose_tot_run*
