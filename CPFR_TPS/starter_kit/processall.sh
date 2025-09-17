#!/bin/bash

for run in {1..50}
do
	 nohup /NFS_homes/software/fluprogfor2020.0.10/flutil/rfluka run_${run}R.inp -e fluka_EF_9MeV.exe -N0 -M1 &
done
