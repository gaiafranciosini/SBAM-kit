#!/usr/bin/sh

file=$1
filename="${file%%.*}"
filenameclamped="${filename}_clamped.mhd"

python3 mhd_clamp.py ${file} -o ${filenameclamped} -vmin -1000 -vmax 1376
python3 fredCTHU2flukaVoxels.py ${filenameclamped} hu2materials_fredEMGPU.txt
mv CTHU.vxl ${filename}.vxl
rm ${filenameclamped}

echo "${filename}.vxl"
