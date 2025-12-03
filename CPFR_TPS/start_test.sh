bash FULL_CHAIN_70mm.sh setup.out imgs/Phantom.mhd imgs/PTV_High.mhd imgs/MARKER.mhd imgs/MUSCOLO.mhd imgs/TESSUTI_MOLLI.mhd >full_log.txt 2>&1 &
#bash FULL_CHAIN_70mm.sh setup.out imgs/Phantom.mhd imgs/PTV_High.mhd imgs/MARKER.mhd imgs/MUSCOLO.mhd imgs/TESSUTI_MOLLI.mhd &
disown
tail -F full_log.txt
