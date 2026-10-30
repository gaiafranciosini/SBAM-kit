#!/bin/bash
# ============================================================
# Usage (from the main folder, where input_log.out is):
#   ./changePulse.sh DOSE_7MeV_1048pulses_CONV_FINAL.mhd <new_pulses>
# ============================================================

die() { echo "Error: $*" >&2; exit 1; }

map=$(realpath "$1" 2>/dev/null)
new="$2"

[ -f "$map" ] && [[ "$new" =~ ^[0-9]+$ ]] \
    || die "usage: $0 DOSE_<E>MeV_<N>pulses_CONV_FINAL.mhd <new_pulses>"

[[ "$(basename "$map")" =~ DOSE_([0-9]+)MeV_([0-9]+)pulses ]] \
    || die "cannot read energy/pulses from the map name."

E=${BASH_REMATCH[1]}
old=${BASH_REMATCH[2]}

[ -f input_log.out ] || die "input_log.out not found (run from the main folder)."
grep -q '^ROI:' input_log.out || die "No 'ROI:' line found in input_log.out."

root=$PWD
kit="$root/starter_kit"
sim_dir=$(dirname "$map")
out="DVH${E}MeV_${new}pulses"

# ROIs (from input_log.out)
read -r _ roi_line < <(grep '^ROI:' input_log.out)
ROIs=($roi_line)
rois_path=(); rois_plan=()

for r in "${ROIs[@]}"; do
    # Evita di duplicare PTV se è già presente in input_log.out
    if [ "$r" != "PTV" ]; then
        rois_path+=("${root}/imgs/${r}_plan.mhd")
        rois_plan+=("${r}_plan")
    fi
done

echo "Original pulses: ${old} | new pulses: ${new} | energy: ${E} MeV"

cd "$sim_dir" || exit 1
mkdir -p "$out"

# Temporary copy of the map (header + data file), removed at exit
tmp="tmp_rescale_$$"
trap 'rm -f "${sim_dir}/${tmp}"*' EXIT

cp "$map" "${tmp}.mhd"

# Estrazione robusta di ElementDataFile (case-insensitive e senza spazi)
raw=$(grep -a -i '^ElementDataFile' "$map" | sed -E 's/^[Ee]lement[Dd]ata[Ff]ile\s*=\s*//i' | tr -d '\r' | xargs)
[ -n "$raw" ] || die "ElementDataFile key not found in $map"

if [ "$raw" != "LOCAL" ]; then
    ext="${raw##*.}"
    data_file="$(dirname "$map")/${raw}"
    [ -f "$data_file" ] || die "Data file ${raw} not found in $(dirname "$map")."
    
    cp "$data_file" "${tmp}.${ext}" || die "Failed to copy data file."
    
    # Compatibile sia con Linux che macOS per sed -i
    if sed --version 2>/dev/null | grep -q GNU; then
        sed -i "s|^[Ee]lement[Dd]ata[Ff]ile.*|ElementDataFile = ${tmp}.${ext}|i" "${tmp}.mhd"
    else
        sed -i '' "s|^[Ee]lement[Dd]ata[Ff]ile.*|ElementDataFile = ${tmp}.${ext}|i" "${tmp}.mhd"
    fi
fi

# Rescale: per pulse, then times the new number of pulses
python3 "$kit/map_rescale.py" "${tmp}.mhd" -divider "$old"    >/dev/null || die "map_rescale (divider) failed."
python3 "$kit/map_rescale.py" "${tmp}.mhd" -multiplier "$new" >/dev/null || die "map_rescale (multiplier) failed."
echo "Map rescaled"

# DVH, plot, PTV summary
"$kit/ComputeDVH/ComputeDVH.py" -Dgoal 5000 -roi "${root}/imgs/PTV_plan.mhd" "${rois_path[@]}" \
    -dose "${tmp}.mhd" -type float -fileLabel "${E}MeV" -dir "$out" >/dev/null || die "ComputeDVH failed."
echo "New DVH curves created"

python3 "$kit/plotDVH.py" -label1 "${E}MeV" -dir1 "$out" -roi PTV_plan "${rois_plan[@]}" >/dev/null || die "plotDVH failed."
python3 "$kit/readDVH.py" -DVH "${out}/PTV_plan${E}MeV.txt" -preD 22.5 > "${out}/DVH_PTV.out" || die "readDVH failed."

[ -f "${E}MeV.png" ] && mv "${E}MeV.png" "${out}/${E}MeV_${new}pulses.png"
echo "DVH plot saved in ${sim_dir}/${out}"
