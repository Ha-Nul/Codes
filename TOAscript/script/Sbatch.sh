#!bin/bash

mapfile -t lines < "../../Param.txt"

# 각 줄을 공백으로 구분하여 배열로 변환
IFS=' ' read -r -a line1 <<< "${lines[0]}"
IFS=' ' read -r -a line2 <<< "${lines[1]}"

unset IFS

for ((i=0; i<${#line1[@]}; i++)); do
    for ((j=0; j<${#line2[@]}; j++)); do
        cd ../g_${line1[i]}_a_${line2[j]}_dat/run
        sbatch json.sh
        cd ../../script
    done
done