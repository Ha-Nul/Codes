#!bin/bash

mapfile -t lines < "../../Param.txt"

# 각 줄을 공백으로 구분하여 배열로 변환
IFS=' ' read -r -a line1 <<< "${lines[0]}"
IFS=' ' read -r -a line2 <<< "${lines[1]}"

unset IFS

for ((i=0; i<${#line1[@]}; i++)); do
    for ((j=0; j<${#line2[@]}; j++)); do        
        #json파일 수정
        cd ../g_${line1[i]}_a_${line2[j]}_dat/run
        mv Gt.dat Gt_g${line1[i]}_a${line2[j]}.dat
        mv StS0.dat StS0_g${line1[i]}_a${line2[j]}.dat
        mv St.dat St_g${line1[i]}_a${line2[j]}.dat

        cp Gt_g${line1[i]}_a${line2[j]}.dat ../../data
        cp StS0_g${line1[i]}_a${line2[j]}.dat ../../data
        cp St_g${line1[i]}_a${line2[j]}.dat ../../data

        cd ../../script
    done
done