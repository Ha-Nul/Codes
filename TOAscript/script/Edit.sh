#!bin/bash

mapfile -t lines < "../../Param.txt"

# 각 줄을 공백으로 구분하여 배열로 변환
IFS=' ' read -r -a line1 <<< "${lines[0]}"
IFS=' ' read -r -a line2 <<< "${lines[1]}"

unset IFS

for ((i=0; i<${#line1[@]}; i++)); do
    for ((j=0; j<${#line2[@]}; j++)); do
        mkdir ../g_${line1[i]}_a_${line2[j]}_dat

        #mkdir g_${line1[i]}_a_${line2[j]}_dat
        cp -r ../run ../g_${line1[i]}_a_${line2[j]}_dat
        cp ../Param/H_loc_g${line1[i]}.dat ../g_${line1[i]}_a_${line2[j]}_dat/run
        cp ../Param/H_N_g${line1[i]}.dat ../g_${line1[i]}_a_${line2[j]}_dat/run
        cp ../Param/INT_Arr_a${line2[j]}.dat ../g_${line1[i]}_a_${line2[j]}_dat/run
        
        #json파일 수정
        cd ../g_${line1[i]}_a_${line2[j]}_dat/run
        sed -i s/Vfunc.txt/INT_Arr_a${line2[j]}.dat/g test.json
        sed -i s/H_loc/H_loc_g${line1[i]}/g test.json
        sed -i s/H_N/H_N_g${line1[i]}/g test.json
        sed -i s/": 1"/": 2"/g test.json #Nscl 값 수정
        cd ../../script
    done
done

mkdir ../data