#!bin/bash

mapfile -t lines < "./beta.txt"

# 각 줄을 공백으로 구분하여 배열로 변환
IFS=' ' read -r -a line1 <<< "${lines[0]}"

unset IFS

for ((i=0; i<${#line1[@]}; i++)); do
        mkdir b${line1[i]}
        cp -r ./run ./b${line1[i]}
        cp -r ./script ./b${line1[i]}     
        cp -r ./Param ./b${line1[i]}

        sed -i s/beta_/${line1[i]}/g ./b${line1[i]}/run/test.json

        cd ./b${line1[i]}/Param
        sed -i s/betavalue/${line1[i]}/g ./Input_Cal.cpp
        sed -i s/tauvalue/20/g ./Input_Cal.cpp #tau_grid 수정
        cd ../run
        sed -i s/200/20/g test.json #tau_grid 수정
        cd ../..
done

python3 Alp_Gam_output.py