#!/opt/homebrew/bin/bash
mapfile -t lines < "Param.txt"

# 각 줄을 공백으로 구분하여 배열로 변환
IFS=' ' read -r -a line1 <<< "${lines[0]}"
IFS=' ' read -r -a line2 <<< "${lines[1]}"

unset IFS

for ((i=0; i<${#line1[@]}; i++)); do
    for ((j=0; j<${#line2[@]}; j++)); do
        echo "Passing values: ${line1[i]} ${line2[j]}"
        echo "${line1[i]} ${line2[j]}" | ./HYBtest
    done
done