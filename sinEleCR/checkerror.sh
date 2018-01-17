#!/bin/bash
array=($(ls -d *.err))
for file in ${array[@]}; do
    if [[ -s ${file} ]]; then
	echo ${file}
	cat ${file}
    fi
done