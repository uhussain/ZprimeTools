#!/bin/bash
category=(085)
signal_mx=(10 50 100)
signal_mv=(100 200 1000 1500 1800 2000 2500 3500)

for cat in ${category[@]}; do
    echo $cat Shape
    for mx in ${signal_mx[@]}; do
	for mv in ${signal_mv[@]}; do
	    ./combined.sh zprimeMx${mx}_Mv${mv}_${cat}_shape.txt > Mx${mx}_Mv${mv}_${cat}.txt
	done
    done
done
