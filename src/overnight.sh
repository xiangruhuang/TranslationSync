#!/bin/bash

stopping=0.0

for sigma in '0.01' '0.04';
do

stopping='0'`echo "0.05+${stopping}" | bc`;
echo ${stopping}

#for g in `seq 1 2`;
#do
#    for i in '2' '6';
#    do
#        make graph${g} noise_ratio=0.${i} sigma=${sigma} stopping=${stopping}
#    done
#done

for g in `seq 3 4`;
do
    for i in '0';
    do
        make graph${g} noise_ratio=0.${i} sigma=${sigma} stopping=${stopping}
    done
done

done

