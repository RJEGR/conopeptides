#!/bin/bash

# Create the sample data which contains four cols
# Factor 1	Factor 2	Forward	Reverse

for file in $(ls *_R1.fq.gz | grep fq)
do
basename=${file##*/}
base="${basename%*_R1.fq.gz}"
bs="${base%*_E_*}"
fwr=${base}_R1.fq.gz
rev=${base}_R2.fq.gz
echo "$bs" "$bs" `printf "$PWD/$fwr"` `printf "$PWD/$rev"`
done > samples.txt
