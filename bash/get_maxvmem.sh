#!/bin/bash

echo "Please enter job numbers separated by space:"
read -a jobnums

for jobnum in "${jobnums[@]}"; do
    # Get job information with qacct
    jobinfo=$(qacct -j $jobnum)

    # Get maxvmem and sample name from job information
    maxvmem=$(echo "$jobinfo" | grep -w maxvmem | awk '{print $2}')
    samplename=$(echo "$jobinfo" | grep -w jobname | awk '{print $2}')

    echo -e "$jobnum:\t$samplename\t$maxvmem"

done


