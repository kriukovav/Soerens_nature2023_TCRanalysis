#!/bin/bash

# this script executes run_prefetch script in sbatch

sra_toolkit_path=$1 #please specify path to sra-toolkit/bin when you execute this script

mkdir -p ./data/sra/stdout_files_fqd
cd ./data/sra

while read i ; do
    mkdir -p /projects/cdr3_common/user/vkriukova/soerens_nature2023/fastq_${i}
    sbatch -J $i --export=i="$i",sra_toolkit_path="$sra_toolkit_path" --output ./stdout_files_fqd/$i.txt ../../run_fastqdump.sh 
done < ../../SRRxxx.tsv

number_submited=$(cat ../../SRRxxx.tsv | grep -c ^)
number_finished=$(ls ./log_finished | grep -c .txt)

while [ $number_submited -gt $number_finished ]; do
  sleep 10
  number_submited=$(cat ../../SRRxxx.tsv | grep -c ^)
  number_finished=$(ls ./log_finished | grep -c  finished.txt)
done

rm ./log_finished/*