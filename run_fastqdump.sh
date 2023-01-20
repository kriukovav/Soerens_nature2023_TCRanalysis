#!/bin/bash

#SBATCH --cpus-per-task=40         # Run on 40 CPU
#SBATCH --mem=16gb                 # Job memory request
#SBATCH --partition=medium

# $i - a line containing SRR number

    
#body of script. sra-toolkit fasterq-dump.
echo $i
$sra_toolkit_path/fastq-dump \
--outdir /projects/cdr3_common/user/vkriukova/soerens_nature2023/fastq_${i} --split-files ${i}/${i}.sra 

echo "${i} finished" > ./log_finished/${i}_finished.txt