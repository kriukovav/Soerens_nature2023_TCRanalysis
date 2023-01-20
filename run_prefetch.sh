#!/bin/bash

#SBATCH --cpus-per-task=40         # Run on 40 CPU
#SBATCH --mem=16gb                 # Job memory request
#SBATCH --partition=medium

# read a line with sample names
# $i - a line containing SRR number

    
#body of script. sra-toolkit prefetch.
echo $i
$sra_toolkit_path/prefetch \
-v $i 

echo "${sample_id} finished" > ./log_finished/${i}_finished.txt