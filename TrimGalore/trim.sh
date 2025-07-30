#!/bin/bash

#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y

# Edit below to adjust time and computing power
#$ -l h_rt=3:00:00,h_data=4G
#$ -pe shared 2

# Email address to notify
#$ -M rubertg@g.ucla.edu
# Notify when
#$ -m bea

# load the job environment:
/u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
module load python/3.9.6
export PATH=$PATH:~/miniconda2/bin/cutadapt

# If you want FastQC to run automatically, need to make sure it’s in same directory as TrimGalore
 
for i in *.fastq.gz 
do
/u/home/r/rubertg/tools/TrimGalore/TrimGalore-0.6.6/trim_galore -q 30 --path_to_cutadapt ~/miniconda2/bin/cutadapt --fastqc -o Trimmed/${i}_trimmed.fastq.gz ${i}
done

