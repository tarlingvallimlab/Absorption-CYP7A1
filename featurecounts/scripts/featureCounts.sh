#!/bin/bash
#$ -cwd
#$ -l h_rt=16:00:00        # Wall time (hh:mm:ss)
#$ -l h_data=24G           # Memory per core
#$ -N counts               # Job name
#$ -o logs/counts.out      # Standard output log
#$ -e logs/counts.err      # Standard error log
#$ -m bea                  # Email on (b)egin, (e)nd, and (a)bort
#$ -M rubertg@g.ucla.edu   # Your email address

# Set the directory where your FASTQ files are and change ETV_ID!
BAM_DIR="/u/scratch/r/rubertg/ETV408/fastqs/merged_fastqs/STAR_output/STAR_aligned_BAMS"
cd "$BAM_DIR"
ETV_ID="ETV408"

mkdir counts_output
mkdir counts_output/counts_files

# Initialize the module command (required for non-interactive jobs)
. /u/local/Modules/default/init/modules.sh

# Loop over unique sample IDs
for file in *.bam; do
	# Extract sample name before "Aligned.sortedByCoord.out.bam"
	sample=$(basename "${file}" Aligned.sortedByCoord.out.bam)
	
	~/tools/subread-2.0.6-Linux-x86_64/bin/featureCounts -p -T 8 -s 2 -t exon -g gene_id \
	-a /u/scratch/r/rubertg/Mus_musculus.GRCm38.102.gtf \
	-o ${sample}.counts-gene.txt ${sample}Aligned.sortedByCoord.out.bam
done
	