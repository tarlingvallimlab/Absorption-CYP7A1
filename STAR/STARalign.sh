#!/bin/bash
#$ -cwd
#$ -l h_rt=16:00:00        # Wall time (hh:mm:ss)
#$ -l h_data=48G           # Memory per core
#$ -N STAR_1_35            # Job name
#$ -o logs/STAR_1.out       # Standard output log
#$ -e logs/STAR_1.err       # Standard error log
#$ -m bea                  # Email on (b)egin, (e)nd, and (a)bort
#$ -M rubertg@g.ucla.edu   # Your email address

# Set the directory where your FASTQ files are and change ETV_ID!
FASTQ_DIR="/u/scratch/r/rubertg/ETV408/fastqs/merged_fastqs"
cd "$FASTQ_DIR"
ETV_ID="ETV408"

mkdir STAR_output
mkdir STAR_output/STAR_aligned_BAMS

# Initialize the module command (required for non-interactive jobs)
. /u/local/Modules/default/init/modules.sh

# Load the STAR module from Hoffman cluster
module load star/2.7.10a

# Loop over unique sample IDs
for sample in {3..35}; do
    # Pad with zero if needed (e.g., 01, 02 ... 09)
    sample_padded=$(printf "%d" "$sample")

    # Construct partial sample ID used in file names
    sample_pattern="${sample_padded}_S${sample_padded}"

	# Align paired reads
        STAR \
	--genomeDir /u/scratch/r/rubertg/mm10_dbStar \
        --runThreadN 4 \
        --readFilesIn ${sample_pattern}_R1_merged.fastq.gz ${sample_pattern}_R2_merged.fastq.gz \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMismatchNmax 5 \
        --outFilterMultimapNmax 1 \
        --outFileNamePrefix /u/scratch/r/rubertg/ETV408/fastqs/merged_fastqs/STAR_output/${ETV_ID}_${sample_pattern}
done

mv /u/scratch/r/rubertg/ETV408/fastqs/merged_fastqs/STAR_output/*.bam /u/scratch/r/rubertg/ETV408/fastqs/merged_fastqs/STAR_output/STAR_aligned_BAMS/
