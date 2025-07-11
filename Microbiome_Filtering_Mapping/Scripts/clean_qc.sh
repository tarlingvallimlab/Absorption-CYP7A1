## portions of this script were specifically written for use in the Hoffman2 High-Performance Compute Cluster at UCLA (https://oarc.ucla.edu/get-help/tools-and-solutions/hoffman2-high-performance-compute-cluster). Please adapt this script to your environment as needed.

set -e 
set -x

mkdir -p QC
mkdir -p QCed

## submit jobs for each sample to run on the Hoffman2 Cluster

SAMPLE=$1

if [[ $SAMPLE == "" ]]; then
  echo "no sample provided. exiting."
  exit;
fi

if [[ $SAMPLE == "gen-jobs" ]]; then
  echo "generating hoffman job scripts"

  cd $SCRATCH/raw_data

  for i in $(cat Metadata.tsv | tail -n+2 | awk '{print $1}'); do
    echo "$SCRATCH/raw_data/clean_qc.sh $i" > $i-clean_qc.sh
    qsub -o joblog.$i-clean_qc.sh -j y -l h_rt=0:45:00,h_data=5G,h_vmem=10G -pe shared 12 -M rochellelai@g.ucla.edu -m bea -cwd $i-clean_qc.sh
  done

exit;
fi

## below lines are executed in the Hoffman2 Cluster job

# prepare job environment

  set +ex
  . /u/local/Modules/default/init/modules.sh
  source /u/home/r/rlai/.bashrc
  conda activate biobakery3.0
  set -ex

  cd $SCRATCH/raw_data

  echo $SAMPLE

# filter and remove host contamination using mus musculus genome GRCm39 (https://genome-idx.s3.amazonaws.com/bt/GRCm39.zip)

  bowtie2 -p 12 -x "/u/scratch/r/rlai/raw_data/GRCm39/GRCm39" \
         -1 ./${SAMPLE}_*/${SAMPLE}_*_L003_R1_001.fastq.gz \
         -2 ./${SAMPLE}_*/${SAMPLE}_*_L003_R2_001.fastq.gz \
	-S ./QC/${SAMPLE}_Mouse_mapped_and_unmapped.sam
  # convert file .sam to .bam
  samtools view -bS ./QC/${SAMPLE}_Mouse_mapped_and_unmapped.sam > ./QC/${SAMPLE}_Mouse_mapped_and_unmapped.bam
  samtools flagstat ./QC/${SAMPLE}_Mouse_mapped_and_unmapped.bam > Mouse_contamination-stat.txt
  # SAMtools SAM-flag filter: get unmapped pairs (both reads R1 and R2 unmapped)
  samtools view -b -f 12 -F 256 ./QC/${SAMPLE}_Mouse_mapped_and_unmapped.bam > ./QC/${SAMPLE}_Mouse_bothReadsUnmapped.bam
  # sort bam file by read name (-n) to have paired reads next to each other
  samtools sort -n -@ 4 ./QC/${SAMPLE}_Mouse_bothReadsUnmapped.bam -o ./QC/${SAMPLE}_Mouse_bothReadsUnmapped_sorted.bam
  samtools fastq -@ 4 ./QC/${SAMPLE}_Mouse_bothReadsUnmapped_sorted.bam \
  -1 ./QC/${SAMPLE}_Mouse_host_removed.R1.fastq.gz \
  -2 ./QC/${SAMPLE}_Mouse_host_removed.R2.fastq.gz \
  -0 /dev/null -s /dev/null -n

# filter and remove index adapter and PhiX (a common Illumina spike-in) sequences

  bbduk.sh -Xmx20g -eoom \
    in1=./QC/${SAMPLE}_Mouse_host_removed.R1.fastq.gz \
    in2=./QC/${SAMPLE}_Mouse_host_removed.R2.fastq.gz \
    out1=./QCed/${SAMPLE}/${SAMPLE}_filtered.R1.fastq.gz \
    out2=./QCed/${SAMPLE}/${SAMPLE}_filtered.R2.fastq.gz \
    ref="adapters,phix" \
    ktrim=r \
    k=23 \
    mink=11 \
    hdist=1 tbo qtrim=rl \
    trimq=25 \
    minlen=100 \
    refstats=./QC/${SAMPLE}_adapter_trimming_stats_per_ref.txt \
    &>> ./QC/${SAMPLE}_bbduk.log.txt