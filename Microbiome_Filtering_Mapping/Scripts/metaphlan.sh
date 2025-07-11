## portions of this script were specifically written for use in the Hoffman2 High-Performance Compute Cluster at UCLA (https://oarc.ucla.edu/get-help/tools-and-solutions/hoffman2-high-performance-compute-cluster). Please adapt this script to your environment as needed.

set -e 
set -x

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
    echo "$SCRATCH/raw_data/metaphlan.sh $i" > $i-metaphlan.sh
    qsub -o joblog.$i-metaphlan.sh -j y -l h_rt=0:60:00,h_data=5G,h_vmem=50G -pe shared 12 -M rochellelai@g.ucla.edu -m bea -cwd $i-metaphlan.sh
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

# make sure there is a clean environment with no redundant folders

rm -f $SCRATCH/raw_data/QCed/${SAMPLE}_filtered.fastq_metagenome.bowtie2.bz2
rm -f $SCRATCH/raw_data/QCed/${SAMPLE}_filtered.fastq_profiled_metagenome.txt

# conduct metagenomic taxonomic profiling using mpa_vOct22_CHOCOPhlAnSGB_202212 database

metaphlan $SCRATCH/raw_data/QC/${SAMPLE}_Mouse_host_removed.R1.fastq.gz,$SCRATCH/raw_data/QC/${SAMPLE}_Mouse_host_removed.R2.fastq.gz --bowtie2out $SCRATCH/raw_data/QCed/${SAMPLE}_filtered.fastq_metagenome.bowtie2.bz2 --nproc 12 --input_type fastq -o $SCRATCH/raw_data/QCed/${SAMPLE}_filtered.fastq_profiled_metagenome.txt --bowtie2db $SCRATCH/metaphlan_databases --index mpa_vOct22_CHOCOPhlAnSGB_202212 -t rel_ab_w_read_stats
