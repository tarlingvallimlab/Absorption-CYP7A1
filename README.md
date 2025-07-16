# Absorption-CYP7A1
Code to analyze bulk liver RNA sequencing and microbiome metagenomic shotgun sequencing with loss of CYP7A1.

## Bulk liver RNA sequencing computational methods

Objective: Conduct gene expression computational analysis for 16 livers obtained from control CRISPR and *Cyp7a1* CRISPR mice. The samples are labeled as follows:

| Sample ID | Group | 
|---|---|
| 7 | Control CRISPR | 
| 12 | Control CRISPR |
| 13 | Control CRISPR |
| 16 | Control CRISPR |
| 20 | Control CRISPR |
| 25 | Control CRISPR |
| 35 | Control CRISPR |
| 38 | Control CRISPR |
| 2 | *Cyp7a1* CRISPR |
| 5 | *Cyp7a1* CRISPR |
| 17 | *Cyp7a1* CRISPR |
| 19 | *Cyp7a1* CRISPR |
| 24 | *Cyp7a1* CRISPR |
| 26 | *Cyp7a1* CRISPR |
| 28 | *Cyp7a1* CRISPR |
| 36 | *Cyp7a1* CRISPR |

### Pipeline
1.	Trim adapters and low-quality bases from raw FASTQ files using Trim Galore.
2.	Map trimmed FASTQ files to Mus musculus genome (from Ensembl) using STAR.
3.	Annotate aligned reads and quantify gene counts using Subread featureCounts.
4.	Conduct differential gene expression analysis using DESeq2 package in R.
5.	Conduct Gene Set Enrichment Analysis (GSEA) using clusterProfiler package in R.

### > Trim Galore
This folder contains input and output for trimgalore.sh

#### 1. trimgalore.sh
This script trims adapter sequences and low-quality bases CYP7A1 FASTQ files (not in folder).
+ Trim Galore v0.6.10
+ Cutadapt v4.0

*Citation:*
+ Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnetJournal*. 2011;17:10-2. https://doi.org/10.14806/ej.17.1.200

### > STAR
This folder contains input and outout for STARalign.sh

#### 1. STARalign.sh
This script aligns Trim Galore-trimmed FASTQ files to the *Mus musculus* genome using Terminal (Mac). 
+ STAR v2.7.11a
+ + *Mus musculus* mm10 reference genome: GCF_000001635.26_GRCm38.p6_genomic.fna

*Citation:*
+ Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, et al. STAR: ultrafast universal RNA-seq aligner. *Bioinforma. Oxf. Engl.* 2013;29:15–21. https://doi.org/10.1093/bioinformatics/bts635. PMID: 23104886

### > featureCounts
This folder contains input and output for featurecounts.sh

#### 1. featurecounts.sh
This script takes in the STAR-mapped reads (INSERT HERE.bam) and the Mus Musculus mm10 reference genome annotation file including chromosomal coordinates and outputs a large matrix of the number of reads assigned to unique gene features. 

+ Subread v3.6.3
+ *Mus musculus* mm10 annotation file: Mus_musculus.GRCm38.102.gtf

*Citation:*
+ Liao, Y., Smyth, G. K. & Shi, W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinforma. Oxf. Engl.* 2014;30:923–930. https://doi.org/10.1093/bioinformatics/btt656. PMID: 24227677

### > DESeq2
This folder contains the input and output for CYP7A1_absorption_DESeq2.R

#### 1. CYP7A1_absorption_DESeq2.R
This script takes the raw count matrix created with featureCounts and processes it through DESeq2 modelling to find differentially expressed genes in livers between control CRISPR and *Cyp7a1* CRISPR mice (FDR 5%). This script also generates and exports a dataframe of normalized counts for all samples (DESeq2's median of ratios scaling) from the raw input. Genes whose counts sum to less than 500 across all samples are excluded from normalized counts matrix. 

+ R Version
+ dplyr v
+ org.Mm.ed.db

*Citation:*
+ Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. 2014;15(12):550. https://doi.org/10.1186/s13059-014-0550-8. PMID: 25516281

### > GSEA
This folder contains the input and output for CYP7A1_absorption_GSEA.R

#### 1. CYP7A1_absorption_GSEA.R
This script performs GSEA at FDR 10% using the clusterProfiler package 

+ R Version 
+ clusterProfiler 3.10.1
+ dplyr
+ org.Mm.eg.db 3.7.0

*Citation:*
+ Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. 2014;15(12):550. https://doi.org/10.1186/s13059-014-0550-8. PMID: 25516281
  
## Microbiome metagenomic shotgun sequencing computational methods

Objective: Conduct microbiome computational analysis for 19 metagenomic shotgun sequencing samples obtained from cecal contents of control CRISPR and *Cyp7a1* CRISPR mice. The samples are labeled as follows:

| Sample ID | Group | 
|---|---|
| ETV364_1 | Control CRISPR | 
| ETV364_3 | Control CRISPR |
| ETV364_8 | Control CRISPR |
| ETV364_15 | Control CRISPR |
| ETV364_17 | Control CRISPR |
| ETV364_21 | Control CRISPR |
| ETV364_29 | Control CRISPR |
| ETV364_33 | Control CRISPR |
| ETV364_35 |Control CRISPR |
| ETV364_36 | Control CRISPR |
| ETV364_7 | *Cyp7a1* CRISPR |
| ETV364_12 | *Cyp7a1* CRISPR |
| ETV364_16 | *Cyp7a1* CRISPR |
| ETV364_20 | *Cyp7a1* CRISPR |
| ETV364_30 | *Cyp7a1* CRISPR |
| ETV364_32 | *Cyp7a1* CRISPR |
| ETV364_37 | *Cyp7a1* CRISPR |
| ETV364_38 | *Cyp7a1* CRISPR |
| ETV364_50 | *Cyp7a1* CRISPR |

### Pipeline
1.	Filter out host contamination from FASTQ files by removing sequences that map to the Mus musculus genome using Bowtie2.
2.	Filter out index adapters and PhiX sequences from FASTQ files using BBDuk.
3.	Perform metagenomic taxonomic profiling on filtered FASTQ files using MetaPhlAn.
4.	Conduct alpha diversity analysis from MetaPhlAn-mapped abundance counts using vegan package in R.
5.	Conduct beta diversity analysis from MetaPhlAn-mapped abundance counts using vegan and tidyverse packages in R.
6.	Conduct differential taxa analysis MetaPhlAn-mapped abundance counts using MaAsLin2 in R.

### > Microbiome_Filtering_Mapping
This folder contains three subfolders, “Scripts,” “Input,” and “Output.”

#### 1. clean_qc.sh
This script takes in the sample manifest (Input > Metadata.tsv) and raw metagenomic shotgun sequencing FASTQ files (not in folder) to remove host contamination, index adapters, and PhiX (a common Illumina spike-in) sequences using Hoffman2 High-Performance Compute Cluster.

*Filter out host contamination*
+ Bowtie2 v2.5.4
+ BioBakery v3.0
+ *Mus musculus* genome reference: GRCm39 (https://genome-idx.s3.amazonaws.com/bt/GRCm39.zip)

*Filter out index adapters and PhiX sequences*
+ BBDuk, which is part of the BBTools suite v37.62

*Citation:*
+ Langmead B, Wilks C, Antonescu V, Charles R. Scaling read aligners to hundreds of threads on general-purpose processors. *Bioinformatics*. 2019;35(3):421-432. doi:10.1093/bioinformatics/bty648. PMID: 30020410
+ Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. *Nat Methods*. 2012;9(4):357-359. doi:10.1038/nmeth.1923. PMID: 22388286
+ Bushnell B, Rood J, Singer E. BBMerge - Accurate paired shotgun read merging via overlap. *PLoS One*. 2017;12(10):e0185056. doi:10.1371/journal.pone.0185056. PMID: 29073143.

#### 2. metaphlan.sh
This script takes in the sample manifest (Input > Metadata.tsv) and filtered FASTQ files (not in folder) produced by the clean_qc.sh script to perform metagenomic taxonomic profiling using Hoffman2 High-Performance Compute Cluster.

+ MetaPhlAn v4.1.0
+ BioBakery v3.0
+ Microbial genome reference: mpa_vOct22_CHOCOPhlAnSGB_20221

*Citation:*
+ Blanco-Míguez A, Beghini F, Cumbo F, et al. Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. *Nat Biotechnol*. 2023;41(11):1633-1644. doi:10.1038/s41587-023-01688-w. PMID: 36823356

#### 3. merge_metaphlan_tables_abs.py
This script parses the abundance tables produced by metaphlan.sh (.txt files, not in folder) into one absolute abundance counts file (Output > merged_metaphlan_filtered_absolute_abundance.txt) using local machine. Adapted from https://github.com/timyerg/Metaphlan-absolute-abundance-merger.

+ Python v3.11.7

#### 4. merge_metaphlan_tables_rel.py
This script parses the abundance tables produced by metaphlan.sh (.txt files, not in folder) into one relative abundance counts file (Output > merged_metaphlan_filtered_absolute_abundance.txt) using local machine. Adapted from https://github.com/timyerg/Metaphlan-absolute-abundance-merger.

+ Python v3.11.7

### > Microbiome_Analysis
This folder contains two subfolders, “Input” and “Output,” which are related to ETV364_Microbiome_Analysis.R.

#### 1. ETV364_Microbiome_Analysis.R
This script takes in the sample manifest (Input > Metadata.tsv) and MetaPhlAn-mapped absolute abundance file (Input > merged_metaphlan_filtered_absolute_abundance.txt) to conduct alpha diversity, beta diversity, and differential taxa analysis using R.

+ R v4.2.1 (2022-06-23)

*Alpha diversity*

Portion of the script creates a filtered MetaPhlAn-mapped absolute abundance file that lists taxa by phylum name (Output > ETV364_phylum_metaphlan_filtered_absolute_abundance.txt) and species name (Output > ETV364_species_metaphlan_filtered_absolute_abundance.txt) for better organization. Additionally, Shannon’s index, Simpson’s index, inverted Simpson’s index, and observed number are calculated at the species level and summarized in one file (Output > ETV364_filtered_alpha_diversity.csv).

+ vegan v2.6.4

*Beta diversity*

Portion of the script performs principal coordinate analysis at the species level with Bray-Curtis dissimilarity and conducts a PERMANOVA test (Output > ETV364_filtered_permanovatest_stats_Cyp7a1.csv). A plot will also be generated using principal coordinate analysis results (Output > ETV364_pcoa_plot.png).

+ tidyverse v2.0.0
+ vegan v2.6.4
+ ggplot2 v3.5.0
+ cowplot v1.1.3
+ ggpubr v0.6.0

*Differential taxa*

Portion of the script determines significant differential taxa at the species level through MaAsLin2 regression modeling (Output > ETV364_differential_taxa_Maaslin2_Cyp7a1).

+ MaAsLin2 v1.12.0

*Citation*
+ Mallick H, Rahnavard A, McIver LJ, et al. Multivariable association discovery in population-scale meta-omics studies. Coelho LP, ed. *PLoS Comput Biol*. 2021;17(11):e1009442. doi:10.1371/journal.pcbi.1009442. PMID: 34784344
