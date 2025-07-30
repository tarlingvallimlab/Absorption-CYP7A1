#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####################################################
#~~ ETV230 RNA-seq Analysis (GER August 2024)  ~##
####################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library("EnhancedVolcano")
library("pheatmap")
library("RColorBrewer")
library("tidyr")
library('apeglm')
library("stringr")
library(tidyverse)
library("ggplot2")
library("DESeq2")
library(data.table)
library(dplyr)
library(matrixStats)
library(RColorBrewer)

~~~#######################################~~~
  ##           1. Load in Data              ~##
  ~~~#######################################~~~
  
  #1.1 Load in raw count files for all smaples
  directory <- "/Volumes/tarling-vallim/ETV Experimental Data/ETV226-250/ETV230 Cyp CRISPR Western Diet/RNAseq/RNAseq Analysis (GER)/ETV230 counts"
  samplefiles <- (list.files(directory))
  
  #1.2 load in the gene key for the Ensembl IDs
  key = read.table('/Users/gabriellarubert/Desktop/Tarling-Vallim Lab/RNASeq/mart_export_ensembl84.txt', header = TRUE, sep = '\t')  
  
  #1.3 Load the sample codes into 
  read.table("/Volumes/tarling-vallim/ETV Experimental Data/ETV226-250/ETV230 Cyp CRISPR Western Diet/RNAseq/RNAseq Analysis (GER)/ETV230_countskey.txt", sep = '\t', header = TRUE)
  ETV230metadata <- read.table("/Volumes/tarling-vallim/ETV Experimental Data/ETV226-250/ETV230 Cyp CRISPR Western Diet/RNAseq/RNAseq Analysis (GER)/ETV230_countskey.txt", sep = '\t', header = TRUE)
  
  #~~~#######################################~~~#
  ##    2. Get norm counts for all samples    ~##
  #~~~#######################################~~~#
  #2.1 Create sample table
  ETV230sampleTable <- data.frame(sampleName = ETV230metadata$filename,
                                  fileName = ETV230metadata$filename,
                                  condition= ETV020metadata$treatment)
  
  ETV230ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = ETV230sampleTable,
                                               directory = directory,
                                               design= ~condition)
  
  #2.2 Filter out reads whose row sums are less than 10 or 100 or 500
  ETV230ddsHTSeq <- ETV230ddsHTSeq[rowSums(counts(ETV230ddsHTSeq))>500,]
  
  #2.3 Run DESeq
  ETV230dds <- DESeq(ETV230ddsHTSeq)
  normcount <-as.data.frame(counts(ETV230dds, normalized=TRUE))
  rawcount <-as.data.frame(counts(ETV230dds, normalized=FALSE))
  
  #Converts the row name to a column named "Gene Name"
  normcount<- rownames_to_column(normcount, var="Ensembl.Gene.ID") %>% left_join(., key, by = "Ensembl.Gene.ID")
  rawcount<- rownames_to_column(rawcount, var="Ensembl.Gene.ID") %>% left_join(., key, by = "Ensembl.Gene.ID")
  
  #Save to computer 
  write.table(normcount, file='/Volumes/tarling-vallim/ETV Experimental Data/ETV226-250/ETV230 Cyp CRISPR Western Diet/RNAseq/RNAseq Analysis (GER)/ETV230_normalized_counts_250424 (GER).txt', sep = '\t', row.names = FALSE)
  write.table(rawcount, file='/Volumes/tarling-vallim/ETV Experimental Data/ETV226-250/ETV230 Cyp CRISPR Western Diet/RNAseq/RNAseq Analysis (GER)/ETV230_raw_counts_250424 (GER).txt', sep = '\t', row.names = FALSE)
