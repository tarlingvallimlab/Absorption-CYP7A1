### VOLCANO PLOT ###
# This script generates a volcano plot based on differentially expressed genes from RNAseq data

# R version 4.2.1 (2022-06-23)
# set working directory and load packages
setwd("/Volumes/tarling-vallim/Alvin/Absorption Paper/Github Repository/RNAseq Enhanced Volcano")
library("EnhancedVolcano") #v1.16.0
library('RColorBrewer') #v1.1.3
library(ggplot2) #v3.5.0

# read differentially expressed genes file
DEgenes <- read.csv("ETV230_DEGS.csv", header = TRUE)

# if you want the plot to highlight specific genes, can uncomment the lines below
# genes <- c("Hmgcs1", "Hmgcr", "Mvk", "Pmvk", "Mvd", "Idi1", "Fdps", "Fdft1", "Sqle", "Lss", "Cyp51a1", "Tm7sf2", "Nsdhl", "Hsd17b7", "Sc5d", "Ebp", "Dhcr7", "Dhcr24", "Srebf2", "Ldlr") # cholesterol biosynthesis
# genes <- c("Fads2", "Elovl2", "Elovl5", "Fads1") # fatty acid elongations

# keyvals <- ifelse(DEgenes$Associated.Gene.Name %in% genes, 'black',
#                   ifelse(DEgenes$log2FoldChange < -1, 'blue',
#                          ifelse(DEgenes$log2FoldChange > 1,'red', 'grey')))

# otherwise, to plot the significantlly up/down volcano plot, use this: 
keyvals <- ifelse(DEgenes$log2FoldChange < -1.5 & DEgenes$padj < 0.05, 'blue',
                         ifelse(DEgenes$log2FoldChange > 1.5 & DEgenes$padj < 0.05,'red', 'grey'))

keyvals[is.na(keyvals)] <- 'grey'
  # names(keyvals)[keyvals == 'black'] <- 'select genes' 
  names(keyvals)[keyvals == 'blue'] <- 'Downregulated genes Log2FC < -1.5' 
  names(keyvals)[keyvals == 'red'] <- 'Upregulated genes Log2FC > 1.5'
  
# counting number of genes changed
sum(DEgenes$log2FoldChange < -1.5 & DEgenes$padj < 0.05) #significantly downregulated genes (11)
sum(DEgenes$log2FoldChange > 1.5 & DEgenes$padj < 0.05) #significantly upregulated genes (18)
  
# variable instructions here: https://github.com/kevinblighe/EnhancedVolcano
  
# create volcano plot 
ETV230DEGVolcano <- EnhancedVolcano::EnhancedVolcano(DEgenes,
                                                    lab=DEgenes$Associated.Gene.Name, selectLab=FALSE,
                                                    labSize = 2,
                                                    xlim = c(-3,3),
                                                    ylim = c(0,25), 
                                                    pointSize = 1,
                                                    labFace = "bold",
                                                    title = 'Cyp7a1 CRISPR vs Control CRISPR',
                                                    x='log2FoldChange',
                                                    y='padj',
                                                    pCutoff=0.05,
                                                    FCcutoff = 1.5,
                                                    colAlpha = 1,
                                                    colCustom = keyvals,
                                                    legendPosition = 'right',
                                                    legendLabSize = 7,
                                                    legendIconSize = 2.0,
                                                    axisLabSize = 8,
                                                    gridlines.minor = FALSE,
                                                    drawConnectors = TRUE,
                                                    widthConnectors = 1,
                                                    boxedLabels = TRUE,
                                                    max.overlaps = Inf
                      
  )
  
ETV230DEGVolcano
ggsave("ETV230Volcano_significant.pdf", height=5, width=7)
