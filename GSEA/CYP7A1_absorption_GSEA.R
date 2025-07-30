#load libraries
if (requireNamespace("clusterProfiler", quietly = TRUE)) {
  print("Package is installed!")
} else {
  BiocManager::install("clusterProfiler")
}
if (requireNamespace("ReactomePA", quietly = TRUE)) {
  print("Package is installed!")
} else {
  BiocManager::install("ReactomePA")
}
if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  print("Package is installed!")
} else {
  BiocManager::install("org.Mm.eg.db")
}
if (requireNamespace("msigdbr", quietly = TRUE)) {
  print("Package is installed!")
} else {
  install.packages('msigdbdf', repos = 'https://igordot.r-universe.dev')
}
library(clusterProfiler)
library(ggplot2)
library(ggnewscale)
library(ReactomePA)
library(msigdbr)
library(tibble)
library(tidyr)
library(dplyr)
library("org.Mm.eg.db")
library(igraph)


#Load in DEGs df
df = read.delim("/Volumes/tarling-vallim/ETV Experimental Data/ETV226-250/ETV230 Cyp CRISPR Western Diet/RNAseq/RNAseq Analysis (GER)/ETV230_Cyp7a1vsCtrl_DEgenes_500.txt", header = TRUE, sep = '\t' )

# Convert Ensembl IDs to Entrez IDs and add Entez ID column to DESeq2 df
df_entrez = bitr(df$Ensembl.Gene.ID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
colnames(df_entrez)[1] <- "Ensembl.Gene.ID"
df <- left_join(df, df_entrez, by = "Ensembl.Gene.ID" )
df <- na.omit(df)

df_sig <- df %>% filter(!(padj > 0.05))

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange 
sig_gene_list <- df_sig$log2FoldChange 


# name the vector
names(original_gene_list) <- df$Associated.Gene.Name
names(sig_gene_list) <- df_sig$Associated.Gene.Name

# omit any NA values 
gene_list<-na.omit(original_gene_list)
sig_gene_list<-na.omit(sig_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
sig_gene_list = sort(sig_gene_list, decreasing = TRUE)

gene_list <- gene_list[!duplicated(names(gene_list))]
sig_gene_list <- sig_gene_list[!duplicated(names(sig_gene_list))]

gene <- names(gene_list)
head(gene)

#GSEA with HALLMARK gene sets
hallmark_sets <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

prep_for_gsea <- function(gene_list) {
  # Remove NA names or NA values
  gene_list <- gene_list[!is.na(names(gene_list)) & !is.na(gene_list)]
  
  # If duplicated gene names, keep the one with the highest absolute ranking
  gene_list <- gene_list[order(-gene_list)]  # Sort decreasing first
  gene_list <- gene_list[!duplicated(names(gene_list))]
  
  # Sort again by decreasing value (GSEA expects this)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}

gsea_result <- GSEA(
  geneList = gene_list,
  TERM2GENE = hallmark_sets,
  pvalueCutoff = 0.05
)
hallmark_gsea_df <- as.data.frame(gsea_result)

p4 <- cnetplot(gsea_result, 
               foldChange=gene_list,
               showCategory = "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
               node_label = "all",
               cex_category = .75,
               cex_gene = 3,
               cex_label_category = 1.5,
               cex_label_gene = .4,
               layout = "circle",
               circular = TRUE, colorEdge = TRUE)
plot(p4)

