### MICROBIOME ANALYSIS ###
## This script takes in the absolute abundance table from MetaPhlAn and runs alpha diversity, beta diversity, and differential taxa analyses ##

# R version 4.2.1 (2022-06-23)
# set working directory and load packages
library(tidyverse) #v2.0.0
library(vegan) #v2.6.4
library(ggplot2) #v3.5.0
library(cowplot) #v1.1.3
library(ggpubr) #v0.6.0
library(Maaslin2) #v1.12.0

options(scipen=999)
setwd("/Volumes/tarling-vallim/Alvin/Absorption Paper/Github Repository/Microbiome_Analysis/Input")

# load input files
input_bug_df = read.table("merged_metaphlan_filtered_absolute_abundance.txt", header = TRUE)
input_bug_df <- input_bug_df[,-2] #remove cladetax_id column
metadata = read_tsv("Metadata.tsv")

### ALPHA DIVERSITY ###

# Transform data for phylum diversity stacked plots to be graphed in Prism Graphpad Software
input_bug_df_phylum <- filter(input_bug_df,!grepl("c__", clade_name)) # only show phylum
input_bug_df_phylum <- filter(input_bug_df_phylum,grepl("p__", clade_name)) # only show phylum
input_bug_df_phylum <- input_bug_df_phylum[,-2] # remove clade_taxid column
input_bug_df_phylum <- rename_at(input_bug_df_phylum, vars(matches("ETV")), ~stringr::str_extract(., "(ETV364)_\\d+")) #clean up column names
input_bug_df_phylum <- as.data.frame(input_bug_df_phylum)
input_bug_df_phylum <- t(input_bug_df_phylum)
colnames(input_bug_df_phylum) <- input_bug_df_phylum[1,]
input_bug_df_phylum <- input_bug_df_phylum[-1,]
input_bug_df_phylum <- as.data.frame(input_bug_df_phylum)
input_bug_df_phylum[,] <- lapply(input_bug_df_phylum[,], as.numeric) # make dataframe numeric

metadata_merge <- metadata %>% remove_rownames %>% column_to_rownames(var="#SampleID")
input_bug_df_phylum_groups <- merge(input_bug_df_phylum, metadata_merge, by = 0) # labels which sample belongs to which groups
input_bug_df_phylum_groups <- input_bug_df_phylum_groups[,-c((ncol(input_bug_df_phylum_groups)-2),(ncol(input_bug_df_phylum_groups)-1))] # remove second to last and third to last columns
rownames(input_bug_df_phylum_groups) <- input_bug_df_phylum_groups[,1]
input_bug_df_phylum_groups <- input_bug_df_phylum_groups[,-1]
write.csv(input_bug_df_phylum_groups, "ETV364_phylum_metaphlan_filtered_absolute_abundance.csv")

# Transform data for alpha and beta diversity by species
input_bug_df <- filter(input_bug_df,grepl("s__", clade_name)) # only show species
input_bug_df <- filter(input_bug_df,!grepl("t__", clade_name)) # only show species
bug_df <- as.data.frame(input_bug_df)
bug_df <- rename_at(bug_df, vars(matches("ETV")), ~stringr::str_extract(., "(ETV364)_\\d+")) #clean up column names
write.csv(bug_df, "ETV364_species_metaphlan_filtered_absolute_abundance.tsv")

bug_df_alpha <- bug_df[,-1] #remove species names
bug_df_alpha[,] <- lapply(bug_df_alpha[,], as.numeric) #make dataframe numeric

# Shannon index
shannon <- diversity(t(bug_df_alpha), index = "shannon")
shannon <- as.data.frame(shannon)
shannon <- rownames_to_column(shannon)
colnames(shannon) <- c("#SampleID", "shannon")

# Simpson index
simpson <- diversity(t(bug_df_alpha), index = "simpson")
simpson <- as.data.frame(simpson)
simpson <- rownames_to_column(simpson)
colnames(simpson) <- c("#SampleID", "simpson")

# Inverted simpson index
invsimpson <- diversity(t(bug_df_alpha), index = "invsimpson")
invsimpson <- as.data.frame(invsimpson)
invsimpson <- rownames_to_column(invsimpson)
colnames(invsimpson) <- c("#SampleID", "invsimpson")

# Observed species (species richness)
observedspec <- specnumber(t(bug_df_alpha))
observedspec <- as.data.frame(observedspec)
observedspec <- rownames_to_column(observedspec)
colnames(observedspec) <- c("#SampleID", "observed_species")

# Merge all diversity indices to metadata file - to be graphed in Prism Graphpad Software
meta_diversity <- merge(metadata, shannon, by = "#SampleID")
meta_diversity <- merge(meta_diversity, simpson, by = "#SampleID")
meta_diversity <- merge(meta_diversity, invsimpson, by = "#SampleID")
meta_diversity <- merge(meta_diversity, observedspec, by = "#SampleID")
write.csv(meta_diversity, "ETV364_filtered_alpha_diversity.csv")

### BETA DIVERSITY ###

# Format species abundance data
bug_mat = bug_df |> 
  column_to_rownames("clade_name") |> 
  as.matrix() |>
  t()

mode(bug_mat) = "numeric" #change from string to numeric values

dist_mat = vegdist(bug_mat, method = "bray")

# Bray-Curtis dissimilarity
cmd_res = cmdscale(dist_mat, 
                   k = (nrow(bug_mat) - 1),
                   eig = TRUE)

str(cmd_res)

pcoa_df = tibble(PC1 = cmd_res$points[,1], 
                 PC2 = cmd_res$points[,2])

# PERMANOVA test
set.seed(10)
permanovastats <- adonis2(dist_mat ~ Group, data = metadata, permutations = 1000)
write.csv(permanovastats, "ETV364_filtered_permanovatest_stats_Cyp7a1.csv")

# Graphing the Principal Coordinates Analysis plot
p = ggplot(pcoa_df, aes(x = PC1, y = PC2)) + 
  geom_point()

pcoa_meta = bind_cols(pcoa_df, metadata)

p_diag = ggplot(pcoa_meta,
                aes(x = PC1, y = PC2, color = Group)) + 
  geom_point() #+ # show dots
  # geom_text(label = pcoa_meta$MouseNumber, nudge_x = 0.02, check_overlap = TRUE) # to show Mouse ID

p = p_diag + theme_light()
p = p + scale_color_manual(values = c("Control CRISPR" = "#808080","Cyp7a1 CRISPR" = "#0080FF"))
p_final = p + labs(color = "Group", 
                   title = "Principal coordinates by group")

p_final
ggsave(filename = "ETV364_pcoa_plot.png", 
       plot = p_final,
       width = 5,
       height = 4)

### DIFFERENTIAL TAXA ###

diff_taxa_df <- as.data.frame(bug_df)
diff_taxa_df <- t(diff_taxa_df)
colnames(diff_taxa_df) <- diff_taxa_df[1,]
diff_taxa_df <- diff_taxa_df[-1,]
diff_taxa_df <- as.data.frame(diff_taxa_df)
diff_taxa_df[,] <- lapply(diff_taxa_df[,], as.numeric) #m ake dataframe numeric

# Add groups column to end of dataframe to sort groups
metadata_merge <- metadata %>% remove_rownames %>% column_to_rownames(var="#SampleID")
diff_taxa_df_groups <- merge(diff_taxa_df, metadata_merge, by = 0)
diff_taxa_df_groups <- diff_taxa_df_groups[,-c((ncol(diff_taxa_df_groups)-2),(ncol(diff_taxa_df_groups)-1))] # remove second to last and third to last columns
rownames(diff_taxa_df_groups) <- diff_taxa_df_groups[,1]
diff_taxa_df_groups <- diff_taxa_df_groups[,-1]

# Control vs Cyp7a1  CRISPR differential species
diff_taxa_df_groups_Cyp7a1 <- filter(diff_taxa_df_groups,grepl("Control CRISPR|Cyp7a1 CRISPR", Group))
diff_taxa_df_groups_Cyp7a1 <- diff_taxa_df_groups_Cyp7a1[,-ncol(diff_taxa_df_groups_Cyp7a1)] # remove groups column

diff_taxa_metadata_Cyp7a1 <- filter(metadata_merge, grepl("Control CRISPR|Cyp7a1 CRISPR", Group))

fit_data = Maaslin2(input_data=diff_taxa_df_groups_Cyp7a1, input_metadata=diff_taxa_metadata_Cyp7a1, 
                    output = "ETV364_differential_taxa_Maaslin2_Cyp7a1", 
                    fixed_effects = c("Group"),
                    reference = c("Group,Control CRISPR"),
                    normalization = "TSS", transform = "log", plot_heatmap = FALSE, plot_scatter = FALSE)
