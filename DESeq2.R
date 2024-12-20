
### This script was written by Marion Hardy for a collaboration with 
## Akshaya Karthikeyan from the Lombard lab

## This sets up your DESeq2 object for further data analysis
## This is the C42B cell line treated with different mono/bi-therapies

library(tidyverse)
library(readxl)
library(DESeq2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)

## Load the raw count matrices, two in this case

counts = read.table("./data/rawcountdata_C42B.txt", header = T)

## Load metadata (sample annotation file)

coldata = read_xlsx("./data/metadata.xlsx", sheet = "Recoded")
rownames(coldata) = coldata$ID

## Clean up counts matrix and ensembl id

strrep = sub(pattern = "\\.(.*)","",counts$X)
counts$X = strrep
rownames(counts) = counts$X
counts = counts %>% select(!c(X,gene_name))

# Create the full model for comparison of samples
# AK said compare DMSO to Ola day 5 and Ola day 9

dds = DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~Treatment) 
# checked experimental design, looks good, no recorded confounding variable
# no need for pairwise comparisons either

# Generate a linear model
dds$Treatment = relevel(dds$Treatment, "DMSO")
dds = DESeq(dds)
resultsNames(dds)
# Generate models for the comparisons between mono and bi-therapies
dds$Treatment = relevel(dds$Treatment, "Abi_Ola")
dds1 = DESeq(dds)
resultsNames(dds1)

dds$Treatment = relevel(dds$Treatment, "Enz_Ola")
dds2 = DESeq(dds)
resultsNames(dds2)

## Checking distribution of counts per sample

as_tibble(assay(dds)) %>%
  gather(sample, value = counts) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ sample) # looks ok

# Checking size factors and dispersion

sizeFactors(dds) # only takes into account the sequencing depth, looks ok
plotDispEsts(dds) # verifies normalization, graph looks a-ok

# Checking PCA

rld = vst(dds)

plotPCA(rld,intgroup="Treatment") + 
  theme_bw()+
  labs(title = 'PCA per treatment')

# Checking sample similarity

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$Treatment, sep="-")
colnames(sampleDistMatrix) <- paste(rld$Treatment, sep="-")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# you can see one of the day9 is clustering with day5
# but the pca shows that very little of the variance is explained
# by the length of Ola treatment so it makes sense that they are similar
# enough to cluster together when computing euclidian distances

## Saving the DESeq object

saveRDS(dds, "./data_output/C42B_DESeq_comp_to_DMSO.Rds")
saveRDS(dds1, "./data_output/C42B_DESeq_comp_to_AbiOla.Rds")
saveRDS(dds2, "./data_output/C42B_DESeq_comp_to_EnzOla.Rds")






