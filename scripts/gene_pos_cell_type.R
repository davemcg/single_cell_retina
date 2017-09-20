# gene+ cell assignment
library(Seurat)
library(tidyverse)
library(here)
# takes the macosko cell classifications and returns, for each gene, the positive counts (for the gene) in each of the 12 major cell types 
load('data/retina_seurat_superSet.Rdata')
source('scripts/macosko_cluster_assignments.R')



all_genes <- retina_superset@data %>% rownames()
gene_type_counts <- data.frame(c('Amacrine cells','Astrocytes','Bipolar cells','Cones','Fibroblasts','Horizontal cells','Microglia','Muller glia','Pericytes','Retinal ganglion cells','Rods','Vascular endothelium'))
colnames(gene_type_counts) <- 'Cell Type'
                               
for (i in all_genes){
  gene_pos <- retina_superset@data[i,] > 0.1
  df <- retina_superset@meta.data[gene_pos,'Cell Type'] %>% table() %>% as.data.frame() %>% select(Freq)
  colnames(df) <- c(i)
  gene_type_counts <- cbind(gene_type_counts, df)
}


row.names(gene_type_counts) <- gene_type_counts$c..Amacrine.cells....Astrocytes....Bipolar.cells....Cones....Fibroblasts...
gene_type_counts <- gene_type_counts[,-1]
save(gene_type_counts, file='data/gene_type_counts.Rdata')
