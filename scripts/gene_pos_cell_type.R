# gene+ cell assignment
library(Seurat)
library(tidyverse)
library(here)
library(RSQLite)
# takes the macosko cell classifications and my processed (Seurat) macosko data and returns:
# 1. for each gene, the positive counts (for the gene) in each of the 12 major cell types (gene_type_counts)
# 2. a sqlite file with the individual sample counts for each of the 12 major cell types
# 3. a file with the tsne coordinates 
load('data/retina_seurat_superSet.Rdata')
load('data/retina_seurat_subSet.Rdata')
source('scripts/macosko_cluster_assignments.R')


## 1
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


## 2
gene_count_long <- retina_superset@data %>% as.matrix() %>% data.frame() %>% rownames_to_column('Gene') %>% gather(`Cell ID`, `Gene Count`, -Gene) %>% filter(`Gene Count` > 0)
metadata_long <-  retina_superset@meta.data %>% as.matrix() %>% data.frame() %>% rownames_to_column('Cell ID') %>% dplyr::select(`Cell ID`, `Cell.Type`, Macosko_Clusters) %>% mutate(`Cell Type` = Cell.Type) %>% select(-Cell.Type)
metadata_short <-  retina@meta.data %>% as.matrix() %>% data.frame() %>% rownames_to_column('Cell ID') %>% dplyr::select(`Cell ID`, `Cell.Type`, Macosko_Clusters) %>% mutate(`Cell Type` = Cell.Type) %>% select(-Cell.Type)
gene_count_long_metadata <- left_join(gene_count_long, metadata_long) %>% filter(!is.na(`Cell Type`))
sqlite_file <- '~/git/Human_eyeIntegration_App/www/single_cell_retina_info.sqlite'
sqldb <- dbConnect(SQLite(), dbname=sqlite_file)
dbWriteTable(sqldb, 'single_cell_gene_counts', gene_count_long_metadata, field.types=NULL)
dbWriteTable(sqldb, 'single_cell_metadata_long', metadata_long, field.types=NULL)
dbWriteTable(sqldb, 'single_cell_metadata_short', metadata_short, field.types=NULL)
dbGetQuery(sqldb, "CREATE INDEX GeneName on single_cell_gene_counts(Gene)")

dbDisconnect(sqldb)

# 3
tsne_coords <-GetDimReduction(object = retina, reduction.type = "tsne", slot = "cell.embeddings") %>% data.frame() %>% rownames_to_column('Cell ID') %>% 
  left_join(metadata_long) %>% filter(!is.na(`Cell Type`))
dbWriteTable(sqldb, 'tsne_coords', tsne_coords, field.types=NULL)
dbDisconnect(sqldb)


# example ggplot
tsne_coords2 <- dbGetQuery(sqldb, 'SELECT * FROM tsne_coords')
ggplot(tsne_coords, aes(x=tSNE_1,y=tSNE_2, colour=`Cell Type`))  + geom_point()
# example labeling of samples with ZFP503 expression
samples_gene_up <- colnames(retina_superset@data[,retina_superset@data['ABCA4',] > 5])
ggplot(tsne_coords, aes(x=tSNE_1,y=tSNE_2, colour=`Cell Type`))  + 
  geom_point(alpha=tsne_coords %>% mutate(Alpha = ifelse(`Cell ID` %in% samples_gene_up, 0.5, 0.03)) %>% pull(Alpha)) + 
  xlab('tSNE 1') + ylab('tSNE 2') 

