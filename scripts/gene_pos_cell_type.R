# gene+ cell assignment
library(Seurat)
library(tidyverse)
library(here)
library(RSQLite)
# takes the macosko cell classifications and my processed (Seurat) macosko data and returns for eyeIntegration:
# 1. a sqlite file with the individual sample counts for each of the 12 major cell types and metadata (subset and superset / short and long)
# 2. a sql table with the tsne coordinates 
# 3. Gene names sql table
load('~/git/single_cell_retina/data/retina_seurat_superSet.Rdata')
load('~/git/single_cell_retina/data/retina_seurat_subSet.Rdata')
source('~/git/single_cell_retina/scripts/macosko_cluster_assignments.R')

## 1
gene_count_long <- retina_superset@data %>% as.matrix() %>% data.frame() %>% rownames_to_column('Gene') %>% gather(`Cell ID`, `Gene Count`, -Gene) %>% filter(`Gene Count` > 0)
metadata_long <-  retina_superset@meta.data %>% as.matrix() %>% data.frame() %>% rownames_to_column('Cell ID') %>% dplyr::select(`Cell ID`, `Cell.Type`, Macosko_Clusters) %>% mutate(`Cell Type` = Cell.Type) %>% select(-Cell.Type)
metadata_short <-  retina@meta.data %>% as.matrix() %>% data.frame() %>% rownames_to_column('Cell ID') %>% dplyr::select(`Cell ID`, `Cell.Type`, Macosko_Clusters) %>% mutate(`Cell Type` = Cell.Type) %>% select(-Cell.Type)
gene_count_long_metadata <- left_join(gene_count_long, metadata_long) %>% filter(!is.na(`Cell Type`))
gene_names <- gene_count_long$Gene %>% unique() %>% sort() %>% data.frame()
colnames(gene_count_long) <- 'Gene'
sqlite_file <- '~/git/Human_eyeIntegration_App/www/single_cell_retina_info.sqlite'
sqldb <- dbConnect(SQLite(), dbname=sqlite_file)
dbWriteTable(sqldb, 'single_cell_gene_counts', gene_count_long_metadata, field.types=NULL)
dbWriteTable(sqldb, 'single_cell_metadata_long', metadata_long, field.types=NULL)
dbWriteTable(sqldb, 'single_cell_metadata_short', metadata_short, field.types=NULL)
dbWriteTable(sqldb, 'gene_names', gene_names, field.types=NULL)
dbGetQuery(sqldb, "CREATE INDEX GeneName on single_cell_gene_counts(Gene)")

#dbDisconnect(sqldb)

# 2
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

# calculate ks p values for each gene against twelve cell types
library(broom)
cell_types <- c('Rods','Bipolar cells','Amacrine cells','Cones','Muller glia','Retinal ganglion cells','Vascular endothelium','Horizontal cells','Microglia','Pericytes','Astrocytes','Fibroblasts')

cell_type_ids <- list()
set.seed(1234)
for (i in cell_types){
  cell_type_ids[[i]] <- retina_superset@meta.data %>% data.frame() %>% rownames_to_column('Cell ID') %>% filter(`Cell.Type` == i) %>% pull(`Cell ID`)
}

gene_test <- 'ZFP503'

cell_type_tester <- function(gene_name, cell_type_ids){
  cell_type_test <- data.frame()
  for (i in cell_types){
    if (sum(retina_superset@data[gene_name,cell_type_ids[[i]]] > 0) > 0){
      test_results <- wilcox.test(retina_superset@data[gene_name,retina_superset@data[gene_name,cell_type_ids[[i]]] > 0], 
                             as.vector(retina_superset@data[,cell_type_ids[[i]]])) %>% 
        tidy() %>% 
        mutate(cell.type=i, gene=gene_name)
        
        cell_type_test <- rbind(test_results, cell_type_test)
    }
  }
  return(cell_type_test)
}

cell_type_tester('ZFP503', cell_type_ids)
