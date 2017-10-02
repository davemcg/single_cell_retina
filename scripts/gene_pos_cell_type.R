# gene+ cell assignment
library(Seurat)
library(tidyverse)
library(here)
library(RSQLite)
# takes the macosko cell classifications and my processed (Seurat) macosko data and returns for eyeIntegration:
# 1. a sqlite file with the individual sample counts for each of the 12 major cell types and metadata (subset and superset / short and long) and Gene names sql table
# 2. a sql table with the tsne coordinates 
# 3. a sql table with summarized values for ratios of each of the 12 cell types by gene
# 4. a sql table with, for each gene, the percentage of cells (split by cell type) expressing the gene
# 5. a sql table with mean expression of each gene by the 12 cell types and the ranks and deciles for each gene

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
dbWriteTable(sqldb, 'single_cell_gene_counts', gene_count_long_metadata, field.types=NULL, overwrite=TRUE)
dbWriteTable(sqldb, 'single_cell_metadata_long', metadata_long, field.types=NULL)
dbWriteTable(sqldb, 'single_cell_metadata_short', metadata_short, field.types=NULL)
dbWriteTable(sqldb, 'gene_names', gene_names, field.types=NULL)
dbGetQuery(sqldb, "CREATE INDEX GeneName on single_cell_gene_counts(Gene)")

#dbDisconnect(sqldb)

# 2
tsne_coords <-GetDimReduction(object = retina, reduction.type = "tsne", slot = "cell.embeddings") %>% data.frame() %>% rownames_to_column('Cell ID') %>% 
  left_join(metadata_long) %>% filter(!is.na(`Cell Type`))
dbWriteTable(sqldb, 'tsne_coords', tsne_coords, field.types=NULL)


# 3
cell_types <- c('Rods','Bipolar cells','Amacrine cells','Cones','Muller glia','Retinal ganglion cells','Vascular endothelium','Horizontal cells','Microglia','Pericytes','Astrocytes','Fibroblasts')
gene_cell_type_perc <-gene_count_long_metadata %>% 
  group_by(Gene, `Cell Type`) %>% 
  summarise(`Cell Count`=n())%>% 
  mutate(`Percentage Cell Types`=(`Cell Count`/sum(`Cell Count`))*100) %>% 
  left_join(data.frame(cell_types) %>% 
              select(`Cell Type`=cell_types),.,by=c('Cell Type'='Cell Type'))
# calc rank and decile
#gene_cell_type_perc <- gene_cell_type_perc %>% arrange(`Cell Type`, -Percentage) %>% group_by(`Cell Type`) %>% mutate(Rank=row_number()) %>% mutate(Decile=ntile(-Rank, 10))
#dbWriteTable(sqldb, 'gene_cell_type_perc', gene_cell_type_perc, field.types=NULL, overwrite=TRUE)

# 4
cell_type_counts <- metadata_long %>% pull(`Cell Type`) %>% table() %>% data.frame()
colnames(cell_type_counts) <- c('Cell Type', 'Total Count')
cell_types <- c('Rods','Bipolar cells','Amacrine cells','Cones','Muller glia','Retinal ganglion cells','Vascular endothelium','Horizontal cells','Microglia','Pericytes','Astrocytes','Fibroblasts')
cell_type_gene_perc <-gene_count_long_metadata %>% 
  group_by(`Cell Type`, Gene) %>% 
  summarise(`Cell Count`=n())%>% 
  left_join(data.frame(cell_types) %>% 
              select(`Cell Type`=cell_types),.,by=c('Cell Type'='Cell Type')) %>% 
  left_join(., cell_type_counts) %>% 
  mutate(`Percentage Cells` = (`Cell Count`/`Total Count`)*100)
# calc rank and decile
#cell_type_gene_perc <- cell_type_gene_perc %>% arrange(`Cell Type`, -Percentage) %>% group_by(`Cell Type`) %>% mutate(Rank=row_number()) %>% mutate(Decile=ntile(-Rank, 10)) 
#dbWriteTable(sqldb, 'cell_type_perc', gene_cell_type_perc, field.types=NULL, overwrite=TRUE)
# 5
cell_type_ids <- list()
set.seed(1234)
for (i in cell_types){
  cell_type_ids[[i]] <- retina_superset@meta.data %>% data.frame() %>% rownames_to_column('Cell ID') %>% filter(`Cell.Type` == i) %>% pull(`Cell ID`)
}
gene_means_by_type <- data.frame()
for (i in cell_types){
  Mean <- apply(retina_superset@data[,cell_type_ids[[i]]], 1, function(x) mean(x))
  Mean <- data.frame(Mean)
  Mean$Gene <- row.names(retina_superset@data)
  Mean$`Cell Type` <- i
  gene_means_by_type <- bind_rows(gene_means_by_type, Mean)
}
#gene_means_by_type <- gene_means_by_type %>% arrange(`Cell Type`, -Mean) %>% group_by(`Cell Type`) %>% mutate(Rank=row_number()) %>% mutate(Decile=ntile(-Rank, 10)) %>% 
#  select(Gene, `Cell Type`, Mean, Rank, Decile)
#dbWriteTable(sqldb, 'gene_means_by_type', gene_means_by_type, field.types=NULL, overwrite=TRUE)
gene_cell_type_stats <- left_join(gene_cell_type_perc, cell_type_gene_perc, by=c('Cell Type', 'Gene')) %>% 
  right_join(. , gene_means_by_type, by=c('Cell Type', 'Gene'))

## left_join 3,4,5 together

## calc stat (rank/decile) for 3,4,5

# disconnect database
dbDisconnect(sqldb)

# example ggplot
tsne_coords2 <- dbGetQuery(sqldb, 'SELECT * FROM tsne_coords')
ggplot(tsne_coords, aes(x=tSNE_1,y=tSNE_2, colour=`Cell Type`))  + geom_point()
# example labeling of samples with ZFP503 expression
samples_gene_up <- colnames(retina_superset@data[,retina_superset@data['ABCA4',] > 5])
ggplot(tsne_coords, aes(x=tSNE_1,y=tSNE_2, colour=`Cell Type`))  + 
  geom_point(alpha=tsne_coords %>% mutate(Alpha = ifelse(`Cell ID` %in% samples_gene_up, 0.5, 0.03)) %>% pull(Alpha)) + 
  xlab('tSNE 1') + ylab('tSNE 2') 


