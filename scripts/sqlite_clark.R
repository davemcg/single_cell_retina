# run on biowulf2 with 128GB of memory
library(Seurat)
library(tidyverse)
library(RSQLite)
# takes the clark blackshaw mouse retina cell classifications and returns for eyeIntegration:
# 1. a sqlite file with the individual sample counts for each of the major cell types and metadata
# 3. a sql table with summarized values for ratios of each of the cell types by gene and age
# 4. a sql table with, for each gene, the percentage of cells (split by cell type) expressing the gene
# 5. a sql table with mean expression of each gene by the 12 cell types and the ranks and deciles for each gene

#load('~/git/single_cell_retina/data/retina_seurat_superSet.Rdata') # NOT USING NOW. DATA IS TOO MESSY.
load('/data/mcgaugheyd/projects/nei/mcgaughey/mouse_SC_clarkShaw/process_clark.Rdata')


## 1
gene_count_long <- retina@data %>% as.matrix() %>% data.frame() %>% rownames_to_column('Gene') %>% gather(`Cell ID`, `Gene Count`, -Gene) %>% filter(`Gene Count` > 0)
metadata_long <-  retina@meta.data %>% as.matrix() %>% data.frame() %>% rownames_to_column('Cell ID') %>% dplyr::select(`Cell ID`, age, umap_CellType, new_CellType, umap_coord1, umap_coord2, umap_coord3) %>% mutate(`Cell Type` = new_CellType) %>% select(-new_CellType)

gene_count_long_metadata <- left_join(gene_count_long, metadata_long) %>% filter(!is.na(`Cell Type`))
gene_names <- gene_count_long$Gene %>% unique() %>% sort() %>% data.frame()
colnames(gene_names) <- 'Gene'
sqlite_file <- '/data/mcgaugheyd/single_cell_retina_info_2.sqlite'
sqldb <- dbConnect(SQLite(), dbname=sqlite_file)
dbWriteTable(sqldb, 'clark__SC_gene_counts', gene_count_long_metadata, field.types=NULL, overwrite=TRUE)
dbWriteTable(sqldb, 'clark__SC_metadata_long', metadata_long %>% select(`Cell ID`, Age = age, `Cell Type`, umap_CellType), field.types=NULL)
#dbWriteTable(sqldb, 'single_cell_metadata_short', metadata_short, field.types=NULL)
dbWriteTable(sqldb, 'clark__gene_names', gene_names, field.types=NULL, overwrite=TRUE)
dbSendQuery(sqldb, "CREATE INDEX clark_GeneName on clark__SC_gene_counts(Gene)")

#dbDisconnect(sqldb)

# 2
dbWriteTable(sqldb, 'clark__tsne_coords', metadata_long %>% 
               mutate(umap_coord1 = as.numeric(as.character(umap_coord1)),
                      umap_coord2 = as.numeric(as.character(umap_coord2)),
                      umap_coord3 = as.numeric(as.character(umap_coord3))), field.types=NULL, overwrite=TRUE)


# 3
cell_types <- retina@meta.data$new_CellType %>% unique()
cell_types <- cell_types[!grepl('Doub', cell_types)]
cell_types <- cell_types[!grepl('Red', cell_types)]

gene_cell_type_perc <- gene_count_long_metadata %>% 
  filter(`Cell Type` != 'Doublets',
         `Cell Type` != 'Red Blood Cells') %>% 
  select(Gene, `Cell Type`, Age = age) %>% 
  group_by(Gene, `Cell Type`, Age) %>% 
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
#cell_types <- c('Rods','Bipolar cells','Amacrine cells','Cones','Muller glia','Retinal ganglion cells','Vascular endothelium','Horizontal cells','Microglia','Pericytes','Astrocytes','Fibroblasts')
cell_type_gene_perc <-gene_count_long_metadata %>% 
  select(Gene, `Cell Type`, Age = age) %>% 
  group_by(`Cell Type`, Gene, Age) %>% 
  summarise(`Cell Count`=n())%>% 
  left_join(data.frame(cell_types) %>% 
              select(`Cell Type`=cell_types),.,by=c('Cell Type'='Cell Type')) %>% 
  left_join(., cell_type_counts) %>% 
  mutate(`Percentage Cells` = (`Cell Count`/`Total Count`)*100)
# calc rank and decile
#cell_type_gene_perc <- cell_type_gene_perc %>% arrange(`Cell Type`, -Percentage) %>% group_by(`Cell Type`) %>% mutate(Rank=row_number()) %>% mutate(Decile=ntile(-Rank, 10)) 
#dbWriteTable(sqldb, 'cell_type_perc', gene_cell_type_perc, field.types=NULL, overwrite=TRUE)
# 5
# by type AND age

gene_means_by_type <- gene_count_long_metadata %>% 
						group_by(Gene, `Cell Type`, age) %>% 
						summarise(Mean = mean(`Gene Count`)) %>%
						select(Gene, `Cell Type`, Age = age, Mean)

# join 3,4,5 together
gene_cell_type_stats <- left_join(gene_cell_type_perc, cell_type_gene_perc, by=c('Cell Type', 'Gene', 'Age')) %>% 
  right_join(. , gene_means_by_type, by=c('Cell Type', 'Gene', 'Age')) %>% 
  select(-`Cell Count.y`, -`Total Count`) %>%
  rename(`Cell Count`='Cell Count.x')
gene_cell_type_stats[is.na(gene_cell_type_stats)] <- 0
gene_cell_type_stats <- gene_cell_type_stats %>% left_join(cell_type_counts) %>% select(Gene, Age, `Cell Type`, `Cell Count`, `Total Count`, Mean, `Percentage Cells`, `Percentage Cell Types`)

# calculate rank/decile for 3,4,5
gene_cell_type_stats <- gene_cell_type_stats %>% 
  arrange(`Cell Type`, -`Percentage Cell Types`) %>% group_by(`Cell Type`, Age) %>% mutate(Rank_cell_types=row_number()) %>% mutate(Decile_cell_types=ntile(-Rank_cell_types, 10)) %>% #3
  arrange(`Cell Type`, -`Percentage Cells`) %>% group_by(`Cell Type`, Age) %>% mutate(Rank_cells=row_number()) %>% mutate(Decile_cells=ntile(-Rank_cells, 10)) %>%  #4
  arrange(`Cell Type`, -Mean) %>% group_by(`Cell Type`, Age) %>% mutate(Rank_mean=row_number()) %>% mutate(Decile_mean=ntile(-Rank_mean, 10)) #5

# add new table to db
dbWriteTable(sqldb, 'clark__gene_cell_type_stats', gene_cell_type_stats, field.types=NULL, overwrite=TRUE)
dbSendQuery(sqldb, "CREATE INDEX clark__GeneNameStats on clark__gene_cell_type_stats(Gene)")
# disconnect database
dbDisconnect(sqldb)

# example ggplot
tsne_coords2 <- dbGetQuery(sqldb, 'SELECT * FROM clark__tsne_coords')
ggplot(z, aes(x=umap_coord1,y=umap_coord2, colour=`Cell Type`))  + geom_point() + scale_color_viridis_d()
# example labeling of samples with ZFP503 expression
samples_gene_up <- colnames(retina@data[,retina@data['ABCA4',] > 5])
ggplot(tsne_coords, aes(x=umap_coord1,y=umap_coord2, colour=`Cell Type`))  + 
  geom_point(alpha=tsne_coords %>% mutate(Alpha = ifelse(`Cell ID` %in% samples_gene_up, 0.5, 0.03)) %>% pull(Alpha)) + 
  xlab('tSNE 1') + ylab('tSNE 2') 


