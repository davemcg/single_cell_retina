# load Macosko cluster assignments for each cell
macosko_cluster = read_tsv(here('data/retina_clusteridentities.txt'), col_names = F) %>% data.frame()
row.names(macosko_cluster) <- macosko_cluster$X1
macosko_cluster <- macosko_cluster %>% select(X1, X2)
colnames(macosko_cluster) <- c('id','Macosko_Clusters')
macosko_cluster$Macosko_Clusters <- as.integer(as.character(macosko_cluster$Macosko_Clusters))
# https://github.com/olgabot/macosko2015/blob/master/notebooks/02_make_celltype_metadata.ipynb
naming <- data.frame(rbind(c(1, 'Horizontal cells'),
                           c(2, 'Retinal ganglion cells'),
                           c(3, 'Amacrine cells'),
                           c(4, 'Amacrine cells'),
                           c(5, 'Amacrine cells'),
                           c(6, 'Amacrine cells'),
                           c(7, 'Amacrine cells'),
                           c(8, 'Amacrine cells'),
                           c(9, 'Amacrine cells'),
                           c(10, 'Amacrine cells'),
                           c(11, 'Amacrine cells'),
                           c(12, 'Amacrine cells'),
                           c(13, 'Amacrine cells'),
                           c(14, 'Amacrine cells'),
                           c(15, 'Amacrine cells'),
                           c(16, 'Amacrine cells'),
                           c(17, 'Amacrine cells'),
                           c(18, 'Amacrine cells'),
                           c(19, 'Amacrine cells'),
                           c(20, 'Amacrine cells'),
                           c(21, 'Amacrine cells'),
                           c(22, 'Amacrine cells'),
                           c(23, 'Amacrine cells'),
                           c(24, 'Rods'),
                           c(25, 'Cones'),
                           c(26, 'Bipolar cells'),
                           c(27, 'Bipolar cells'),
                           c(28, 'Bipolar cells'),
                           c(29, 'Bipolar cells'),
                           c(30, 'Bipolar cells'),
                           c(31, 'Bipolar cells'),
                           c(32, 'Bipolar cells'),
                           c(33, 'Bipolar cells'),
                           c(34, 'Muller glia'),
                           c(35, 'Astrocytes'),
                           c(36, 'Fibroblasts'),
                           c(37, 'Vascular endothelium'),
                           c(38, 'Pericytes'),
                           c(39, 'Microglia')))
colnames(naming) <- c('Macosko_Clusters', 'Cell Type')
naming$Macosko_Clusters <- as.integer(as.character(naming$Macosko_Clusters))
macosko_cluster <- macosko_cluster %>% left_join(.,naming)
row.names(macosko_cluster) <- macosko_cluster$id
macosko_cluster <- macosko_cluster %>% select(Macosko_Clusters, `Cell Type`)