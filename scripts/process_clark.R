# process counts data from https://www.biorxiv.org/content/early/2018/07/30/378950
# blackshaw / goff data
# https://github.com/gofflab/developing_mouse_retina_scRNASeq
# renamed "10x_mouse_retina_development.mtx" to "matrix.mtx"

library(cowplot)
library(Seurat)
library(tidyverse)
library(here)


phenotype_data <- read_csv('/Volumes/Arges/PROJECTS/mcgaughey/SC_RNAseq_outside/clark_blackshaw__mouse_retina/10x_mouse_retina_development_phenotype.csv')
write_tsv(phenotype_data$barcode %>% data.frame(), col_names = F, path = '/Volumes/Arges/PROJECTS/mcgaughey/SC_RNAseq_outside/clark_blackshaw__mouse_retina/barcodes.tsv')
gene_data <- read_csv('/Volumes/Arges/PROJECTS/mcgaughey/SC_RNAseq_outside/clark_blackshaw__mouse_retina/10x_mouse_retina_development_feature.csv')
gene_data %>% head()
write_tsv(gene_data %>% select(id, gene_short_name), col_names = F, path = '/Volumes/Arges/PROJECTS/mcgaughey/SC_RNAseq_outside/clark_blackshaw__mouse_retina/genes.tsv')
retina_raw <- Read10X('/Volumes/Arges/PROJECTS/mcgaughey/SC_RNAseq_outside/clark_blackshaw__mouse_retina/') #, min.cells = 0.001 * ncol(dge_mat), min.genes = 200)
retina <- CreateSeuratObject(retina_raw, min.cells = 0.001 * ncol(retina_raw), min.genes = 200)

mito.genes <- grep(pattern = "^MT-", x = rownames(x = retina@data), value = TRUE)
percent.mito <- colSums(as.matrix(retina@data[mito.genes, ])) / colSums(as.matrix(retina@data))
retina <- AddMetaData(object = retina, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = retina, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)


# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@data.info, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = retina, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = retina, gene1 = "nUMI", gene2 = "nGene")

# We filter out cells that have unique gene counts over 6,000 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate' -Inf and Inf should be used if you don't want a lower or upper
# threshold.
retina <- FilterCells(object = retina, subset.names = c("nGene", "percent.mito"), 
                      low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.05))

# Normalize the data
retina <- NormalizeData(object = retina, normalization.method = "LogNormalize")


# load metadata into seurat object
seurat_data <- phenotype_data %>% data.frame() %>% select(-contains('Patt'))
row.names(seurat_data) <- gsub('-1','',seurat_data$barcode)
seurat_data <- seurat_data[,c(-1,-2)]
retina <- AddMetaData(object = retina, metadata = seurat_data)

#print with macosko clustering
#TSNEPlot(object = retina, do.label = TRUE, group.by = 'Macosko_Clusters')


# save processed data
save(retina, file = '/Volumes/Arges/PROJECTS/mcgaughey/SC_RNAseq_outside/clark_blackshaw__mouse_retina/process_clark.Rdata')
