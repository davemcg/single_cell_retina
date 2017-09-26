# process counts data from Macosko et al to

library(cowplot)
library(Seurat)
library(tidyverse)
library(here)


# raw data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63472
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63472&format=file&file=GSE63472%5FP14Retina%5FlogDGE%2Etxt%2Egz
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63472&format=file&file=GSE63472%5FP14Retina%5Fmerged%5Fdigital%5Fexpression%2Etxt%2Egz
# http://mccarrolllab.com/wp-content/uploads/2015/05/retina_clusteridentities.txt

dge <- read_tsv(here('data/GSE63472_P14Retina_merged_digital_expression.txt.gz'))
dge_mat <- dge[,2:ncol(dge)] %>% as.matrix() %>%  Matrix(.,sparse = TRUE)
row.names(dge_mat) <- dge$gene
retina <- CreateSeuratObject(raw.data = dge_mat, min.cells = 0.001 * ncol(dge_mat), min.genes = 200)

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
# 900 Note that low.thresholds and high.thresholds are used to define a
# 'gate' -Inf and Inf should be used if you don't want a lower or upper
# threshold.

# this dramatically reduces the number of cells from ~45,000 to ~13,0000, but clustering is much better
# macosko et al clustered on ~13,000 initially, then did some fancy stuff to project the remaining 35,000 (mostly rod cells)
# onto the mapping
retina_superset <- FilterCells(object = retina, subset.names = c("nGene", "percent.mito"), 
                               low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.05))
retina <- FilterCells(object = retina, subset.names = c("nGene", "percent.mito"), 
                      low.thresholds = c(900, -Inf), high.thresholds = c(6000, 0.05))

# Normalize the data
retina <- NormalizeData(object = retina, normalization.method = "LogNormalize")
retina_superset <- NormalizeData(object = retina_superset, normalization.method = 'LogNormalize')
# Choose gene outliers on mean-variability plot
retina <- FindVariableGenes(object = retina, x.low.cutoff = 0, y.cutoff = 2)
length(x = retina@var.genes)

# Perform negative-binomial regression on the variable genes, this sets their value in retina@scale.data, which is used for PCA/clustering
# We only do this on the variable genes to save time, but you can do this genome-wide
# We treat mitochondrial percentage, batch, and nUMI as confounding variables,
# I am using the variable gene set calculated from the 'tight' FilterCell threshold, as these are a more stable set
retina <- ScaleData(object = retina, vars.to.regress = c("percent.mito", "nUMI","orig.ident"), genes.use = retina@var.genes, model.use = "negbinom")
# also scale on the superset for later classification
# retina_superset <- ScaleData(object = retina_superset, vars.to.regress = c("percent.mito", "nUMI","orig.ident"), genes.use = retina@var.genes, model.use = "negbinom")

# Run PCA with the IRLBA package (iteratively computes the top dimensions, dramatic increase in speed since we only use a fraction of the PCs anyways)
# if you see the warning "did not converge–results might be invalid!; try increasing maxit or fastpath=FALSE", try increasing maxit
retina <- RunPCA(object = retina, pc.genes = retina@var.genes, pcs.compute = 40, pcs.print = 1:10, maxit = 2000, weight.by.var = FALSE)

# elbow at 30
#PCElbowPlot(object = retina, num.pc = 40)
#PrintPCA(object = retina, pcs.print = 1:36)
#PCHeatmap(object = retina, pc.use = 1:12,100)
#PCHeatmap(object = retina, pc.use = 13:24,100)
#PCHeatmap(object = retina, pc.use = 25:36,100)

#VizPCA(object = retina, pcs.use = 1:2)

# save.SNN means that you can easily re-run with different resolution values.
retina <- FindClusters(object = retina, reduction.type = "pca", dims.use = 1:30, resolution = 1.2, save.SNN = TRUE)

# match PCs used to FindClusters above
retina <- RunTSNE(object = retina, dims.use = 1:30, do.fast = TRUE, perplexity=30, check_duplicates=F)
#TSNEPlot(object = retina, do.label = TRUE)

# load macosko cluster assignments into seurat object
source('scripts/macosko_cluster_assignments.R')
retina <- AddMetaData(object = retina, metadata = macosko_cluster, col.name = c('Macosko_Clusters', `Cell Type`))
retina <- AddMetaData(object = retina, metadata = macosko_cluster$`Cell Type`, col.name = 'Cell_Type')
retina_superset <- AddMetaData(object = retina_superset, metadata = macosko_cluster, col.name = c('Macosko_Clusters', `Cell Type`))
retina_superset <- AddMetaData(object = retina_superset, metadata = macosko_cluster$`Cell Type`, col.name = 'Cell_Type')

#print with macosko clustering
#TSNEPlot(object = retina, do.label = TRUE, group.by = 'Macosko_Clusters')


# save processed data
save(retina, file = here("data/retina_seurat_subSet.Rdata"))
save(retina_superset, file = here("data/retina_seurat_superSet.Rdata"))

