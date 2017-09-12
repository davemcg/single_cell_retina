# Below is a very simple example of using Seurat to analyze a beautiful human
# pancreas dataset generated using inDrop, from the Yanai lab at NYUMC.

library(Seurat)
library(Matrix)

# Sample workflow for leading in files from GEO
######################################################################
all.files <- list.files("~/Downloads/YanaiData/data/")
all.data <- data.frame()
for (i in all.files[1:4]) {
  dat <- read.csv(paste("~/Downloads/YanaiData/data/", i, sep = ""))
  all.data <- rbind(all.data,data.frame(dat))
  print(i)
}

new.data <- t(all.data[, c(-1, -2, -3)])
colnames(new.data) <- all.data[, 1]
pancreas.data <- new.data
pancreas.md <- all.data[, 2:3]
rownames(pancreas.md) <- all.data[, 1]

######################################################################
pancreas.data <- Matrix(pancreas.data, sparse = T)
pancreas <- CreateSeuratObject(raw.data = pancreas.data, min.cells = 3)
pancreas <- AddMetaData(pancreas, metadata = pancreas.md)
pancreas <- FilterCells(pancreas, subset.names = "nGene", low.thresholds = 500,  high.thresholds = Inf)
pancreas <- NormalizeData(pancreas)
pancreas <- FindVariableGenes(pancreas, x.low.cutoff = 0.1)
pancreas <- ScaleData(pancreas, latent.vars = c("orig.ident", "nUMI"), genes.use = pancreas@var.genes, model.use = "negbinom")
pancreas <- RunPCA(pancreas, pcs.compute = 30, weight.by.var = FALSE)
pancreas <- RunTSNE(pancreas, dims.use = 1:19, do.fast = T)
pancreas <- FindClusters(pancreas, reduction.type = "pca", dims.use = 1:19, save.SNN = T)

# color by cluster ID, annotated cluster from the manuscript, or batch
# Can switch the identity class using SetAllIdent if desired
TSNEPlot(pancreas, do.label = T)
TSNEPlot(pancreas, group.by = "assigned_cluster")
TSNEPlot(pancreas, group.by = "orig.ident")

# Find Markers of ductal cell subcluster, using the negative binomial test
# only test genes with a 20% difference in detection rate to speed-up (optional)
ductal.markers <- FindMarkers(pancreas, ident.1 = 5, ident.2 = 12, test.use = "negbinom", min.diff.pct = 0.2)

# Visualize canonical and new markers
FeaturePlot(pancreas, c("GCG", "INS","TFF1","PRSS1","VGF","TRIB3","DDR1","CRYBA2","SLC30A8"),
            cols.use = c("lightgrey","blue"), nCol = 3)

# Can save the object for future loading
save(pancreas, file = "~/Projects/datasets/pancreas.Robj")
