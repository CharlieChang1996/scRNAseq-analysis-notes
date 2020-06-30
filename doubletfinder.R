library(Seurat)
library(dplyr)
library(patchwork)
library(umap)
library(DoubletFinder)

##Doublet Removal##
#pK Identification (no ground-truth)
#select pc according to elbow plot
sweep.res.list_seurat_obj <- paramSweep_v3(seurat_obj, PCs = 1:15, sct = FALSE)
sweep.stats_seurat_obj <- summarizeSweep(sweep.res.list_seurat_obj, GT = FALSE)
bcmvn_seurat_obj <- find.pK(sweep.stats_seurat_obj)
# set pk according to the pk plot
nExp_poi <- round(0.05*length(seurat_obj$orig.ident))  ## Assuming 5% doublet formation rate(set your rate according to 10X manual)
seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:15, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
png("ump_doublet.png")
show(DimPlot(seurat_obj, reduction = "umap",group.by = "DF.classifications_0.25_0.06_540"))
dev.off()
# filter doublets
seurat_obj <- subset(seurat_obj, subset = DF.classifications_0.25_0.06_540 =="Singlet")
