library(Seurat)
nk.counts <- ReadMtx(mtx = "./data/raw/matrix.mtx", cells = "./data/raw/barcodes.tsv", 
                  features = "./data/raw/genes.tsv", cell.sep = "\t", feature.sep = "\t")
seuratobj <- CreateSeuratObject(counts = nk.counts, min.cells = 3, min.features = 200, project= "10X_NK")

seuratobj[["percent.mt"]] <- PercentageFeatureSet(seuratobj, pattern = "^MT-")
rm(nk.counts)
    
# QC Visualization
VlnPlot(seuratobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(seuratobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seuratobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filtering
seuratobj <- subset(seuratobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# NORMALIZATION
seuratobj <- NormalizeData(seuratobj)

# HVG Selection
seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst", nfeatures = 1500)

top10 <- head(VariableFeatures(seuratobj), 10)

# Plot HVGs
p1 <- VariableFeaturePlot(seuratobj)
LabelPoints(plot = p1, points = top10, repel = TRUE)

# Scaling
all.genes <- rownames(seuratobj)
seuratobj <- ScaleData(seuratobj, features = all.genes)

# PCA
seuratobj <- RunPCA(seuratobj, features = VariableFeatures(object = seuratobj))

# PCA loadings with elbow plot
ElbowPlot(seuratobj)

seuratobj <- FindNeighbors(seuratobj, dims = 1:15)
seuratobj <- FindClusters(seuratobj)

# UMAP
seuratobj <- RunUMAP(seuratobj, dims = 1:15)
DimPlot(seuratobj, reduction = "umap")

# Check NCAM1 expression
VlnPlot(seuratobj, features = "NCAM1")

ncol(seuratobj)
## 8302

# Number of NCAM1 expressing cells
ncol(subset(seuratobj, subset = NCAM1 > 1))
## 322 out of 8302 ~ 4 %
ncol(subset(seuratobj, subset = NCAM1 > 0))
## 322 out of 8302 


FeaturePlot(seuratobj, "GNLY")
FeaturePlot(seuratobj, "NKG7")
ncol(subset(seuratobj, subset = GNLY > 0))
## 8249
ncol(subset(seuratobj, subset = NKG7 > 0))
## 8234
