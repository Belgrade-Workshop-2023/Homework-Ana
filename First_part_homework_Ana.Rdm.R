library(dplyr)
library(Seurat)
library(patchwork)

pbmc.data <- Read10X(data.dir = "C:/Users/Ana/Desktop/Emil/ORG_D45_filtered_feature_bc_matrix")

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "retina", min.cells = 3, min.features = 100)
pbmc

data <- data.frame(nCount_RNA = c(100, 200, 300), percent.mt = c(5, 10, 15))
library(ggplot2)  
ggplot(data, aes(x = nCount_RNA, y = percent.mt)) +
  geom_point()

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
                        
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 30)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 30)
                        
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
                        
top10 <- head(VariableFeatures(pbmc), 5)
plot2 <- LabelPoints(plot = plot1, points = top5, repel = TRUE)
plot1 + plot2
                        
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
                        
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 100, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:10, cells = 100, balanced = TRUE)
                        
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:17)
ElbowPlot(pbmc)
                        
pbmc <- FindNeighbors(pbmc, dims = 1:17)
pbmc <- FindClusters(pbmc, resolution = 0.5)
                        
head(Idents(pbmc), 5)
                        
pbmc <- RunUMAP(pbmc, dims = 1:17)
DimPlot(pbmc, reduction = "umap")
                        
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.8, logfc.threshold = 0.8)
pbmc.markers 

VlnPlot(pbmc, features = c("CYP1B1", "POU4F2"))
VlnPlot(pbmc, features = c("GRM6", "CALB1"))

