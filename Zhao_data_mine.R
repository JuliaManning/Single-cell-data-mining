# Data downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169332
# from this paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8170046/

# Install packages
library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)

# Download the data
Chow <- Read10X(data.dir = "C:/Users/julia/OneDrive - University of Birmingham/For_others/Franzi/Zhao_data/Chow_filtered_feature_bc_matrix")
DDC <- Read10X(data.dir = "C:/Users/julia/OneDrive - University of Birmingham/For_others/Franzi/Zhao_data/DDC_filtered_feature_bc_matrix")

# Create seurat objects and trimm for cells with <200 genes and genes expresssed in fewer than 5 cells
Chow <- CreateSeuratObject(counts = Chow, min.cells = 5, project = "Chow")
DDC <- CreateSeuratObject(counts = DDC, min.cells = 5, project = "DDC")

# Merge the two objects
combined <- merge(Chow, y = DDC, add.cell.ids = c("Chow", "DDC"), project = "Combined")

# QC
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Subset as described in paper, gene counts >3000 or <200, and those expressing >10% mitochondrial unique molecular identifier (UMI) 
combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)

# Split the dataset into a list of two seurat objects (stim and CTRL)
Data.list <- SplitObject(combined, split.by = "orig.ident")
Data.list <- lapply(X = Data.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Data.list)

EC.anchors <- FindIntegrationAnchors(object.list = Data.list, anchor.features = features)

EC.combined <- IntegrateData(anchorset = EC.anchors)

DefaultAssay(EC.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
EC.combined <- ScaleData(EC.combined, verbose = FALSE)
EC.combined <- RunPCA(EC.combined, npcs = 100, verbose = FALSE)

ElbowPlot(EC.combined, ndims = 100)

# UMAP and Clustering
EC.combined <- RunUMAP(EC.combined, reduction = "pca", dims = 1:15)
EC.combined <- FindNeighbors(EC.combined, reduction = "pca", dims = 1:15)
EC.combined <- FindClusters(EC.combined, resolution = 0.7)

# Visualization
p1 <- DimPlot(EC.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(EC.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

# To see how the two differnent cell types cluster
DimPlot(EC.combined, reduction = "umap", split.by = "orig.ident", label = TRUE)

# Which markers contibute to differences in clusters
# For cluster 1, can run for all clusters, but know what to use from paper
EC.markers.1 <- FindConservedMarkers(EC.combined, ident.1 =1, grouping.var = "orig.ident", verbose = FALSE)
head(EC.markers.1)

# Downstream analysis with all genes
DefaultAssay(EC.combined) <- "RNA"

# Make dot plot as in paper
classify_genes <- c("Cdh5", "Pecam1", "Tie1", "Tek", "Vwf", "Col4a1", "Cd68", "Adgre1","Lgals3", "Dcn", "Col1a1", "Col3a1", "Myh11", "Acta2", "Tagln")
DotPlot(object = EC.combined, features = classify_genes)

# Feature plots
a <- FeaturePlot(EC.combined, features = c("Cdh5") )
b <- FeaturePlot(EC.combined, features = c("Pecam1") )
c <- FeaturePlot(EC.combined, features = c("Tek") )
d <- FeaturePlot(EC.combined, features = c("Vwf") )
cowplot::plot_grid(a,b,c,d,
                   nrows =2, ncol = 2)

# Remove clusters 6 (fibros), 14 (macros), 15 (SMC), 11 and 13 (non-ECs)

# First, rename the clusters
new.cluster.ids <- c("EC_1", "EC_2", "EC_3", "EC_4", "EC_5", "EC_6", "Fibros",
                     "EC_7", "EC_8", "EC_9", "EC_10", "non-EC", "EC_11", "non-EC", 
                     "Macro", "SMC")
names(new.cluster.ids) <- levels(EC.combined)
EC.combined <- RenameIdents(EC.combined, new.cluster.ids)
DimPlot(EC.combined, reduction = "umap", split.by = "orig.ident", label = TRUE)

# Remove all clusters that aren't ECs
EC <- subset(EC.combined, idents = c("EC_1", "EC_2", "EC_3", "EC_4", "EC_5", "EC_6", 
                                     "EC_7", "EC_8", "EC_9", "EC_10","EC_11"))

DimPlot(EC, reduction = "umap", split.by = "orig.ident", label = TRUE)
FeaturePlot(EC, features = c("Lgals9"), split.by = "orig.ident" )
DimPlot(EC, reduction = "umap",label = TRUE)
FeaturePlot(EC, features = c("Lgals9"))
VlnPlot(EC, features = c("Lgals9"), split.by = "orig.ident",
        pt.size = 0, )
VlnPlot(EC, features = c("Lgals9"),  group.by= "orig.ident",
        pt.size = 0)





