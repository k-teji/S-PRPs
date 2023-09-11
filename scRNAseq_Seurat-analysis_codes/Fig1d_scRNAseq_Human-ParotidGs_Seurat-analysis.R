##' scRNAseq of ParotidGs from Human
##' 
##' 

##=== library ===##
library(tidyverse)
library(data.table)
library(Seurat)

##=== import reference data ===##
## The count matrix and the annotation metadata were extracted from Tabula sapiens h5ad files 
## Import raw count of SalivaryG in Tablua sapiens
Hs_PGs_ref_count <- data.table::fread('Tabula-sapiens__raw_count_Salivary_Gland.csv',
                                      header = T, stringsAsFactors = F, check.names = F)
## Import all metadata in Tabula sapiens
Hs_PGs_ref_metadata <- data.table::fread('Tabula-sapiens_meta_data.csv',
                                         header = T, stringsAsFactors = F, check.names = F)
## Filter metadata by only Salivary Gland data
Hs_PGs_ref_metadata <- Hs_PGs_ref_metadata %>%
  dplyr::filter(., organ_tissue == 'Salivary_Gland' & method == '10X')
head(Hs_PGs_ref_metadata)
unique(Hs_PGs_ref_metadata$anatomical_information)

## Further select metadata of only ParotidGs
Hs_PGs_ref_metadata <- Hs_PGs_ref_metadata %>%
  dplyr::filter(., str_detect(anatomical_information, '[Pp]arotid'))

## Extract count data in only ParotidGs
Hs_PGs_ref_count <- Hs_PGs_ref_count %>%
  dplyr::select(., all_of(c('index', Hs_PGs_ref_metadata$cell_id)))
Hs_PGs_ref_count <- Hs_PGs_ref_count %>%
  tibble::column_to_rownames('index')

## Create Seurat object
Hs_PGs_ref_seu <- CreateSeuratObject(Hs_PGs_ref_count)

##=== reference data prep for mapping ===##
Hs_PGs_ref_metadata <- Hs_PGs_ref_metadata %>%
  tibble::column_to_rownames('cell_id')
Hs_PGs_ref_seu <- AddMetaData(Hs_PGs_ref_seu, metadata = Hs_PGs_ref_metadata)

Hs_PGs_ref_seu <- NormalizeData(Hs_PGs_ref_seu, normalization.method = 'LogNormalize', scale.factor = 10000)
Hs_PGs_ref_seu <- FindVariableFeatures(Hs_PGs_ref_seu, selection.method = 'vst', nfeatures = 2000)
Hs_PGs_ref_seu <- ScaleData(Hs_PGs_ref_seu, features = rownames(Hs_PGs_ref_seu))

Hs_PGs_ref_seu <- RunPCA(Hs_PGs_ref_seu, features = VariableFeatures(Hs_PGs_ref_seu), npcs = 200)
Hs_PGs_ref_seu@reductions
Hs_PGs_ref_seu <- JackStraw(Hs_PGs_ref_seu, num.replicate = 100, dims = 200)
Hs_PGs_ref_seu <- ScoreJackStraw(Hs_PGs_ref_seu, dims = 1:200)
JackStrawPlot(Hs_PGs_ref_seu, dims = 1:200)

npcs <- 56 # pvalue < 0.01

##########################################

## Query data analysis
##=== import query data ===##
## The file 'data/outs/filtered_feature_bc_matrix/' is ours.
seu_count <- Read10X('data/outs/filtered_feature_bc_matrix/')

##=== create seurat object ===##
seu_obj <- CreateSeuratObject(counts = seu_count, min.cells = 3, min.features = 100)
seu_obj

##=== calculate mitochondrial percent ===##
seu_obj[['percent.mt']] <- PercentageFeatureSet(seu_obj, pattern = '^MT-')

##=== QC ===##
VlnPlot(seu_obj, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
plot1 <- FeatureScatter(seu_obj, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
plot2 <- FeatureScatter(seu_obj, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
plot1+plot2
hist(seu_obj@meta.data$nCount_RNA)
hist(seu_obj@meta.data$nFeature_RNA)
hist(seu_obj@meta.data$percent.mt)
summary(seu_obj@meta.data$nCount_RNA)
summary(seu_obj@meta.data$nFeature_RNA)

##=== subset data ===##
seu_obj <- subset(seu_obj, subset = nFeature_RNA >= 100 & percent.mt <= 25)
seu_obj
filter_vln <- VlnPlot(seu_obj, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)

plot1 <- FeatureScatter(seu_obj, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
plot2 <- FeatureScatter(seu_obj, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
plot1+plot2

##=== Normalise data ===##
seu_obj <- NormalizeData(seu_obj, scale.factor = 10000)

##=== Find HVGs ===##
seu_obj <- FindVariableFeatures(seu_obj, selection.method = 'vst', nfeatures = 2000)
top10 <- head(VariableFeatures(seu_obj), 10)
plot1 <- VariableFeaturePlot(seu_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)
plot2

##=== scale data ===##
all_genes <- rownames(seu_obj)
seu_obj <- ScaleData(seu_obj, features = all_genes)

## Mapping onto the reference data to annotate cell types
hs_PGs_anchors <- FindTransferAnchors(reference = Hs_PGs_ref_seu, query = seu_obj, dims = 1:npcs)
prediction_table <- TransferData(anchorset = hs_PGs_anchors, refdata = Hs_PGs_ref_seu@meta.data$free_annotation,
                                 dims = 1:npcs)

##=== Query data analysis ===##
seu_obj <- AddMetaData(seu_obj, metadata = prediction_table)

celltype_label_list <- unique(prediction_table$predicted.id)
celltype_label_list <- celltype_label_list[c(3,11,10,15,5,8,13,7,17,2,1,4,16,9,18,6,19,14,12)]
celltype_label_list
seu_obj@active.ident <- factor(seu_obj@meta.data$predicted.id, levels = celltype_label_list)
levels(seu_obj)

## Run PCA
seu_obj <- RunPCA(seu_obj, features = VariableFeatures(seu_obj), npcs = 200)
DimHeatmap(seu_obj, dims = 1:5, cells = 500, balanced = T)

## JackStrawPlot
seu_obj <- JackStraw(seu_obj, num.replicate = 100, dims = 200)
seu_obj <- ScoreJackStraw(seu_obj, dims = 1:200)
JackStrawPlot(seu_obj, dims = 1:200)

npcs_q <- 49 # p val < 0.01

## UMAP
seu_obj <- RunUMAP(seu_obj, reduction = 'pca', dims = 1:npcs_q)
DimPlot(seu_obj, label = T) + NoLegend()

## save 
saveRDS(seu_obj, 'Hs_PGs_seurat-primary-analysis.rds')
