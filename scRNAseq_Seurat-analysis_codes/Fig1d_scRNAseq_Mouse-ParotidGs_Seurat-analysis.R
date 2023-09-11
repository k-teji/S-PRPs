##' scRNAseq of ParotidGs from 4 weeks old male mouse
##' 

##=== library ===##
library(tidyverse)
library(data.table)
library(Seurat)

#=== Reference data analysis ===##
## Import dataset
## The file "Parotid scRNAseq CTRL vs IR integrated (annotated).rds" was obtained from 
## the paper 'scRNAseq of healthy and irradiated mouse parotid glands highlights crosstalk between immune and secretory cells during chronic injury'
mm_PGs_ref <- readRDS('Parotid scRNAseq CTRL vs IR integrated (annotated).rds')
mm_PGs_ref@meta.data$Treatment %>% unique()

## subset data for use of CTRL data 
mm_PGs_ref <- subset(mm_PGs_ref, subset = Treatment == 'CTRL')

## data preprocessing
mm_PGs_ref <- NormalizeData(mm_PGs_ref)
mm_PGs_ref <- FindVariableFeatures(mm_PGs_ref)
mm_PGs_ref <- ScaleData(mm_PGs_ref, features = rownames(mm_PGs_ref))
mm_PGs_ref <- RunPCA(mm_PGs_ref, npcs = 200)

mm_PGs_ref <- JackStraw(mm_PGs_ref, dims = 200)
mm_PGs_ref <- ScoreJackStraw(mm_PGs_ref, dims = 1:200)
JackStrawPlot(mm_PGs_ref, dims = 1:200)

npcs <- 40 # pval < 0.01

## import query
## The data of 'data/outs/filtered_feature_bc_matrix/' is ours.
## Import data and create Seurat objects
mm_PGs_query <- Read10X('data/outs/filtered_feature_bc_matrix/')
mm_PGs_query <- CreateSeuratObject(mm_PGs_query)
mm_PGs_query[['percent.mt']] <- PercentageFeatureSet(mm_PGs_query, pattern = '^[Mm][Tt]-')
VlnPlot(mm_PGs_query, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)

## Subset data for QC
mm_PGs_query <- subset(mm_PGs_query, subset = nFeature_RNA > 200 & percent.mt < 5)
VlnPlot(mm_PGs_query, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
mm_PGs_query

## data preprocessing
mm_PGs_query <- NormalizeData(mm_PGs_query)
mm_PGs_query <- FindVariableFeatures(mm_PGs_query)

## mapping onto reference data to annotate cell types
mm_PGs_anchor <- FindTransferAnchors(reference = mm_PGs_ref, query = mm_PGs_query, dims = 1:npcs)
predicton_table <- TransferData(anchorset = mm_PGs_anchor, refdata = mm_PGs_ref@meta.data$CellType)
hist(predicton_table$prediction.score.max)
table(predicton_table$prediction.score.max >= 0.5)

## add metadata of cell type prediction tables
mm_PGs_query <- AddMetaData(object = mm_PGs_query, metadata = predicton_table)

## UMAP analysis
mm_PGs_query <- ScaleData(mm_PGs_query, features = rownames(mm_PGs_query))
mm_PGs_query <- RunPCA(mm_PGs_query, npcs = 200)
mm_PGs_query <- JackStraw(mm_PGs_query, dims = 80)
mm_PGs_query <- ScoreJackStraw(mm_PGs_query, dims = 1:80)
JackStrawPlot(mm_PGs_query, dims = 1:80)

npcs_q <- 58 # p value < 0.01

mm_PGs_query <- RunUMAP(mm_PGs_query, dims = 1:npcs_q)
mm_PGs_query@active.ident <- factor(mm_PGs_query@meta.data$predicted.id, 
                                    levels = c('Acinar', 'Striated duct', 'Intercalated duct', 'Stromal',
                                               'Etv1+', 'Endothelial', 'Myoepithelial', 'B-cells',
                                               'T-cells CD4+', 'T-cells CD8+', 'T-cells CD4+CD8+', 'T-cells FoxP3+', 
                                               'T-cells Cxcr6+', 'NK cells', 'Macrophages', 'DCs'))

umap_p <- DimPlot(mm_PGs_query, label.size = 5) +
  NoLegend()
umap_p <- LabelClusters(umap_p, id = "ident", color = 'black', size = 5, repel = T,  box.padding = 1)

mm_PGs_query@meta.data$predicted.id <- factor(mm_PGs_query@meta.data$predicted.id,
                                              levels = c('Acinar', 'Striated duct', 'Intercalated duct', 'Stromal',
                                                         'Etv1+', 'Endothelial', 'Myoepithelial', 'B-cells',
                                                         'T-cells CD4+', 'T-cells CD8+', 'T-cells CD4+CD8+', 'T-cells FoxP3+', 
                                                         'T-cells Cxcr6+', 'NK cells', 'Macrophages', 'DCs'))

## save data
saveRDS(mm_PGs_query, 'Mm_PGs_4wks_scRNAseq_seurat-analysis.rds')

