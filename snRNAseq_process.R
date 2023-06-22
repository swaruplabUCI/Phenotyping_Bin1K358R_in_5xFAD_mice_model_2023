#================================== NOTE ==================================
# Project from Vivek
# NOTE MODEL_AD
# Bin1_5xFAD_snRNAseq_Nov2022

# conda activate scRNAnATAC_R

###==============================
# NOTE Dec 05, 2022
###==============================
library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
theme_set(theme_cowplot())


# NOTE loading data
data_dir <- '/dfs7/swaruplab/shared_lab/Collaborations/MODEL_AD/Bin1_5xFAD_snRNAseq/cellranger_output/'
out_dir <- '/Collaborations/MODEL_AD/Bin1_5xFAD_snRNAseq_Nov2022/'

NucSeq <- readRDS(paste0(data_dir, 'NucSeq_final.rds'))

head(NucSeq@meta.data)

table(NucSeq@meta.data$orig.ident)

table(NucSeq$Genotype)

table(NucSeq@meta.data$Sample.ID[NucSeq$Genotype == "5xFAD_WT"])
table(NucSeq@meta.data$Sample.ID[NucSeq$Genotype == "BIN1HO"])
table(NucSeq@meta.data$Sample.ID[NucSeq$Genotype == "5xFAD_HEMI"])
table(NucSeq@meta.data$Sample.ID[NucSeq$Genotype == "5xFAD_HEMI-BIN1HO"])


# NOTE checking QC
dir.create(paste0(out_dir, 'Figures/'))
fig_dir <- paste0(out_dir, 'Figures/')

pdf(paste0(fig_dir, 'BeforeAnnotation/', 'ViolinPlot_feature_Cluster.pdf'),height=10,width=8)
VlnPlot(NucSeq, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="Sample.ID", ncol = 1,pt.size=0)
dev.off()

# NOTE checking UMAP
png(paste0(fig_dir, 'BeforeAnnotation/', 'basic_umap_clusters.png'), width=10, height=8, res=200, units='in')
DimPlot(NucSeq, reduction = "umap", group.by='seurat_clusters', label=TRUE) + ggtitle('UMAP colored by seurat clusters')
  # + umap_theme + NoLegend()
dev.off()


png(paste0(fig_dir, 'BeforeAnnotation/', 'basic_umap_clusters_groupbySample.ID.png'), width=10, height=8, res=200, units='in')
DimPlot(NucSeq, reduction = "umap", group.by='Sample.ID', label=TRUE) + ggtitle('UMAP colored by Sample.ID')
  # + umap_theme + NoLegend()
dev.off()

png(paste0(fig_dir, 'BeforeAnnotation/', 'basic_umap_clusters_sepbySample.ID.png'), width=80, height=8, res=200, units='in')
DimPlot(NucSeq, reduction = "umap", split.by='Sample.ID', label=TRUE) + ggtitle('UMAP colored split.by Sample.ID')
  # + umap_theme + NoLegend()
dev.off()


####################################
# Normalization and Harmony on batch information
####################################
# Normalize and FindVariableFeatures
NucSeq <- NormalizeData(NucSeq, normalization.method = "LogNormalize", scale.factor = 10000)

# NucSeq <- FindVariableFeatures(NucSeq, selection.method = "vst", nfeatures = 2000)
NucSeq <- FindVariableFeatures(
  NucSeq,
  selection.method = "vst",
  nfeatures = 4000
)
# scale data:
NucSeq <- ScaleData(NucSeq, features = rownames(NucSeq))
# seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
#?ScaleData

NucSeq <- RunPCA(
  NucSeq,
  features = VariableFeatures(object = NucSeq),
  npcs=50
)

# Harmony
library(harmony)
# Harmony with two or more covariates
NucSeq <- NucSeq %>%
  RunHarmony("sample", assay.use="RNA")

# UMAP and clustering with harmonized PCs
NucSeq <- RunUMAP(NucSeq, reduction='harmony', dims = 1:30)
NucSeq <- FindNeighbors(NucSeq, reduction='harmony')
NucSeq <- FindClusters(NucSeq, resolution = 0.8)

#=====================
# NOTE UMAPs after harmony
# NOTE checking UMAP
png(paste0(fig_dir, 'BeforeAnnotation/', 'basic_umap_clusters_harmony.png'), width=20, height=16, res=600, units='in')
DimPlot(NucSeq, reduction = "umap", group.by='seurat_clusters', label=TRUE) + ggtitle('UMAP colored by seurat clusters after harmony')
  # + umap_theme + NoLegend()
dev.off()


png(paste0(fig_dir, 'BeforeAnnotation/', 'basic_umap_clusters_groupbySample_harmony.png'), width=10, height=8, res=200, units='in')
DimPlot(NucSeq, reduction = "umap", group.by='sample', label=TRUE) + ggtitle('UMAP colored by Sample.ID after harmony')
  # + umap_theme + NoLegend()
dev.off()

png(paste0(fig_dir, 'BeforeAnnotation/', 'basic_umap_clusters_sepbySample_harmony.png'), width=80, height=8, res=200, units='in')
DimPlot(NucSeq, reduction = "umap", split.by='sample', label=TRUE) + ggtitle('UMAP colored split.by Sample.ID after harmony')
  # + umap_theme + NoLegend()
dev.off()

# NOTE saving harmonized NucSeq
saveRDS(NucSeq, file = paste0(out_dir, 'Data/', 'NucSeq_harmonized.rds'))


# set up list of canonical cell type markers
canonical_markers <- list(
  'Astrocyte' = c('Gfap', 'Aqp4', 'Slc1a2'),
  'Pan-neuronal' = c('Snap25', 'Syt1'),
  'Excitatory Neuron' = c('Slc17a7', 'Satb2'),
  'Inhibitory Neuron' = c('Gad1', 'Gad2'),
  'Microglia' = c('Csf1r', 'Cd74', 'P2ry12'),
  'Oligodendrocyte' = c('Mobp', 'Mbp', 'Mog'),
  'OPC' = c('Pdgfra', 'Cspg4'),
  'Pericytes_Endothelia'=c('Pdgfrb', 'Cldn5', 'Flt1')
)

pdf(paste0(fig_dir, 'BeforeAnnotation/', 'ViolinPlot_gene_Clusters_harmonized.pdf'), width=30, height=15)
VlnPlot(NucSeq, features = unlist(canonical_markers), pt.size=0)
dev.off()

# plot heatmap:
library(viridis)
# create feature plots, cutoff expression values for the 98th and 99th percentile
plot_list <- FeaturePlot(
  NucSeq,
  features=unlist(canonical_markers),
  combine=TRUE,
  # max.cutoff='q98'
) + NoLegend()

png(paste0(fig_dir, 'BeforeAnnotation/', 'basic_canonical_marker_featurePlot_harmonized.png'), width=25, height=15, units='in', res=300)
plot_list
# cowplot::plot_grid(plotlist=plot_list, ncol=4)
dev.off()


####################################
# Cluster annotation
####################################
cluster_annotations <- list(
    "0" = "Oligodendrocytes",
    "1" = "Inhibitory_Neuron",
    "2" = "Inhibitory_Neuron",
    "3" = "Excitatory_Neuron",
    "4" = "Excitatory_Neuron",
    "5" = "Inhibitory_Neuron",
    "6" = "Inhibitory_Neuron",
    "7" = "Oligodendrocytes",
    "8" = "Excitatory_Neuron",
    "9" = "Excitatory_Neuron",
    "10" = "Astrocyte",
    "11" = "Inhibitory_Neuron",
    "12" = "Inhibitory_Neuron",
    "13" = "Excitatory_Neuron",
    "14" = "Inhibitory_Neuron",
    "15" = "Microglia",
    "16" = "OPC",
    "17" = "Excitatory_Neuron",
    "18" = "Inhibitory_Neuron",
    "19" = "Pericytes_Endothelia",
    "20" = "Excitatory_Neuron",
    "21" = "Excitatory_Neuron",
    "22" = "Inhibitory_Neuron"
)


NucSeq$Celltypes <- unlist(cluster_annotations[NucSeq$seurat_clusters])

table(NucSeq$Celltypes)


# Idents(seurat_obj) <- seurat_obj$seurat_clusters
p_celltypes <- DimPlot(NucSeq, reduction = "umap", group.by='Celltypes', label = TRUE, repel = TRUE, raster=FALSE) + ggtitle("Celltypes")
png(paste0(fig_dir, 'AfterAnnotation/', 'umap_NucSeq_Harmony_Celltypes_0120_2023.png'), width=12, height=10, units='in', res=600)
p_celltypes
dev.off()

p_celltypes <- DimPlot(NucSeq, reduction = "umap", group.by='Celltypes', label = FALSE, repel = TRUE, raster=FALSE) + ggtitle("Celltypes")
png(paste0(fig_dir, 'AfterAnnotation/', 'umap_NucSeq_Harmony_Celltypes_nolable_0120_2023.png'), width=12, height=10, units='in', res=600)
p_celltypes #+ NoLegend()
dev.off()

NucSeq$previous_Celltypes_clustering <- NucSeq$pervious_Celltypes_clustering
NucSeq$pervious_Celltypes_clustering <- NULL

table(NucSeq$previous_Celltypes_clustering)
head(NucSeq@meta.data)


# NOTE Saving without c("Doublets", "Unknown")
saveRDS(NucSeq, file = paste0(out_dir, 'Data/', 'NucSeq_Updated_0120_2023.rds'))
