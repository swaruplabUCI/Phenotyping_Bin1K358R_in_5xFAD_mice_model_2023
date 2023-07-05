#================================== NOTE ==================================
# Project from Vivek
# NOTE MODEL_AD
# Bin1_5xFAD_snRNAseq_Nov2022

# conda activate scRNAnATAC_R

library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(UCell)
library(scCustomize)
library(viridis)
theme_set(theme_cowplot())
library(ggplot2)
library(ggrastr)
#
library(RColorBrewer)


umap_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank()
)


###=================================
## NOTE Date: 06/24/2023
###=================================
##########################################################################
# NOTE OPC + ODC
####################################
###=================================
# NOTE loading PATH
data_dir <- '/dfs7/swaruplab/zechuas/Collaborations/MODEL_AD/Bin1_5xFAD_snRNAseq_Nov2022/'
out_dir <- '/dfs7/swaruplab/zechuas/Collaborations/MODEL_AD/Bin1_5xFAD_snRNAseq_Nov2022/Analysis/Subcluster_analysis/'
LT_fig_dir <- paste0(data_dir, 'Analysis/Label_transfer_Rosenberg_etal_2018/')
#
fig_dir <- paste0(data_dir, 'Figures/')

NucSeq <- readRDS(paste0(data_dir, 'Data/', 'NucSeq_Updated_w_rosenberg_subclusterNewNames_06222023.rds'))

##===============================
## NOTE change DefaultAssay
# DefaultAssay(NucSeq) <- 'prediction.score.celltype'
DefaultAssay(NucSeq) <- 'RNA'


#=================
# NOTE loading
table(NucSeq@meta.data$Sample.ID)

NucSeq$Genotype <- factor(
  as.character(NucSeq$Genotype),
  levels = c('5xFAD_WT', 'BIN1HO', '5xFAD_HEMI', '5xFAD_HEMI-BIN1HO') # An example
)

##########################################################################
# NOTE Oligodendrocytes lineage
## NOTE ## Park_NatureComm_2023
library(dplyr)
library(stringr)

Park_NatureComm_markers <- read.csv("/dfs7/swaruplab/zechuas/Collaborations/MODEL_AD/Bin1_5xFAD_snRNAseq_Nov2022/Data/Park_NatureComm_2023_data/FromGitHub_mouse_cell_type_markers.csv")
# Park_NatureComm_markers <- sapply(Park_NatureComm_markers[c("COP", "MFOL", "MOL1", "MOL2", "MOL3", "MOL4", "DAO")], tolower)
Park_NatureComm_markers <- as.data.frame(Park_NatureComm_markers)

colnames(Park_NatureComm_markers)
head(Park_NatureComm_markers)
class(Park_NatureComm_markers)

# Function to convert uppercase letters to lowercase, excluding the first character
convertToLowercase <- function(string) {
  first_char <- substr(string, 1, 1)  # Preserve the first character
  rest_chars <- substring(string, 2)  # Get the remaining characters
  lowercase_rest <- tolower(rest_chars)  # Convert the remaining characters to lowercase
  return(paste(first_char, lowercase_rest, sep = ""))  # Concatenate the preserved and converted characters
}

COPtoODC_markers <- as.character()

# Apply the conversion function to the 'Name' column
col_names <- colnames(Park_NatureComm_markers)
for (i in col_names){
  print(i)
  print(paste("Working on", i))
  Park_NatureComm_markers[[i]] <- sapply(Park_NatureComm_markers[[i]], convertToLowercase)

  COPtoODC_markers <- c(COPtoODC_markers, Park_NatureComm_markers[[i]])

}

length(COPtoODC_markers)
# [1] 210

# COPtoODC_markers <-c(COPtoODC_markers, "Bin1")

##########################################################################
# NOTE Oligodendrocytes + OPC lineage
# Oligo  MFOL -- Oligo MOL -- Oligo NFOL -- OPC

# NucSeq_OPC_ODC <- NucSeq %>% subset(subclusters %in% c('OPC', 'COP', 'Oligo NFOL', 'Oligo  MFOL', 'Oligo MOL'))
NucSeq_ODC <- NucSeq %>% subset(subclusters %in% c('COP', 'Oligo NFOL', 'Oligo  MFOL', 'Oligo MOL'))
table(NucSeq_ODC@meta.data$subclusters)
#

cur_RNA <- NucSeq_ODC
# cur_RNA <- NucSeq_ODC

library(hdWGCNA)
# get expressed_genes at least 25%
cur_RNA <- SetupForWGCNA(
  cur_RNA,
  gene_select = "fraction",
  fraction = 0.05,
  group.by = 'subclusters',
  wgcna_name = "OPC_ODC"
)
length(GetWGCNAGenes(cur_RNA))
expressed_genes <- GetWGCNAGenes(cur_RNA)

# expressed_genes <- rownames(cur_RNA)

valid_COPtoODC_markers <- COPtoODC_markers[COPtoODC_markers %in% expressed_genes]
length(valid_COPtoODC_markers)

length(rownames(NucSeq_ODC))

NucSeq_ODC_clusters <- NormalizeData(NucSeq_ODC, normalization.method = "LogNormalize", scale.factor = 10000)
NucSeq_ODC_clusters <- FindVariableFeatures(
  NucSeq_ODC_clusters,
  selection.method = "vst",
  nfeatures = 500
)

# VariableFeatures(object = NucSeq_ODC_clusters)
length(unique(c(valid_COPtoODC_markers, VariableFeatures(object = NucSeq_ODC_clusters))))

NucSeq_ODC_clusters <- ScaleData(NucSeq_ODC_clusters, features = unique(c(valid_COPtoODC_markers, VariableFeatures(object = NucSeq_ODC_clusters))))
NucSeq_ODC_clusters <- RunPCA(
  NucSeq_ODC_clusters,
  features = valid_COPtoODC_markers,
  npcs=30
)

#=====================
NucSeq_ODC_clusters <- FindNeighbors(NucSeq_ODC_clusters, reduction='pca')
NucSeq_ODC_clusters <- FindClusters(NucSeq_ODC_clusters, resolution = 0.5)
NucSeq_ODC_clusters <- RunUMAP(NucSeq_ODC_clusters, reduction='pca', dims = 1:30)

table(NucSeq_ODC_clusters$seurat_clusters)
# length(rownames(NucSeq_ODC_clusters))


DimPlot(NucSeq_ODC_clusters, reduction = "umap", group.by='seurat_clusters', label=TRUE, raster = TRUE) + ggtitle('UMAP colored by seurat clusters') # + rasterise(geom_point(), dpi = 600, scale = 2)

DimPlot(NucSeq_ODC_clusters, reduction = "umap", group.by='subclusters', label=TRUE, raster = TRUE) + ggtitle('UMAP colored by subclusters') # + rasterise(geom_point(), dpi = 600, scale = 2)


#=====================
# NOTE # add gene score
#=====================
signatures <- list(
  COP = Park_NatureComm_markers$COP,
  MOL1 = Park_NatureComm_markers$MOL1,
  MOL2 = Park_NatureComm_markers$MOL2,
  MOL3 = Park_NatureComm_markers$MOL3,
  MOL4 = Park_NatureComm_markers$MOL4,
  DAO = Park_NatureComm_markers$DAO
)

NucSeq_ODC_clusters <- AddModuleScore_UCell(NucSeq_ODC_clusters, features = signatures, name = NULL, ncores = 4)

pal <- viridis(n = 5, option = "D")

plot_list <- FeaturePlot_scCustom(
  NucSeq_ODC_clusters,
  features=names(signatures),
  colors_use = pal,
  max.cutoff='q95',
  min.cutoff='q05'
)

png(paste0(out_dir, 'ODC_COP/ModuleScore/', 'ModuleScore_UCell_ODC_subcluster_06242023.png'), width=10, height=8, units='in', res=300)
# FeaturePlot(NucSeq_Oligodendrocytes, reduction = "umap", features = names(signatures), ncol = 3, order = T)
# cowplot::plot_grid(plotlist=plot_list, ncol=4)
plot_list
dev.off()



# seurat_clusters
VlnCOP <- VlnPlot(object = NucSeq_ODC_clusters, features = "COP", group.by = 'seurat_clusters') #
VlnMFOL <- VlnPlot(object = NucSeq_ODC_clusters, features = "MFOL", group.by = 'seurat_clusters')
VlnMOL1 <- VlnPlot(object = NucSeq_ODC_clusters, features = "MOL1", group.by = 'seurat_clusters')
VlnMOL2 <- VlnPlot(object = NucSeq_ODC_clusters, features = "MOL2", group.by = 'seurat_clusters')
VlnMOL3 <- VlnPlot(object = NucSeq_ODC_clusters, features = "MOL3", group.by = 'seurat_clusters')
VlnMOL4 <- VlnPlot(object = NucSeq_ODC_clusters, features = "MOL4", group.by = 'seurat_clusters')
VlnDAO <- VlnPlot(object = NucSeq_ODC_clusters, features = "DAO", group.by = 'seurat_clusters')


# #=====================
cluster_annotations <- list(
    "0" = "DAO",
    "1" = "MOL",
    "2" = "MOL",
    "3" = "Homeostatic MOL",
    "4" = "MOL",
    "5" = "COP/NFOL/MFOL"
)

NucSeq_ODC_clusters$ODC_Status <- unlist(cluster_annotations[NucSeq_ODC_clusters$seurat_clusters])

#
# saveRDS
saveRDS(NucSeq_ODC_clusters, file = paste0(out_dir, 'ODC_COP/', 'NucSeq_ODC_clusters_ODC_Status_Updates_06242023.rds'))
