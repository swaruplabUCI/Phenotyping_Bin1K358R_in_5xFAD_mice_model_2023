#================================== NOTE ==================================
# Project from Vivek
# NOTE BIN1
# MODLE_AD
# Time received: 2022
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
library(CellChat)
library(patchwork)
#
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
## NOTE Date: 06/27/2023
###=================================
##########################################################################
# NOTE

####################################
###=================================
# NOTE loading PATH
loading_dir <- '/dfs7/swaruplab/zechuas/Collaborations/MODEL_AD/Bin1_5xFAD_snRNAseq_Nov2022/'
subclusters_dir <- '/dfs7/swaruplab/zechuas/Collaborations/MODEL_AD/Bin1_5xFAD_snRNAseq_Nov2022/Analysis/Subcluster_analysis/'
# LT_fig_dir <- paste0(data_dir, 'Analysis/Label_transfer_Rosenberg_etal_2018/')
out_dir <- '/dfs7/swaruplab/zechuas/Collaborations/MODEL_AD/Bin1_5xFAD_snRNAseq_Nov2022/Analysis/CellChat/'

# NOTE setting working directory
setwd(out_dir)

data_dir <- 'data/'
fig_dir <- 'figures/FADvsFAD_BIN1/'

# dir.create(paste0())

#=================
# NOTE loading
NucSeq <- readRDS(paste0(out_dir, 'data/', 'NucSeq_wODC_subclusters_for_CellChat_06262023.rds'))


#########################################################################
# ####=================================
# # NOTE CellChat
# ####=================================
# library(CellChat)
# library(patchwork)

# #=================
# # NOTE change the Idents
table(Idents(NucSeq))
Idents(NucSeq) <- NucSeq$celltype_subcluster

seurat_obj <- NucSeq


data.input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
labels <- Idents(seurat_obj)
meta <- data.frame(group = labels, row.names = names(labels))
conditions <- seurat_obj$Genotype %>% unique %>% as.character

# make a list of cellchat objects:
CellChatDB <- CellChatDB.mouse

cellchat_list <- list()
for(cond in conditions){
  print(cond)
  cellchat_list[[cond]] <- createCellChat(
    object = data.input[,seurat_obj$Genotype == cond],
    meta = seurat_obj@meta.data %>% subset(Genotype == cond),
    group.by = "celltype_subcluster"
  )
  cellchat_list[[cond]]@DB <- CellChatDB
}

################################################################################
# process data
################################################################################

future::plan("multiprocess", workers = 8)

for(cond in conditions){
  print(cond)
  cellchat_list[[cond]] <- subsetData(cellchat_list[[cond]])
  cellchat_list[[cond]] <- identifyOverExpressedGenes(cellchat_list[[cond]])
  cellchat_list[[cond]] <- identifyOverExpressedInteractions(cellchat_list[[cond]])
  cellchat_list[[cond]] <- projectData(cellchat_list[[cond]], PPI.mouse)
  cellchat_list[[cond]] <- computeCommunProb(cellchat_list[[cond]], raw.use = TRUE)
  cellchat_list[[cond]] <- filterCommunication(cellchat_list[[cond]], min.cells = 50)
  df.net <- subsetCommunication(cellchat_list[[cond]])
  cellchat_list[[cond]] <- computeCommunProbPathway(cellchat_list[[cond]])
  cellchat_list[[cond]] <- aggregateNet(cellchat_list[[cond]])
  cellchat_list[[cond]] <- netAnalysis_computeCentrality(cellchat_list[[cond]], slot.name = "netP")

  # save individual cellchat object:
  print(paste('Finishing filtering, processing and saving ', cond))
  saveRDS(cellchat_list[[cond]], file=paste0(data_dir, gsub(' ', '_', cond), '_cellchat_filtered_n_processed_06202023.rds'))
}

#
# merge datasets into one cellchat object:
cellchat <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

saveRDS(cellchat, file = paste0(data_dir, 'cellchat_merged_MODEL_AD_BIN1_06272023.rds'))


# # # #=================
# # NOTE re-load cellchat objects:
# # # #=================
cellchat_list <- list()
conditions <- seurat_obj$Genotype %>% unique %>% as.character
# > conditions
conditions <- conditions[-c(2, 4)]

for(cond in conditions){
  cellchat_list[[cond]] <- readRDS(paste0(data_dir, gsub(' ', '_', cond), '_cellchat_filtered_n_processed_06202023.rds'))
}

# checking the levels
group_new_5xFAD_HEMI = levels(cellchat_list[['5xFAD_HEMI']]@idents)
group_new_5xFAD_HEMI_BIN1HO = levels(cellchat_list[['5xFAD_HEMI-BIN1HO']]@idents)
# NOTE

cellchat <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

# # #=================
# NOTE change for later plotting -- make it easier
object.list <- cellchat_list

gg1 <- netVisual_heatmap(cellchat,
  sources.use = c('Astrocyte', 'Microglia', 'Oligo DAO', 'Oligo Homeostatic MOL', 'Oligo MOL'), remove.isolate = TRUE)
gg2 <- netVisual_heatmap(cellchat,
  sources.use = c('Astrocyte', 'Microglia', 'Oligo DAO', 'Oligo Homeostatic MOL', 'Oligo MOL'), measure = "weight", remove.isolate = TRUE)

pdf(paste0(fig_dir, 'cellchat_5xFAD_vs_5xFAD_BIN1HO_netVisual_heatmap_Ast_ODC_06282023.pdf'), width=14, height=8)
gg1 + gg2
dev.off()



gg1 <- netVisual_heatmap(cellchat,
  sources.use = c('Astrocyte', 'Microglia', 'Oligo DAO', 'Oligo Homeostatic MOL', 'Oligo MOL')) # , remove.isolate = TRUE
gg2 <- netVisual_heatmap(cellchat,
  sources.use = c('Astrocyte', 'Microglia', 'Oligo DAO', 'Oligo Homeostatic MOL', 'Oligo MOL'), measure = "weight") # , remove.isolate = TRUE

pdf(paste0(fig_dir, 'cellchat_5xFAD_vs_5xFAD_BIN1HO_netVisual_heatmap_Ast_ODC_06302023.pdf'), width=14, height=8)
gg1 + gg2
dev.off()

#================================== NOTE ==================================
