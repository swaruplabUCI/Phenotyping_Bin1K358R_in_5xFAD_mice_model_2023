###==================================
### NOTE Project BIN1
###==================================
# @rootze
# conda activate scRNAnATAC_R

###==================================
### NOTE Date: 06/21/2023
###==================================
# NOTE
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(viridis)
library(cowplot)
library(ggrepel)

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

## loading data
## PATH

# NOTE loading data
data_dir <- '/Collaborations/MODEL_AD/Bin1_5xFAD_snRNAseq_Nov2022/'
out_dir <- '/Collaborations/MODEL_AD/Bin1_5xFAD_snRNAseq_Nov2022/'
fig_dir <- paste0(out_dir, 'Analysis/Label_transfer_Rosenberg_etal_2018/')

# NOTE Saving without c("Doublets", "Unknown")
NucSeq <- readRDS(paste0(out_dir, 'Data/', 'NucSeq_Updated_0120_2023.rds'))


table(NucSeq@meta.data$Sample.ID)
#

#======================================
# NOTE rosenberg reference
reference_PATH <- "/resources/Public_dataset/rosenberg_2018/data/"

rosenberg <- readRDS(paste0(reference_PATH, 'rosenberg_brain_seurat_processed.rds'))
# An object of class Seurat

head(rosenberg@meta.data)

prosenberg <- DimPlot(rosenberg, reduction = "umap", group.by = "cluster_assignment", label = TRUE, label.size = 3,
    repel = TRUE) + ggtitle("Reference cluster_assignment annotations") #+ NoLegend()

#
pdf(paste0(fig_dir, 'UMAP_cluster_assignment_rosenberg_brain_06212023.pdf'), height=10, width=25)
prosenberg
dev.off()

# NOTE keeping the same features from rosenberg as NucSeq
rosenberg <- rosenberg[rownames(rosenberg)[rownames(rosenberg) %in% rownames(NucSeq)],]
# An object of class Seurat

#======================================
# NOTE label transfer

transfer.anchors_spinal <- FindTransferAnchors(
  reference = rosenberg, # rosenberg not working in the later
  query = NucSeq,
  dims = 1:30,
  reference.reduction = "pca"
)

saveRDS(transfer.anchors_spinal, paste0(fig_dir, 'transfer_anchors_rosenberg_06212023.rds'))


#======================================
# NOTE MapQuery
rosenberg <- RunUMAP(
  rosenberg,
  reduction = 'pca',
  dims = 1:30, n.neighbors=30L, min.dist=0.10, return.model=TRUE
)

#======================================
# NOTE MapQuery
NucSeq <- MapQuery(anchorset = transfer.anchors_spinal,
    reference = rosenberg,
    query = NucSeq,
    refdata = list(celltype = "cluster_assignment"),
    reference.reduction = "pca",
    reduction.model = "umap"
)

# UMAP
p1 <- DimPlot(rosenberg, reduction = "umap", group.by = "cluster_assignment", label = TRUE, label.size = 3,
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(NucSeq, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")

pdf(paste0(fig_dir, 'UMAP_BIN1_cluster_Rosenberg_brain_labeltransfer_reference_query_06212023.pdf'), height=6, width=14)
p1 + p2
dev.off()


p2_only <- DimPlot(NucSeq, reduction = "umap", group.by = "predicted.celltype", label = TRUE,
    label.size = 3, repel = TRUE) #+ NoLegend() #+ ggtitle("Query transferred labels")

pdf(paste0(fig_dir, 'UMAP_BIN1_only_06212023.pdf'), height=12, width=28)
p2_only
dev.off()

table(NucSeq$predicted.celltype)


quantile(NucSeq$predicted.celltype.score)

NucSeq_ref <- NucSeq

table(NucSeq_ref$predicted.celltype)

NucSeq$predictedcelltype <- NucSeq$predicted.celltype
head(NucSeq@meta.data)
table(NucSeq$predicted.celltype)
table(NucSeq$predictedcelltype)

#=================================================
# NOTE make changes to simply cluster names
# NOTE skip here
#=================================================

# NucSeq$predictedcelltype[NucSeq$Celltypes == 'Excitatory_Neuron' & NucSeq$predictedcelltype == 'Unresolved'] <- 'CTX EX'

p <- DimPlot(NucSeq, reduction = "umap", group.by = "predictedcelltype", label = TRUE,
    label.size = 3, repel = TRUE) #+ NoLegend() #+ ggtitle("Query transferred labels")

pdf(paste0(fig_dir, 'UMAP_BIN1_with_newSubclusters_06212023.pdf'), height=12, width=20)
p
dev.off()


p_nolabel <- DimPlot(NucSeq, reduction = "umap", group.by = "predictedcelltype", label = FALSE,
    label.size = 3, repel = TRUE) #+ NoLegend() #+ ggtitle("Query transferred labels")

png(paste0(fig_dir, 'UMAP_BIN1_with_newSubclusters_06212023.png'), height=8, width=16, units='in', res=600)
p_nolabel
dev.off()


# NOTE Saving without c("Doublets", "Unknown")
# NucSeq <- readRDS(paste0(out_dir, 'Data/', 'NucSeq_reprocessed_harmonized_Updated_0120_2023.rds'))
# saveRDS
saveRDS(NucSeq, paste0(out_dir, 'Data/', 'NucSeq_reprocessed_harmonized_Updated_w_rosenberg_predictedcelltype_06212023.rds'))



# ##============================================
# NOTE reloading for plotting the prediction figures
# NucSeq <- readRDS(paste0(out_dir, 'Data/', 'NucSeq_reprocessed_harmonized_Updated_w_rosenberg_predictedcelltype_06212023.rds'))
# NOTE Prediction Plotting
# ##============================================

# NOTE
seurat_obj <- NucSeq
DefaultAssay(seurat_obj) <- 'prediction.score.celltype'

library(RColorBrewer)
# need to load library(RColorBrewer) before Seurat
# Or NOTE library(RColorBrewer)
# Error in value[[3L]](cond) :
#   Package ‘RColorBrewer’ version 1.1.2 cannot be unloaded:
#  Error in unloadNamespace(package) : namespace ‘RColorBrewer’ is imported by ‘Seurat’ so cannot be unloaded
colfunc <- colorRampPalette(c(rev(brewer.pal(9, 'Purples' )[2:9])))

# prediction_matrix <- GetAssayData(seurat_obj, assay='predictions')
prediction_matrix <- GetAssayData(seurat_obj, assay='prediction.score.celltype')

# p1 <- DimPlot(seurat_obj, label=TRUE, group.by='leiden', reduction='umap') + umap_theme + NoLegend()
# p1 <- DimPlot(seurat_obj, label=TRUE, reduction = "ref.umap") + umap_theme + NoLegend()
p1 <- DimPlot(seurat_obj, label=TRUE, reduction = "umap") + umap_theme + NoLegend()


# NOTE I changed the the following code reduction = "ref.umap" to reduction = "umap"
for(label in rownames(seurat_obj)[rowSums(prediction_matrix) > 0]){
  name <- gsub(' ', '_', label)
  name <- gsub('/', '_', label)
  print(name)
  # umap feature plot
  p2 <- FeaturePlot(seurat_obj, features=label, order=TRUE, reduction = "umap") + # reduction = "ref.umap"
    scale_color_gradientn(colors=rev(colfunc(256)), guide = guide_colorbar(barwidth=15, barheight=0.5, ticks=FALSE)) +
    umap_theme + theme(legend.position='bottom')
  # cluster violin plot:
  # p3 <- VlnPlot(seurat_obj, features=label, pt.size=0, group.by='leiden') +
  p3 <- VlnPlot(seurat_obj, features=label, pt.size=0) +
   NoLegend() + ggtitle('') +
   ylab(paste(label, 'score')) + xlab('clusters')
  # patchwork
  patch <- (p1 + p2) / p3
  pdf(paste0(fig_dir, 'Prediction_figures/', name, '.pdf'), width=10, height=10, useDingbats=FALSE)
  print(patch + plot_layout(heights=c(2,1)))
  dev.off()
}
