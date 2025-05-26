library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(ggpubr)
library(AUCell)
library(GSEABase)
################ NKT subset ###################
GC_all_immune <- readRDS("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/GC_ALL_IMMUNE.rds")
table(GC_all_immune$orig.ident)
GC_all_immune@meta.data$"Cell_type" <- plyr::revalue(as.character(GC_all_immune$seurat_clusters),
                                                     c("0" = "CD4_Tcells",
                                                       "1" = "B_cells",
                                                       "2" = "CD8_Tcells",
                                                       "3" = "Macrophages", 
                                                       "4" = "Macrophages/DCs/Monocytes",
                                                       "5" = "NKT",
                                                       "6" = "CD8_Tcells",
                                                       "7" = "Tregs",
                                                       "8" = "CD8_Tcells_cytotoxic",
                                                       "9" = "Cycling_cells",
                                                       "10" = "B_cells",
                                                       "11" = "Melanoma_immune_like",
                                                       "12" = "pDC",
                                                       "13" = "Mast_cells",
                                                       "14" = "B_cells",
                                                       "15" = "Plasma_cells",
                                                       "16" = "B_cells",
                                                       "17" = "CD4_Tcells"
                                                     ))
table(GC_all_immune$Cell_type)
DimPlot(GC_all_immune, label = T)
NKT <-GC_all_immune %>% subset(Cell_type %in% c("CD8_Tcells","CD8_Tcells_cytotoxic", "NKT"))
NKT
table(NKT$orig.ident, NKT$`GC number`)
NKT <- SCTransform(NKT, verbose = FALSE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'))
NKT <- RunPCA(NKT , verbose = TRUE)
NKT <- RunHarmony(NKT, group.by.vars = "orig.ident", assay.use="SCT")
#pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Tcells/Harmony_Heatmap_all.pdf", width = 7, height = 7)
#harmony_embeddings <- Embeddings(Macro, 'harmony')
#harmony_embeddings[1:5, 1:5]
#col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
#Heatmap(harmony_embeddings, 
#        cluster_rows = TRUE, 
#        cluster_columns = FALSE,  
#        clustering_distance_columns = "euclidean",
#        clustering_method_columns = "complete", 
#        show_column_names = TRUE,
#        show_row_names = FALSE,
#        name = "Hramony_embeedding",
        #heatmap_legend_param = list(legend_direction = "vertical", title_position = "leftcenter"),
#        col = col_fun,
#        row_title_rot = 0)
#dev.off()

NKT  <- FindNeighbors(NKT , dims = 1:5, reduction = "harmony") #1:5
NKT  <- FindClusters(NKT , resolution = 0.2, reduction = "harmony") #0.2
NKT  <- RunUMAP(NKT , dims=1:5, reduction = "harmony")
DimPlot(NKT, group.by = c('seurat_clusters'), label=T)
DimPlot(NKT, group.by = "BT/OT")

DimPlot(NKT, group.by = c('Cell_type_CD8_NK'), pt.size = 2, cols = "Dark2") +NoAxes()
dev.copy2pdf(file="/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/CD8_Tcells/figures_for_paper/UMAP_T_NK_cells.pdf",useDingbats=FALSE,family="sans")

################################## Marker genes ##############################
NKT_markers <- FindAllMarkers(NKT, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(NKT_markers, "/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/CD8_Tcells/NKTs_subset_markers.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)
top10 <- NKT_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/CD8_Tcells/NKTs_subset_markers_top10.pdf", width = 15, height = 15)
DoHeatmap(NKT, features = top10$gene) + NoLegend()
dev.off()
Idents(NKT) <- "Cell_type_CD8_NK"
cluster.averages <- AverageExpression(NKT, return.seurat = TRUE)

pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/CD8_Tcells/NKTs_subset_markers_top10_avarage.pdf", width = 10, height = 10)
DoHeatmap(cluster.averages, features = top10$gene) + NoLegend()
dev.off()

colourCount = length(unique(NKT$Cell_type_CD8_NK))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
new_levels <- c("Memory_CD8_Tcells", "Activated_CD8_Tcells", "Dysfunctional_CD8_Tcells","Cytotoxic_CD8_Tcells", "NK")
levels(cluster.averages) <- new_levels
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/CD8_Tcells/NKTs_subset_markers_selected_markers_avarage.pdf", width = 5, height = 5)
DoHeatmap(cluster.averages, features = c("IL7R","FOS","TCF7","SELL","CCL4L2","CCL4","CCL3","IFNG","CXCL13","LAG3","CTLA4","HAVCR2","GZMB","LGALS1","PRF1","GNLY","TYROBP","FCER1G","TRDC","KLRB1","KLRF1","NCAM1", "CD8A", "CD4", "CD3E"), size = 4, 
          draw.lines = FALSE, group.colors  = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"))   + scale_fill_gradientn(colours = rev(mapal)) + NoLegend() + theme(plot.margin=unit(c(2,1.5,1.5,1.5),"cm")) 
dev.off()
pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/CD8_Tcells/NKTs_subset_markers_selected_markers_dotplot.pdf", width = 5, height = 7)
DotPlot(NKT, features = c("CD8A", "CD4", "CD3E","TRAC", "CCL4L2","CCL4","CCL3","IFNG", "GNLY", "GZMB","LGALS1", "PRF1", "CXCL13","LAG3","CTLA4","HAVCR2","PDCD1", "CD274","TIGIT", "IL7R", "FOS","TCF7","SELL","TYROBP","FCER1G","TRDC","KLRB1","KLRF1","NCAM1", "NKG7", "NCR1"), cols = c("blue", "red"), scale = F, group.by = "Cell_type_CD8_NK")  +  RotatedAxis() + coord_flip()
dev.off()

#saveRDS(NKT,"/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/CD8_Tcells/TCells_NKs.rds")

