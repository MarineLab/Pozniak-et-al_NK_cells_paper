library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(ggpubr)
library(AUCell)
library(GSEABase)

Macrophages<-readRDS("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages_Monocytes.rds")
Macro <-Macrophages %>% subset(seurat_clusters %in% c("0"))
DimPlot(Macrophages)
Macro
table(Macro$orig.ident)

#GC_all_immune <- subset(GC_all_immune, orig.ident != "sc5rCMA188")
Macro <- SCTransform(Macro, verbose = FALSE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'))
Macro <- RunPCA(Macro , verbose = TRUE)
Macro <- RunHarmony(Macro, group.by.vars = "orig.ident", assay.use="SCT")

#pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages/Harmony_Heatmap_all.pdf", width = 7, height = 7)
#harmony_embeddings <- Embeddings(Macro, 'harmony')
#harmony_embeddings[1:5, 1:5]
#col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
#ComplexHeatmap::Heatmap(harmony_embeddings, 
#     cluster_rows = TRUE, 
#      cluster_columns = FALSE,  
#       clustering_distance_columns = "euclidean",
#       clustering_method_columns = "complete", 
#       show_column_names = TRUE,
#      show_row_names = FALSE,
 #     name = "Hramony_embeedding",
      #heatmap_legend_param = list(legend_direction = "vertical", title_position = "leftcenter"),
 #    col = col_fun,
 #    row_title_rot = 0)
#dev.off()

Macro  <- FindNeighbors(Macro , dims = 1:5, reduction = "harmony")
Macro  <- FindClusters(Macro , resolution = 0.2, reduction = "harmony")
Macro  <- RunUMAP(Macro , dims=1:5, reduction = "harmony")
DimPlot(Macro, group.by = c('seurat_clusters'))

########################## markers of macrophages subtypes
Macro_markers <- FindAllMarkers(Macro, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
write.table(Macro_markers, "/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages/Macro_subset_markers.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)
top20 <- Macro_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages/Macro_subset_markers_top20.pdf", width = 15, height = 15)
DoHeatmap(Macro, features = top20$gene) + NoLegend()
dev.off()


Macro@meta.data$"Macro_type" <- plyr::revalue(as.character(Macro$seurat_clusters),
                                              c("0" = "Macrophages_CCL3",
                                                "1" = "Macrophages_CXCL9",
                                                "2" = "Macrophages_PLA2G2D",
                                                "3" = "Macrophages_LYVE1",
                                                "4" = "Melanophages"
                                              ))

DimPlot(Macro, group.by = "Macro_type", pt.size = 1)+NoAxes()
dev.copy2pdf(file="/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages/Macrophages_reclustered.pdf", width = 6, height = 4,useDingbats=FALSE,family="sans")

#saveRDS(Macro, "/Users/u0128760/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages/Macro_subtypes.rds")
