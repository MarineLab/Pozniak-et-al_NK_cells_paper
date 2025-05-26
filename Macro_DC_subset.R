library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(ggpubr)
library(AUCell)
library(GSEABase)
################ Macro subset ###################
GC_all_immune <- readRDS("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/GC_ALL_IMMUNE.rds")
#GC_all_immune@meta.data$"Tissue" <- plyr::revalue(as.character(GC_all_immune$`Tissue`),
 #                                                    c("Lymph node" = "LN"
 #                                                    ))

DimPlot(GC_all_immune, label = T)
Macrophages <-GC_all_immune %>% subset(seurat_clusters %in% c("3","4", "12"))
Macrophages
table(Macrophages$orig.ident)

Macrophages <- SCTransform(Macrophages, verbose = FALSE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'))
Macrophages <- RunPCA(Macrophages , verbose = TRUE)
Macrophages <- RunHarmony(Macrophages, group.by.vars = "orig.ident", assay.use="SCT")
#pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages_DCs/Harmony_Heatmap_all.pdf", width = 7, height = 7)
#harmony_embeddings <- Embeddings(Macrophages, 'harmony')
#harmony_embeddings[1:5, 1:5]
#col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
#ComplexHeatmap::Heatmap(harmony_embeddings, 
#        cluster_rows = TRUE, 
#        cluster_columns = FALSE,  
#        clustering_distance_columns = "euclidean",
#        clustering_method_columns = "complete", 
#        show_column_names = TRUE,
#        show_row_names = FALSE,
 #       name = "Hramony_embeedding",
        #heatmap_legend_param = list(legend_direction = "vertical", title_position = "leftcenter"),
 #       col = col_fun,
 #       row_title_rot = 0)
#dev.off()

Macrophages  <- FindNeighbors(Macrophages , dims = 1:5, reduction = "harmony")
Macrophages  <- FindClusters(Macrophages , resolution = 0.15, reduction = "harmony")
Macrophages  <- RunUMAP(Macrophages , dims=1:5, reduction = "harmony")

Macrophages@meta.data$"Cell_type_macro_dc" <- plyr::revalue(as.character(Macrophages$seurat_clusters),
                                                            c("0" = "Macrophages",
                                                              "1" = "Monocytes",
                                                              "2" = "DCs",
                                                              "3" = "Macrophages_necrosis", #melano
                                                              "4" = "pDC"
                                                            ))
DimPlot(Macrophages, pt.size = 1, group.by = "Cell_type_macro_dc")+ NoAxes()
dev.copy2pdf(file="/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages_DCs/Myeloid_reclustered.pdf", width = 6, height = 4,useDingbats=FALSE,family="sans")


Macrophages_markers <- FindAllMarkers(Macrophages, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
write.table(Macrophages_markers, "/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages_DCs_subset_markers.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)
top10 <- Macrophages_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages_DCs_subset_markers_top10.pdf", width = 15, height = 15)
DoHeatmap(Macrophages, features = top10$gene) + NoLegend()
dev.off()


#saveRDS(Macrophages,"/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages_Monocytes.rds")


############################## DC subset ##########################################################
###################################################################################################
##################################################################################################
Macrophages <- readRDS("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages_CYTOTRACE/Macrophages_Monocytes.rds")
DimPlot(Macrophages)
DCs <-Macrophages %>% subset(seurat_clusters %in% c("2", "1"))
DCs <- subset(DCs, subset = CD3D < 0.0001)
DimPlot(DCs)
table(DCs$orig.ident)

DCs <- SCTransform(DCs, verbose = FALSE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'))
DCs <- RunPCA(DCs , verbose = TRUE)
DCs <- RunHarmony(DCs, group.by.vars = "orig.ident", assay.use="SCT")
pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages_DCs/DCs_Harmony_Heatmap_all.pdf", width = 7, height = 7)
harmony_embeddings <- Embeddings(DCs, 'harmony')
harmony_embeddings[1:5, 1:5]
col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
ComplexHeatmap::Heatmap(harmony_embeddings, 
                        cluster_rows = TRUE, 
                        cluster_columns = FALSE,  
                        clustering_distance_columns = "euclidean",
                      clustering_method_columns = "complete", 
                        show_column_names = TRUE,
                       show_row_names = FALSE,
                       name = "Hramony_embeedding",
                       #heatmap_legend_param = list(legend_direction = "vertical", title_position = "leftcenter"),
                        col = col_fun,
                        row_title_rot = 0)
dev.off()

DCs  <- FindNeighbors(DCs , dims = 1:6, reduction = "harmony")
DCs  <- FindClusters(DCs , resolution = 0.15, reduction = "harmony")
DCs  <- RunUMAP(DCs , dims=1:6, reduction = "harmony")
DimPlot(DCs, group.by = c('seurat_clusters'))
#saveRDS(DCs, "/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages_DCs/DC.rds")
DCs <- readRDS( "/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages_DCs/DC.rds")
DCs@meta.data$"Cell_type_DC" <- plyr::revalue(as.character(DCs$seurat_clusters),
                                              c("0" = "Monocytes_CD14",
                                                "1" = "cDC2",
                                                "2" = "Monocytes_CD16",
                                                "3" = "Monocytes_TNFSF13",
                                                "4" = "cDC1"
                                              ))
DimPlot(DCs, pt.size = 1, group.by = "Cell_type_DC")+ NoAxes()
dev.copy2pdf(file="/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages_DCs/DC_reclustered.pdf", width = 6, height = 4,useDingbats=FALSE,family="sans")

DCs_markers <- FindAllMarkers(DCs, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
write.table(DCs_markers, "/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages_DCs/DCs_subset_markers.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)
top20 <- DCs_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages_DCs/DCs_subset_markers_top20.pdf", width = 15, height = 15)
DoHeatmap(DCs, features = top20$gene) + NoLegend()
dev.off()

