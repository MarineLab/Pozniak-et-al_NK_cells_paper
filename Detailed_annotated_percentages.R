library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(ggpubr)
library(AUCell)
library(GSEABase)
library(RColorBrewer)

GC_all_immune <- readRDS("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/GC_ALL_IMMUNE.rds")
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

NKT <- readRDS("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/CD8_Tcells/TCells_NKs.rds")
NKT@meta.data$"Cell_type_CD8_NK" <- plyr::revalue(as.character(NKT$seurat_clusters),
                                                  c("0" = "Dysfuntional_CD8_Tcells",
                                                    "1" = "Memory_CD8_Tcells",
                                                    "2" = "Activated_CD8_Tcells",
                                                    "3" = "Cytotoxic_CD8_Tcells", 
                                                    "4" = "NK"
                                                  ))


GC_all_immune$NEW_cluster <- NKT$Cell_type_CD8_NK
GC_all_immune@meta.data$NEW_cluster <- ifelse(is.na(GC_all_immune$NEW_cluster) ,GC_all_immune$Cell_type, GC_all_immune$NEW_cluster)
DimPlot(GC_all_immune, group.by = "NEW_cluster")

Macrophages <- readRDS("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages_DCs/Macrophages.rds")
GC_all_immune$NEW_cluster_all <- Macrophages$NEW_cluster_Macro_DC
GC_all_immune@meta.data$NEW_cluster_all <- ifelse(is.na(GC_all_immune$NEW_cluster_all) ,GC_all_immune$NEW_cluster, GC_all_immune$NEW_cluster_all)
DimPlot(GC_all_immune, group.by = "NEW_cluster_all")

Macro <- readRDS("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Macrophages/Macro_subtypes.rds")
GC_all_immune$NEW_cluster_all_1 <- Macro$Macro_type
GC_all_immune@meta.data$NEW_cluster_all_1 <- ifelse(is.na(GC_all_immune$NEW_cluster_all_1), GC_all_immune$NEW_cluster_all, GC_all_immune$NEW_cluster_all_1)
DimPlot(GC_all_immune, group.by = "NEW_cluster_all_1")

#saveRDS(GC_all_immune, "/Users/Desktop/GC_ALL_immune_NEW_cluster_all_1.rds")

GC_all_immune <- readRDS("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/GC_ALL_immune_NEW_cluster_all_1.rds")
FeaturePlot(GC_all_immune, features = c("CD3D", "CD8A", "CD4","MS4A1", "NCAM1", "ITGAM", "ITGAX", "CD14", "FCGR3A", "CD68","CD163", "CD63"), ncol = 6)  & Seurat::NoAxes() 
#dev.copy2pdf(file="/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/UMAP_selected_marker.pdf", width = 18, height = 7,useDingbats=FALSE,family="sans")
VlnPlot(GC_all_immune, "mean_cnv", group.by = "NEW_cluster_all_FINAL") +NoLegend()
GC_all_immune@meta.data$"NEW_cluster_all_FINAL" <- plyr::revalue(as.character(GC_all_immune$NEW_cluster_all_1),
                                                                c( "DC1" = "cDC1",
                                                                   "DC2" = "cDC2",
                                                                   "Dysfuntional_CD8_Tcells" = "Dysfunctional_CD8_Tcells",
                                                                   "Melanoma_immune_like" = "Melanocytes",
                                                                   "Cycling_CD8_Tcells" = "Cycling_cells",
                                                                   "Macrophages_necrosis" = "Macrophages_SPP1",
                                                                   "Macrophages_CXCL10" = "Macrophages_CXCL9",
                                                                   "Macrophages_FOS" = "Macrophages_CCL3",
                                                                   "M2 Macrophages" = "Macrophages_LYVE1"
                                                                   
                                                                ))

GC_all_immune$NEW_cluster_all_FINAL <- factor(GC_all_immune$NEW_cluster_all_FINAL, levels = c("Activated_CD8_Tcells",
                                                                                            "Cytotoxic_CD8_Tcells",
                                                                                            "Memory_CD8_Tcells",
                                                                                            "Dysfunctional_CD8_Tcells",
                                                                                            "CD4_Tcells",
                                                                                            "Tregs",
                                                                                            "NK",
                                                                                            "B_cells",
                                                                                            "Plasma_cells",
                                                                                            "Mast_cells",
                                                                                            "Monocytes_TNFSF13",
                                                                                            "Monocytes_CD14",
                                                                                            "Monocytes_CD16",
                                                                                            "Macrophages_CXCL9",
                                                                                            "Macrophages_CCL3",
                                                                                            "Macrophages_SPP1",
                                                                                            "Macrophages_PLA2G2D",
                                                                                            "Macrophages_LYVE1",
                                                                                            "Melanophages",
                                                                                            "Melanocytes",
                                                                                            "cDC1",
                                                                                            "cDC2",
                                                                                            "pDC",
                                                                                            "Cycling_cells",
                                                                                            "Doublets"))

DimPlot(GC_all_immune, group.by = "NEW_cluster_all_FINAL")  + NoAxes()

GC_all_immune@meta.data$"NEW_cluster_all_FINAL_low_res" <- plyr::revalue(as.character(GC_all_immune$seurat_clusters),
                                                                 c( "0" = "CD4_Tcells",
                                                                    "1" = "B_cells",
                                                                    "2" = "CD8_Tcells",
                                                                    "3" = "Myeloid_cells",
                                                                    "4" = "Myeloid_cells",
                                                                    "5" = "CD8_Tcells/NKcells",
                                                                    "6" = "CD8_Tcells",
                                                                    "7" = "Tregs",
                                                                    "8" = "CD8_Tcells",
                                                                    "9" = "Cycling_cells",
                                                                    "10" = "B_cells",
                                                                    "11" = "Melanocytes",
                                                                    "12" = "pDC",
                                                                    "13" = "Mast_cells",
                                                                    "14" = "B_cells",
                                                                    "15" = "Plasma_cells",
                                                                    "16" = "B_cells",
                                                                    "17" = "CD4_Tcells"
                                                                 ))

GC_all_immune$NEW_cluster_all_FINAL_low_res <- factor(GC_all_immune$NEW_cluster_all_FINAL_low_res, levels = c("Myeloid_cells",
                                                                                              "pDC",
                                                                                              "B_cells",
                                                                                              "Plasma_cells",
                                                                                              "CD4_Tcells",
                                                                                              "Tregs",
                                                                                              "CD8_Tcells",
                                                                                              "CD8_Tcells/NKcells",
                                                                                              "Mast_cells",
                                                                                              "Melanocytes",
                                                                                              "Cycling_cells"))


colourCount = length(unique(GC_all_immune$NEW_cluster_all_FINAL_low_res))
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
DimPlot(GC_all_immune, group.by =  "NEW_cluster_all_FINAL_low_res", cols = getPalette(colourCount)) + NoAxes()
#dev.copy2pdf(file="/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/UMAP_annotated.pdf", width = 6, height = 4,useDingbats=FALSE,family="sans")

#colourCount = length(unique(GC_all_immune$NEW_cluster_all_FINAL))
getPalette1 = colorRampPalette(brewer.pal(8,"Set1"))(15) #6
#"#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33"
getPalette2 = colorRampPalette(brewer.pal(12, "Set3"))(10)
rampcols <- c(getPalette1, getPalette2)
rampcols
DimPlot(GC_all_immune, group.by =  "NEW_cluster_all_FINAL", cols = rampcols) + NoAxes() +   guides(color = guide_legend(override.aes = list(size=4), ncol=1) )
#dev.copy2pdf(file="/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/UMAP_annotated_detailed.pdf", width = 9, height = 6,useDingbats=FALSE,family="sans")

################# add updated response info
library(readxl)
Meta_data <- read_excel("/Users/Documents/PROJECTS/Grand_Challenge/Meta_data/Meta_data_GEX_RESPONSE_2.0-25-08-2021.xlsx")
Meta_data <- as.data.frame(Meta_data)
GC_all_immune@meta.data$row_names <- rownames(GC_all_immune@meta.data)
dim(GC_all_immune@meta.data)
GC_all_immune@meta.data<-GC_all_immune@meta.data %>% left_join(Meta_data, by="orig.ident")
row.names(GC_all_immune@meta.data) <- GC_all_immune@meta.data$row_names
GC_all_immune@meta.data$"Response" <- plyr::revalue(as.character(GC_all_immune$Response),
                                                    c( "0" = "NR",
                                                       "1" = "R"
                                                    ))
#pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/UMAP_annotated_response.pdf", width = 8, height = 5)
#DimPlot(GC_all_immune, group.by =  "Response")+NoAxes()
#dev.off()
#pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/UMAP_annotated_patient.pdf", width = 10, height = 5)
#DimPlot(GC_all_immune, group.by =  "orig.ident")+NoAxes()
#dev.off()

############################### Marker genes + pathway enrichemnts
####################################### Cluster Identity #########################

Idents(GC_all_immune) <- "NEW_cluster_all_FINAL"
GC_all_immune_markers_cell_type <- FindAllMarkers(GC_all_immune, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
#write.table(GC_all_immune_markers_cell_type, "/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/GC_all_immune_markers_cell_type.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)
top5 <- GC_all_immune_markers_cell_type %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
GC_all_immune_sub <- GC_all_immune %>% subset(NEW_cluster_all_FINAL %not_in% c("Melanocytes","Doublets"))

pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Heatmap_markers_immune_top10.pdf", width = 10, height = 12)
DoHeatmap(GC_all_immune_sub, features = top5$gene,  group.colors  = rampcols) + NoLegend() + theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
dev.off()

cluster.averages <- AverageExpression(GC_all_immune_sub, return.seurat = TRUE)

#new_levels <- c("B_cells", "CD4_Tcells", "CD8_Tcells", "CD8_Tcells_cytotoxic", "Cycling_cells", "Macrophages", "Macrophages/DCs/Monocytes", "Mast_cells", "Melanoma_Immune", "pDC", "Plasma_cells", "Tcell/NK", "Tregs")

#levels(cluster.averages) <- new_levels
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Heatmap_markers_immune_top10_average.pdf", width = 9, height = 14)
DoHeatmap(cluster.averages, features = top5$gene, size = 4, 
          draw.lines = FALSE, group.colors  = rampcols) + scale_fill_gradientn(colours = rev(mapal)) + NoLegend() + theme(plot.margin=unit(c(1.5,1.5,1.5,1.5),"cm")) 
dev.off()

table(GC_all_immune$NEW_cluster_all_FINAL,GC_all_immune$orig.ident)

############ Calculate percentages
GC_all_immune@meta.data$"Treatment" <- plyr::revalue(as.character(GC_all_immune$`Treatment`),
                                           c("Nivolumab + Ipilimumab" = "NNivolumabIpilimumab"
                                           ))
GC_all_immune@meta.data$"TILs" <- plyr::revalue(as.character(GC_all_immune$orig.ident),
                                                c("scrCMA036"="Brisk",#was NA
                                                  "scrCMA041" = "NonBrisk",#was NA
                                                  "scrCMA046" = "Brisk", #was NA
                                                  "scrCMA038" = "NonBrisk",
                                                  "scrCMA044" = "Brisk",
                                                  "scrCMA040" = "Absent",
                                                  "scrCMA048" = "NonBrisk",
                                                  "scrCMA049" = "NonBrisk",
                                                  "scrCMA050" = "NonBrisk",
                                                  "scrCMA054" = "Absent",
                                                  "scrCMA063" = "Absent",
                                                  "scrCMA055" = "Absent",
                                                  "scrCMA064" = "NonBrisk",
                                                  "sc5rCMA061" = "Absent", #was NA
                                                  "sc5rCMA066" = "Absent",
                                                  "scrCMA068" = "NonBrisk",
                                                  "scrCMA072" = "NonBrisk",
                                                  "sc5rCMA070" = "NonBrisk",
                                                  "sc5rCMA074" = "NA",
                                                  "scrCMA076" = "NonBrisk",
                                                  "scrCMA088" = "Brisk",
                                                  "scrCMA077" = "NonBrisk",
                                                  "scrCMA087" = "NA",
                                                  "scrCMA089" = "Absent",#was NA
                                                  "scrCMA090" = "Absent",
                                                  "scrCMA091" = "NonBrisk",
                                                  "scrCMA093" = "Necrosis",
                                                  "scrCMA112" = "Necrosis",
                                                  "scrCMA094" = "Brisk", #was NA
                                                  "scrCMA109" = "NoTumor",
                                                  "scrCMA119" = "Brisk",
                                                  "scrCMA129" = "Necrosis",
                                                  "scrCMA120" = "Absent",
                                                  "scrCMA121" = "NonBrisk",
                                                  "scrCMA130" = "Absent",#was NA
                                                  "scrCMA131" = "Brisk",
                                                  "sc5rCMA136" = "NonBrisk",
                                                  "sc5rCMA149" = "NonBrisk",
                                                  "sc5rCMA141" = "NonBrisk",
                                                  "sc5rCMA152" = "NA",
                                                  "sc5rCMA144" = "Brisk",
                                                  "sc5rCMA155" = "NA" ))
GC_all_immune@meta.data$"Tissue" <- plyr::revalue(as.character(GC_all_immune$`Tissue`),
                                                  c("Lymph node" = "LN"
                                                  ))
GC_all_immune@meta.data$"Mutation" <- plyr::revalue(as.character(GC_all_immune$`Mut type`),
                                                    c("N/A" = "NA"
                                                    ))

GC_all_immune@meta.data$"TILs_bin" <- plyr::revalue(as.character(GC_all_immune$TILs),
                                                    c("Absent" = "Absent",
                                                      "NonBrisk" = "Present",
                                                      "Brisk" = "Present"
                                                    ))


############################
#####################################
##############################################
###########################################################################
`%not_in%` <- purrr::negate(`%in%`)
GC_all_immune_sub <- GC_all_immune %>% subset(NEW_cluster_all_FINAL %not_in% c("Melanocytes","Doublets"))
table(GC_all_immune_sub$NEW_cluster_all_FINAL)
TEST <- GC_all_immune_sub@meta.data
cell_num <- TEST %>%
  mutate(sample_id = as.factor(paste(orig.ident, `GC number`, `TILs_bin`,`TILs`, `Tissue`,`Mutation`,`Response`,`BT/OT`,`Treatment`, sep="_"))) %>%
  mutate(NEW_cluster_all_FINAL = as.factor(NEW_cluster_all_FINAL)) %>%
  group_by(sample_id, NEW_cluster_all_FINAL, .drop=FALSE) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::separate(sample_id, c("orig.ident", "GC_number", "TILs_bin", "TILs", "Tissue", "Mutation", "Response", "Timepoint", "Treatment"))
cell_num
TEST1 <- GC_all_immune_sub@meta.data
total_cells<- TEST1 %>%
  group_by(orig.ident) %>%
  dplyr::summarise(total = n())
total_cells
cell_percentage<- left_join(cell_num, total_cells) %>%
  mutate(percentage = n/total*100)
cell_percentage

cell_percentage <- cell_percentage %>% subset(NEW_cluster_all_FINAL %not_in% c("Melanocytes","Doublets")) 

#Response
ggboxplot(cell_percentage, x = "Response", y = "percentage",
          color = "Response", palette =c("#00AFBB", "#FC4E07"),
          shape = "Response",
          add = "jitter")+
  stat_compare_means(aes(group = Response), label = "p.format", method =  "wilcox.test") +
  facet_wrap(~NEW_cluster_all_FINAL, ncol = 13) + theme_classic() + RotatedAxis()  +  theme(text=element_text(size=20))




########################################## TILS #################################
cell_percentage_tils <- cell_percentage %>% subset(TILs %in% c("Brisk", "Absent", "NonBrisk"))
#cell_percentage_tils <- cell_percentage_tils %>% subset(cell_percentage_tils %not_in% c("Cycling_cells", "Doublets"))
pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Cell_percentages_TILs.pdf", width = 25, height = 10)
my_comp <- list(c("Brisk", "Absent"), c("Brisk", "NonBrisk"), c("NonBrisk", "Absent"))
ggboxplot(cell_percentage_tils, x = "TILs", y = "percentage",
          fill = "TILs",
          shape = "TILs",
          add = "jitter")+
  stat_compare_means(comparisons =  my_comp, label = "p.format", method =  "wilcox.test") + NoLegend()+
  facet_wrap(~NEW_cluster_all_FINAL, ncol = 11, scales = "free") + theme_classic() + RotatedAxis() +  theme(text=element_text(size=12),axis.text.x = element_text(angle = 45))
dev.off()

cell_percentage_tils_main <- cell_percentage_tils %>% subset(NEW_cluster_all_FINAL %in% c("Activated_CD8_Tcells",
                                                                                    "Cytotoxic_CD8_Tcells",
                                                                                    "Memory_CD8_Tcells",
                                                                                    "Dysfunctional_CD8_Tcells",
                                                                                    "NK"))

pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Cell_percentages_annotated_detailed_TILS_main.pdf",  width = 18, height = 7)
cell_percentage_tils_main$var <- paste(cell_percentage_tils_main$Timepoint,cell_percentage_tils_main$Response, sep = "_" )
table(cell_percentage_tils_main$var)
my_comp <- list(c("Brisk", "Absent"), c("Brisk", "NonBrisk"), c("NonBrisk", "Absent"))
ggboxplot(cell_percentage_tils_main, x = "TILs", y = "percentage",
          fill = "TILs",
          shape = "TILs",
          palette = c("#3366FF","#996666",  "#CC3333"), 
          add = "jitter")+
  stat_compare_means(comparisons =  my_comp, label = "p.format", method =  "wilcox.test") + NoLegend()+
  facet_wrap(~NEW_cluster_all_FINAL, ncol = 12, scales = "free") + theme_classic() + RotatedAxis() +  theme(text=element_text(size=20),axis.text.x = element_text(angle = 45))
dev.off()


cell_percentage_tils_main <- cell_percentage_tils %>% subset(NEW_cluster_all_FINAL %in% c("NK"))

pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Cell_percentages_annotated_detailed_TILS_main_NK_Response.pdf",  width = 7, height = 7)
cell_percentage_tils_main$var <- paste(cell_percentage_tils_main$Timepoint,cell_percentage_tils_main$Response, sep = "_" )
table(cell_percentage_tils_main$var)
my_comp <- list(c("Brisk", "Absent"), c("Brisk", "NonBrisk"), c("NonBrisk", "Absent"))
#my_comp1 <- list(c("Brisk", "NonBrisk"))
ggboxplot(cell_percentage_tils_main, x = "TILs", y = "percentage",
          fill = "TILs",
          shape = "TILs",
          palette = c("#3366FF","#996666",  "#CC3333"), 
          add = "jitter")+
  stat_compare_means(comparisons =  my_comp, label = "p.format", method =  "wilcox.test") + NoLegend()+
  facet_wrap(~Response, ncol = 12, scales = "free") + theme_classic() + RotatedAxis() +  theme(text=element_text(size=20),axis.text.x = element_text(angle = 45))
dev.off()

############### CD8 T cells among TILs and response
cell_percentage_tils_main <- cell_percentage_tils %>% subset(NEW_cluster_all_FINAL %in% c("Activated_CD8_Tcells",
                                                                                          "Cytotoxic_CD8_Tcells",
                                                                                          "Memory_CD8_Tcells",
                                                                                          "Dysfunctional_CD8_Tcells"))

table(cell_percentage_tils_main$var)
my_comp <- list(c("Brisk", "Absent"), c("Brisk", "NonBrisk"), c("NonBrisk", "Absent"))
#my_comp1 <- list(c("Brisk", "NonBrisk"))

ggboxplot(cell_percentage_tils_main, x = "TILs", y = "percentage",
          fill = "TILs",
          shape = "TILs",
          palette = c("#3366FF","#996666",  "#CC3333"), 
          add = "jitter")+
  stat_compare_means(comparisons =  my_comp, label = "p.format", method =  "wilcox.test") + NoLegend()+
  facet_wrap(~Response+NEW_cluster_all_FINAL, ncol = 4, scales = "free") + theme_classic() + RotatedAxis() +  theme(text=element_text(size=20),axis.text.x = element_text(angle = 45))



cell_percentage_tils_sup <- cell_percentage_tils %>% subset(NEW_cluster_all_FINAL %not_in% c("Activated_CD8_Tcells",
                                                                                        "Cytotoxic_CD8_Tcells",
                                                                                        "Memory_CD8_Tcells",
                                                                                        "Dysfunctional_CD8_Tcells",
                                                                                        "NK", "NA", "Cycling_cells"))

pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Cell_percentages_annotated_detailed_TILS_supplement.pdf",  width = 22, height = 20)
cell_percentage_tils_sup$var <- paste(cell_percentage_tils_sup$Timepoint,cell_percentage_tils_sup$Response, sep = "_" )
table(cell_percentage_tils_sup$var)
my_comp <- list(c("Brisk", "Absent"), c("Brisk", "NonBrisk"), c("NonBrisk", "Absent"))
ggboxplot(cell_percentage_tils_sup, x = "TILs", y = "percentage",
          fill = "TILs",
          shape = "TILs",
          palette = c("#3366FF","#996666",  "#CC3333"), 
          add = "jitter")+
  stat_compare_means(comparisons =  my_comp, label = "p.format", method =  "wilcox.test") + NoLegend()+
  facet_wrap(~NEW_cluster_all_FINAL, ncol = 6, scales = "free") + theme_classic() + RotatedAxis() +  theme(text=element_text(size=20),axis.text.x = element_text(angle = 45))
dev.off()


########## CD8 Tcells TILS only in NR
my_comp <- list(c("Brisk", "Absent"), c("Brisk", "NonBrisk"), c("NonBrisk", "Absent"))

ggboxplot(cell_percentage_tils, x = "TILs", y = "percentage",
          fill = "TILs",
          shape = "TILs",
          add = "jitter")+
  stat_compare_means(comparisons =  my_comp, label = "p.format", method =  "wilcox.test") + NoLegend()+
  facet_wrap(~Response, ncol = 11, scales = "free") + theme_classic() + RotatedAxis() +  theme(text=element_text(size=12),axis.text.x = element_text(angle = 45))
####################

pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Cell_percentages_Tissue.pdf", width = 19, height = 13)
my_comp2 <- list(c("Skin", "LN"), c("Subcutis", "LN"), c("Subcutis", "Skin"))

ggboxplot(cell_percentage, x = "Tissue", y = "percentage",
          fill = "Tissue",
          shape = "Tissue",
          add = "jitter")+
  stat_compare_means(comparisons =  my_comp2, label = "p.format", method =  "wilcox.test") + NoLegend()+
  facet_wrap(~NEW_cluster_all_FINAL, ncol = 8, scales = "free") + theme_classic() + RotatedAxis() +  theme(text=element_text(size=12),axis.text.x = element_text(angle = 45))
dev.off()


pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Cell_percentages_Mutation.pdf", width = 18, height = 7)
ggboxplot(cell_percentage, x = "Mutation", y = "percentage",
          color = "Mutation",
          shape = "Mutation",
          add = "jitter")+
  stat_compare_means(aes(group = Mutation), label = "p.format", method =  "kruskal.test") +
  facet_wrap(~NEW_cluster_all_FINAL, ncol = 13) + theme_classic() + RotatedAxis() +  theme(text=element_text(size=20),axis.text.x = element_text(angle = 45))
dev.off()

pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Stack_per_patient.pdf", width = 19, height = 6)
ggplot(cell_percentage, aes(orig.ident, y=percentage, fill=NEW_cluster_all_FINAL))+geom_bar(stat="identity", position = "stack") +  scale_fill_manual(values = rampcols) +  facet_grid(~GC_number+Timepoint, scale="free", drop = TRUE)+  theme(text=element_text(size=40)) + theme_classic()+RotatedAxis()
dev.off()

pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Stack_per_timepoint_site.pdf", width = 19, height = 6)
ggplot(cell_percentage, aes(orig.ident, y=percentage, fill=NEW_cluster_all_FINAL))+geom_bar(stat="identity", position = "stack") +  scale_fill_manual(values = rampcols) +  facet_grid(~Tissue+Timepoint, scale="free", drop = TRUE)+  theme(text=element_text(size=40)) + theme_classic()+RotatedAxis()
dev.off()



cell_percentage$Mutation <-  factor(cell_percentage$Mutation, levels = c("BRAF", "NRAS", "WT", "NA"))
pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Stack_per_Mutation.pdf", width = 13, height = 6)
ggplot(cell_percentage, aes(orig.ident, y=percentage, fill=NEW_cluster_all_FINAL))+geom_bar(stat="identity", position = "stack", color="white") + facet_grid(~Mutation, scale="free", drop = TRUE) +  theme(text=element_text(size=40)) + theme_classic()+RotatedAxis()
dev.off()

pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Stack_per_Tissue.pdf", width = 19, height = 6)
ggplot(cell_percentage, aes(orig.ident, y=percentage, fill=NEW_cluster_all_FINAL))+geom_bar(stat="identity", position = "stack")   +  scale_fill_manual(values = getPalette(colourCount)) +  facet_grid(~Tissue+Response, scale="free", drop = TRUE) +  theme(text=element_text(size=40)) + theme_classic()+RotatedAxis()
dev.off()


pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Cell_percentages_annotated_detailed_Response.pdf",  width = 20, height = 20)
ggboxplot(cell_percentage, x = "Response", y = "percentage",
          fill = "Response",
          shape = "Response",
          add = "jitter") + 
  stat_compare_means(aes(group = Response), label = "p.format", method =  "wilcox.test") +
  facet_wrap(~NEW_cluster_all_FINAL, ncol = 6, scales = "free") + theme_classic() + RotatedAxis() +  theme(text=element_text(size=20),axis.text.x = element_text(angle = 45)) +NoLegend()
dev.off()

pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Cell_percentages_annotated_detailed_Timepoint.pdf",  width = 20, height = 20)
ggboxplot(cell_percentage, x = "Timepoint", y = "percentage",
          fill = "Timepoint",
          shape = "Timepoint",
          add = "jitter") + 
  stat_compare_means(aes(group = Timepoint), label = "p.format", method =  "wilcox.test") +
  facet_wrap(~NEW_cluster_all_FINAL, ncol = 6, scales = "free") + theme_classic() + RotatedAxis() +  theme(text=element_text(size=13),axis.text.x = element_text(angle = 45)) +NoLegend()
dev.off()


pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Cell_percentages_annotated_detailed_Response_timepoint.pdf",  width = 25, height = 20)
cell_percentage$var <- paste(cell_percentage$Timepoint,cell_percentage$Response, sep = "_" )
table(cell_percentage$var)
my_comp <- list(c("OT_NR", "OT_R"), c("BT_NR", "BT_R"), c("BT_R", "OT_R"), c("BT_NR", "OT_NR"))
ggboxplot(cell_percentage, x = "var", y = "percentage",
          fill = "Timepoint",
          shape = "Response",
          palette = c("#377eb8","#4daf4a"), 
          add = "jitter")+
  stat_compare_means(comparisons =  my_comp, label = "p.format", method =  "wilcox.test") +
  facet_wrap(~NEW_cluster_all_FINAL,ncol = 6, scales = 'free') + theme_classic() + RotatedAxis()  +  theme(text=element_text(size=20))
dev.off()


cell_percentage_main_fig <- cell_percentage %>% subset(NEW_cluster_all_FINAL %in% c("Activated_CD8_Tcells",
                                                                                    "Cytotoxic_CD8_Tcells",
                                                                                    "Memory_CD8_Tcells",
                                                                                    "Dysfunctional_CD8_Tcells",
                                                                                    "NK"))

pdf("/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Cell_percentages_annotated_detailed_Response_timepoint_main_new.pdf",  width = 20, height = 7)
cell_percentage_main_fig$var <- paste(cell_percentage_main_fig$Timepoint,cell_percentage_main_fig$Response, sep = "_" )
table(cell_percentage_main_fig$var)
my_comp <- list(c("OT_NR", "OT_R"), c("BT_NR", "BT_R"), c("BT_R", "OT_R"), c("BT_NR", "OT_NR"))
cell_percentage_main_fig$var <- factor(cell_percentage_main_fig$var, 
                                       levels = c("BT_NR", "OT_NR", "BT_R", "OT_R"))
ggplot(cell_percentage_main_fig, aes(x = var, y = percentage, fill = Timepoint)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group = GC_number), color = "gray", linetype = "dashed") +
  geom_jitter(aes(shape = Response), 
              size = 3, # Increase point size
              width = 0.2, height = 0) +
  stat_compare_means(comparisons = my_comp, label = "p.format", method = "wilcox.test") +
  facet_wrap(~NEW_cluster_all_FINAL, ncol = 6, scales = 'free') +
  scale_fill_manual(values = c("#377eb8", "#4daf4a")) +
  theme_classic() +
  RotatedAxis() +
  labs(y = "log10(percentage)", fill = "Timepoint") 
dev.off()


cell_percentage_main_sup <- cell_percentage %>% subset(NEW_cluster_all_FINAL %not_in% c("Activated_CD8_Tcells",
                                                                                        "Cytotoxic_CD8_Tcells",
                                                                                        "Memory_CD8_Tcells",
                                                                                        "Dysfunctional_CD8_Tcells",
                                                                                        "NK"))

pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Cell_percentages_annotated_detailed_Response_timepoint_supplement_new.pdf",  width = 20, height = 17)
cell_percentage_main_sup <- cell_percentage_main_sup %>% subset(NEW_cluster_all_FINAL %not_in% c("Cycling_cells"))
cell_percentage_main_sup$var <- paste(cell_percentage_main_sup$Timepoint,cell_percentage_main_sup$Response, sep = "_" )
table(cell_percentage_main_sup$var)
my_comp <- list(c("OT_NR", "OT_R"), c("BT_NR", "BT_R"), c("BT_R", "OT_R"), c("BT_NR", "OT_NR"))
cell_percentage_main_sup$var <- factor(cell_percentage_main_sup$var, 
                                       levels = c("BT_NR", "OT_NR", "BT_R", "OT_R"))
ggplot(cell_percentage_main_sup, aes(x = var, y = percentage, fill = Timepoint)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group = GC_number), color = "gray", linetype = "dashed") +
  geom_jitter(aes(shape = Response), 
              size = 3, # Increase point size
              width = 0.2, height = 0) +
  stat_compare_means(comparisons = my_comp, label = "p.format", method = "wilcox.test") +
  facet_wrap(~NEW_cluster_all_FINAL, ncol = 6, scales = 'free') +
  scale_fill_manual(values = c("#377eb8", "#4daf4a")) +
  theme_classic() +
  RotatedAxis() +
  labs(y = "log10(percentage)", fill = "Timepoint") 
dev.off()



###################### PAIRED TEST#######################
library(plyr)
table(cell_percentage$Timepoint,  cell_percentage$GC_number)
`%not_in%` <- purrr::negate(`%in%`)
subset_for_matched <- cell_percentage %>% subset(GC_number %not_in% c("25",
                                                                      "29"))
###drop duplicates samples
subset_for_matched <- subset_for_matched %>% subset(orig.ident %not_in% c("scrCMA048",
                                                                       "scrCMA050",
                                                                     "scrCMA041"))
table(subset_for_matched$Timepoint, subset_for_matched$GC_number)
pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Cell_percentages_NR_paired_supplement.pdf",  width = 25, height = 10)
cell_percentage_NRes <- subset(subset_for_matched, Response == "NR")
cell_percentage_NRes <- cell_percentage_NRes %>% subset(NEW_cluster_all_FINAL %not_in% c("Cycling_cells", "Doublets"))

ggpaired(cell_percentage_NRes, x = "Timepoint", y = "percentage",
         fill = "Timepoint", palette =c("#00AFBB", "#FC4E07"),
         # shape = "Timepoint",
         #add = "jitter", 
         id = "GC_number",
         line.color = "gray", line.size = 0.4)+
  stat_compare_means(label = "p.format", method = "wilcox.test", paired = T) + NoLegend() +
  facet_wrap(~NEW_cluster_all_FINAL, scales = "free", ncol = 12) + theme_classic() + RotatedAxis()  +  theme(text=element_text(size=12))
dev.off()

pdf("/Users/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/All_immune_cells_annotated_deteiled/Cell_percentages_R_paired_supplement.pdf",  width = 25, height = 10)
cell_percentage_Res <- subset(subset_for_matched, Response == "R")
cell_percentage_Res <- cell_percentage_Res %>% subset(NEW_cluster_all_FINAL %not_in% c("Cycling_cells", "Doublets"))

ggpaired(cell_percentage_Res, x = "Timepoint", y = "percentage",
         fill = "Timepoint", palette =c("#00AFBB", "#FC4E07"),
         # shape = "Timepoint",
         #add = "jitter",
         id = "GC_number",
         line.color = "gray", line.size = 0.4)+
  stat_compare_means(label = "p.format", method ="wilcox.test", paired = T) + # Add significance levels + 
  facet_wrap(~NEW_cluster_all_FINAL, scales = "free", ncol = 12) + theme_classic() + RotatedAxis()  +  theme(text=element_text(size=12))
dev.off()


