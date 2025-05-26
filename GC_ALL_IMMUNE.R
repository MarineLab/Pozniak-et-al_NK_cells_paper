library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(ggpubr)
library(AUCell)
library(GSEABase)

GC_ALL_harmony_singlet <- readRDS("/Users//Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/GC_ALL_harmony_singlet_meta_anno.rds")
GC_all_immune <- subset(GC_ALL_harmony_singlet, subset = Immune > 0.1 & mean_cnv < 0.1)
GC_all_immune


table(GC_all_immune$orig.ident, GC_all_immune$`GC number`)
GC_all_immune <- subset(GC_all_immune, orig.ident != "sc5rCMA188") #sample with less than 10 cells
table(GC_all_immune$orig.ident)
GC_all_immune
GC_all_immune <- SCTransform(GC_all_immune, verbose = FALSE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'))
GC_all_immune <- RunPCA(GC_all_immune , verbose = TRUE)
GC_all_immune <- RunHarmony(GC_all_immune, group.by.vars = "orig.ident", assay.use="SCT")

options(future.globals.maxSize = 4000 * 1024^2)
#pdf("/Users//Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Harmony_Heatmap_Immune.pdf", width = 7, height = 7)
#harmony_embeddings <- Embeddings(GC_all_immune, 'harmony')
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
#        #heatmap_legend_param = list(legend_direction = "vertical", title_position = "leftcenter"),
#        col = col_fun,
#        row_title_rot = 0)
#dev.off()

GC_all_immune  <- FindNeighbors(GC_all_immune , dims = 1:18, reduction = "harmony")
GC_all_immune  <- FindClusters(GC_all_immune , resolution = 0.55, reduction = "harmony")
GC_all_immune  <- RunUMAP(GC_all_immune , dims=1:18, reduction = "harmony")
#saveRDS(GC_all_immune, "/Users//Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/GC_ALL_IMMUNE.rds")
#pdf("/Users//Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/UMAP_VLN_IMUNNE.pdf", width = 7, height = 7)
DimPlot(GC_all_immune, group.by = c('seurat_clusters'), label=T)
DimPlot(GC_all_immune, group.by = c("Response")) +theme(legend.text = element_text(size = 5))
DimPlot(GC_all_immune, group.by = c("BT/OT")) +theme(legend.text = element_text(size = 5))

pdf("/Markers_UMAP.pdf", width = 20, height = 20)
FeaturePlot(GC_all_immune, features = c("CD3D", "CD8A", "CD4","FOXP3", "GZMB","MS4A1", "NCAM1", "ITGAM", "ITGAX", "CD14", "FCGR3A", "CD68","CD163", "MLANA"), ncol = 7)  & Seurat::NoAxes() 
dev.off()

pdf("/Users//Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/nCount_Immune_per_sample.pdf", width = 14, height = 7)
give.nmedian <- function(x){
  return(c(y = median(x)*1.1, label = length(x))) 
}
VlnPlot(object = GC_all_immune, features = c("nCount_RNA"), group.by = "orig.ident", pt.size = 0) + 
  #geom_hline(yintercept=20, linetype='dashed') +
  stat_summary(fun.data = give.nmedian, geom = "text", fun = median, size = 3) +
  theme(legend.position = "none", axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 5))
dev.off()
pdf("/Users//Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/nCount_Immune_per_patient.pdf", width = 14, height = 7)
give.nmedian <- function(x){
  return(c(y = median(x)*1.1, label = length(x))) 
}
VlnPlot(object = GC_all_immune, features = c("nCount_RNA"), group.by = "GC number", pt.size = 0) + 
  #geom_hline(yintercept=20, linetype='dashed') +
  stat_summary(fun.data = give.nmedian, geom = "text", fun = median, size = 3) +
  theme(legend.position = "none", axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 5))
dev.off()
####################################### Cluster Identity #########################
Idents(GC_all_immune) <- "seurat_clusters"
GC_all_immune_markers <- FindAllMarkers(GC_all_immune, only.pos = TRUE, min.pct = 0.4, logfc.threshold = 0.4)
write.table(GC_all_immune_markers, "/Users//Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/GC_all_immune_markers.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)
top50 <- GC_all_immune_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
pdf("/Users//Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Heatmap_markers_immune_top50.pdf", width = 15, height = 35)
DoHeatmap(GC_all_immune, features = top50$gene) + NoLegend()
dev.off()
pdf("/Users//Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/Immune_cell_types_IMMUNE.pdf", width = 7, height = 7)

###############################%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Jerby-Arnon (JA) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###############################
#######################################################################################################################################################
# Genes were taken from the summary tab - combined signature from Tirosh and Jerby-Arnon
############################### "Bcells_JA" #################################################
genes<- c("ADAM19", "ADAM28", "AFF3", "ATF7IP", "BACH2", "BANK1", "BCL11A", "BLK", "BLNK", "BTLA", "CCR6", "CD19", "CD1C", "CD22", "CD24", "CD37", "CD52", "CD79A", "CD79B", "CHMP7", "CIITA", "CLEC17A", "CNR2", "COL19A1", "COL4A3", "CR2", "CXCR5", "CYBASC3", "EEF1B2", "EEF1G", "EIF2S3", "ELK2AP", "FAIM3", "FAM129C", "FAM65B", "FCER2", "FCRL1", "FCRL2", "FCRL5", "FCRLA", "GGA2", "GNB2L1", "HLA-DOB", "HLA-DQA2", "HVCN1", "IGJ", "IGLL1", "IGLL3P", "IGLL5", "IRF8", "KBTBD8", "KIAA0125", "KIAA0226L", "LOC283663", "LRMP", "MS4A1", "NAPSB", "NCF1C", "NCOA3", "P2RX5", "PAX5", "PLEKHF2", "PNOC", "POLD4", "POU2AF1", "POU2F2", "PRKCB", "QRSL1", "RALGPS2", "RHOH", "SEL1L3", "SELL", "SMIM14", "SNX29", "SNX29P1", "SP110", "SP140", "SPIB", "ST6GAL1", "STAG3", "STAP1", "STRBP", "TCL1A", "TLR10", "TLR9", "TMEM154", "TNFRSF13B", "TP53INP1", "VPREB3", "WDFY4", "ZCCHC7", "CD2", "CD40", "CD5", "CD69", "CD70", "CD80", "CD86", "CD93", "PDCD1", "SDC1", "TNFRSF13C", "TNFRSF9", "TNFSF4")
geneSets <- GeneSet(genes, setName="Bcells_JA")
geneSets
cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Bcells_JA<-getAUC(cells_AUC)
Bcells_JA<-t(Bcells_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Bcells_JA)
FeaturePlot(GC_all_immune, features = "Bcells_JA", label = T)
VlnPlot(GC_all_immune, features = "Bcells_JA", pt.size = 0)

############################### "Basophils_JA" #################################################
genes<- c("ANPEP", "CCR3", "CD44", "CD63", "CD69", "ENPP3", "ICAM1", "IL3RA", "LAMP1", "TLR4")
geneSets <- GeneSet(genes, setName="Basophils_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Basophils_JA<-getAUC(cells_AUC)
Basophils_JA<-t(Basophils_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Basophils_JA)
FeaturePlot(GC_all_immune, features = "Basophils_JA", label = T)
VlnPlot(GC_all_immune, features = "Basophils_JA", pt.size = 0)
VlnPlot(GC_all_immune, features = "Basophils_JA", group.by = "orig.ident")

############################### "Eosinophils_JA" #################################################
genes<- c("C3AR1", "C5AR1", "CCR1", "CCR3", "CD244", "CD52", "CD53", "CXCR3", "FCER2", "FUT4", "IL9R", "ITGA4", "LAIR1", "PTGDR2", "S100A9", "SIGLEC10", "SIGLEC8")
geneSets <- GeneSet(genes, setName="Eosinophils_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Eosinophils_JA<-getAUC(cells_AUC)
Eosinophils_JA<-t(Eosinophils_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Eosinophils_JA)
FeaturePlot(GC_all_immune, features = "Eosinophils_JA", label = T)
VlnPlot(GC_all_immune, features = "Eosinophils_JA", pt.size = 0)
VlnPlot(GC_all_immune, features = "Eosinophils_JA", group.by = "orig.ident")


############################### "Mast_cells_JA" #################################################
genes<- c("ENPP3", "KIT")
geneSets <- GeneSet(genes, setName="Mast_cells_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Mast_cells_JA<-getAUC(cells_AUC)
Mast_cells_JA<-t(Mast_cells_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Mast_cells_JA)
FeaturePlot(GC_all_immune, features = "Mast_cells_JA", label = T)
VlnPlot(GC_all_immune, features = "Mast_cells_JA", pt.size = 0)

############################### "MDSC_JA" #################################################
genes<- c("CCR7", "CD1A", "CD1B", "CD1C", "CD207", "CD209", "CD4", "CD40", "CD80", "CD83", "CD86", "CMKLR1", "HLA-DOA", "HLA-DOB", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DRB6", "ITGA4", "ITGAM", "ITGAX", "LY75", "NRP1", "PDCD1LG2")
geneSets <- GeneSet(genes, setName="MDSC_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
MDSC_JA<-getAUC(cells_AUC)
MDSC_JA<-t(MDSC_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, MDSC_JA)
FeaturePlot(GC_all_immune, features = "MDSC_JA", label = T)
VlnPlot(GC_all_immune, features = "MDSC_JA", pt.size = 0)

############################### "Megakaryocyte_JA" #################################################
genes<- c("CD9", "GP1BA", "ITGA2B", "ITGAV", "ITGB3", "PECAM1", "SELP")
geneSets <- GeneSet(genes, setName="Megakaryocyte_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Megakaryocyte_JA<-getAUC(cells_AUC)
Megakaryocyte_JA<-t(Megakaryocyte_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Megakaryocyte_JA)
FeaturePlot(GC_all_immune, features = "Megakaryocyte_JA", label = T)
VlnPlot(GC_all_immune, features = "Megakaryocyte_JA", pt.size = 0)

############################### "Myeloid_DCs_JA" #################################################
genes<- c("CCR7", "CD1A", "CD1B", "CD1C", "CD207", "CD209", "CD4", "CD40", "CD80", "CD83", "CD86", "CMKLR1", "DCX", "ITGA4", "ITGAM", "ITGAX", "LY75", "NRP1", "PDCD1LG2")
geneSets <- GeneSet(genes, setName="Myeloid_DCs_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Myeloid_DCs_JA<-getAUC(cells_AUC)
Myeloid_DCs_JA<-t(Myeloid_DCs_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Myeloid_DCs_JA)
FeaturePlot(GC_all_immune, features = "Myeloid_DCs_JA", label = T)
VlnPlot(GC_all_immune, features = "Myeloid_DCs_JA", pt.size = 0)


############################### "Neutrophils_JA" #################################################
genes<- c("ANPEP", "C5AR1", "CD14", "CD33", "CEACAM8", "CSF3R", "CXCR1", "CXCR2", "FCGR1A", "FUT4", "ITGAM", "ITGAX", "MME", "PECAM1", "SELL", "TLR2")
geneSets <- GeneSet(genes, setName="Neutrophils_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Neutrophils_JA<-getAUC(cells_AUC)
Neutrophils_JA<-t(Neutrophils_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Neutrophils_JA)
FeaturePlot(GC_all_immune, features = "Neutrophils_JA", label = T)
VlnPlot(GC_all_immune, features = "Neutrophils_JA", pt.size = 0)


############################### "NKs_JA" #################################################
genes<- c("ALOX5AP", "APMAP", "CALM1", "CD160", "CD244", "CD247", "CLIC3", "CTSW", "FCRL6", "FGFBP2", "GNLY", "GZMB", "GZMM", "HOPX", "ID2", "IL18RAP", "IL2RB", "KIR2DL3", "KIR3DL2", "KLRB1", "KLRC1", "KLRD1", "KLRF1", "MATK", "MYBL1", "NCAM1", "NCR1", "NCR3", "NKG7", "NMUR1", "PRF1", "PTGDR", "PTPN4", "SAMD3", "SH2D1B", "TXK", "XCL1", "XCL2", "B3GAT1", "CD69", "ITGA2", "ITGAM", "ITGAX", "KLRA1", "KLRK1", "NKG2", "SIGLEC7", "SLAMF6", "SLAMF7")
geneSets <- GeneSet(genes, setName="NKs_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
NKs_JA<-getAUC(cells_AUC)
NKs_JA<-t(NKs_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, NKs_JA)
FeaturePlot(GC_all_immune, features = "NKs_JA", label = T)
VlnPlot(GC_all_immune, features = "NKs_JA", pt.size = 0)


############################### "pDCs_JA" #################################################
genes<- c("CCR7", "CD1A", "CD1B", "CD1C", "CD4", "CD40", "CD80", "CD83", "CD86", "CD8A", "CLEC4C", "CMKLR1", "IL3RA", "ITGA4", "ITGAM", "ITGAX", "NRP1", "PDCD1LG2", "TLR9")
geneSets <- GeneSet(genes, setName="pDCs_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
pDCs_JA<-getAUC(cells_AUC)
pDCs_JA<-t(pDCs_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, pDCs_JA)
FeaturePlot(GC_all_immune, features = "pDCs_JA", label = T)
VlnPlot(GC_all_immune, features = "pDCs_JA", pt.size = 0)

############################### "Macrophages_JA" #################################################
genes<- c("ABCA1", "ABI1", "ACAA1", "ACER3", "ACP2", "ACP5", "ACSL1", "ADAMDEC1", 
          "ADAP2", "ADORA3", "ADPGK", "AIF1", "AKR1A1", "ALDH2", "ALDH3B1", "AMICA1", "AMPD3", 
          "ANKRD22", "AP1B1", "APOC1", "AQP9", "ARAP1", "ARHGAP18", "ARHGAP27", "ARHGEF10L", "ARPC1B", 
          "ARRB2", "ASAH1", "ATF5", "ATG3", "ATG7", "ATP6AP1", "ATP6V0B", "ATP6V0D1", "ATP6V1B2", "ATP6V1F", 
          "BACH1", "BCKDHA", "BCL2A1", "BID", "BLOC1S1", "BLVRA", "BLVRB", "C10ORF54", "C11ORF75", "C15ORF48", "C19ORF38", 
          "C1ORF162", "C1QA", "C1QB", "C1QC", "C2", "C3AR1", "C5AR1", "C9ORF72", "CAPG", "CAPZA2", "CARD9", "CASP1", "CAT", 
          "CCDC88A", "CCR1", "CCR2", "CCRL2", "CD14", "CD163", "CD1D", "CD274", "CD300C", "CD300E", "CD300LB", "CD300LF", 
          "CD302", "CD33", "CD68", "CD80", "CD86", "CECR1", "CFD", "CFP", "CLEC10A", "CLEC12A", "CLEC4A", "CLEC4E", "CLEC5A", 
          "CLEC7A", "CMKLR1", "CMTM6", "CNDP2", "CNPY3", "CORO7", "CPVL", "CREG1", "CSF1R", "CSF2RA", "CSF3R", "CST3", "CSTA", "CTSA", 
          "CTSB", "CTSC", "CTSD", "CTSH", "CTSL1", "CTSS", "CXCL10", "CXCL16", "CXCL9", "CXCR2P1", "CYB5R4", "CYBA", "CYBB", "CYP2S1", 
          "DAPK1", "DBNL", "DENND1A", "DHRS9", "DMXL2", "DNAJC5B", "DOK1", "DOK3", "DPYD", "EBI3", "EMR2", "EPSTI1", "ETV6", "EVI2A", "F13A1", 
          "FAM105A", "FAM157B", "FAM26F", "FAM49A", "FAM96A", "FBP1", "FCER1G", "FCGR1A", "FCGR1B", "FCGR1C", "FCGR2A", "FCGR2C", "FCGR3B", "FCGRT", "FCN1", 
          "FERMT3", "FES", "FGL2", "FKBP15", "FLVCR2", "FOLR2", "FPR1", "FPR2", "FPR3", "FTH1", "FTL", "FUCA1", "FUOM", "GAA", "GABARAP", "GALC", "GATM", "GBP1", "GCA", "GGTA1P", "GK", "GLA", "GLB1", "GLRX", "GLUL", "GM2A", "GNA13", "GNA15", "GPBAR1", "GPR34", "GPR84", "GPX1", "GRN", "GSTO1", "H2AFY", "HCAR2", "HCAR3", "HCK", "HEIH", "HERPUD1", "HIST2H2BF", "HK2", "HK3", "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DPB2", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DRB6", "HMOX1", "HN1", "HPS1", "HSPA6", "HSPA7", "HSPBAP1", "IDH1", "IFI30", "IFI35", "IFIT2", "IFNGR1", "IFNGR2", "IGFLR1", "IGSF6", "IL10RB", "IL18", "IL1B", "IL1RN", "IL4I1", "IL8", "IRF5", "IRF7", "ITGAX", "JAK2", "KCNMA1", "KCNMB1", "KYNU", "LAIR1", "LAP3", "LGALS2", "LGALS9", "LGMN", "LILRA1", "LILRA2", "LILRA3", "LILRA4", "LILRA5", "LILRA6", "LILRB1", "LILRB2", "LILRB3", "LILRB4", "LILRB5", "LIPA", "LOC338758", "LOC729737", "LRRC25", "LST1", "LTA4H", "LYN", "LYZ", "M6PR", "MAFB", "MAN2B1", "MAPKAPK3", "MARCO", "MERTK", "MFSD1", "MGAT1", "MIF4GD", "MIIP", "MILR1", "MKNK1", "MNDA", "MOB1A", "MPEG1", "MPP1", "MRC1", "MS4A4A", "MS4A6A", "MS4A7", "MSR1", "MTHFD2", "MTMR14", "MX1", "MX2", "MXD1", "MYD88", "N4BP2L1", "NAAA", "NADK", "NAGA", "NAGK", "NAIP", "NCF2", "NCF4", "NCKAP1L", "NCOA4", "NFAM1", "NFKBID", "NINJ1", "NLRC4", "NLRP3", "NMI", "NOD2", "NPC2", "NPL", "NR1H3", "OAS1", "OAZ1", "OLR1", "OSCAR", "P2RX4", "P2RY12", "P2RY13", "P2RY14", "P2RY6", "PAK1", "PCK2", "PFKFB3", "PGD", "PILRA", "PLA2G15", "PLA2G7", "PLAUR", "PLBD1", "PLEK", "PLEKHO1", "PLEKHO2", "PLIN2", "PLXDC2", "PPM1M", "PPT1", "PRAM1", "PRKCD", "PSAP", "PSME2", "PTAFR", "PTPRE", "PYCARD", "RAB20", "RAB4B", "RAB8A", "RASGEF1B", "RASSF4", "RBM47", "RBPJ", "REEP4", "RELT", "RGS10", "RGS18", "RGS19", "RGS2", "RHBDF2", "RHOG", "RILPL2", "RIPK2", "RNASE6", "RNASEK", "RNASET2", "RNF13", "RNF130", "RNF144B", "RNF149", "RTN1", "S100A11", "S100A8", "S100A9", "SAMHD1", "SAT1", "SCAMP2", "SCIMP", "SCO2", "SCPEP1", "SDS", "SECTM1", "SEMA4A", "SERPINA1", "SERPINB1", "SFT2D1", "SGPL1", "SH3BGRL", "SHKBP1", "SIGLEC1", "SIGLEC14", "SIGLEC5", "SIGLEC7", "SIGLEC9", "SIRPA", "SIRPB1", "SIRPB2", "SKAP2", "SLAMF8", "SLC11A1", "SLC15A3", "SLC16A3", "SLC1A3", "SLC25A19", "SLC29A3", "SLC2A5", "SLC2A8", "SLC2A9", "SLC31A2", "SLC43A2", "SLC46A3", "SLC7A7", "SLC9A9", "SLCO2B1", "SMPDL3A", "SNX10", "SNX6", "SOD2", "SPI1", "SPINT2", "SQRDL", "SRC", "STX11", "STXBP2", "TALDO1", "TBXAS1", "TFRC", "TGFBI", "THEMIS2", "TIFAB", "TLR1", "TLR2", "TLR4", "TLR5", "TLR8", "TMEM106A", "TMEM144", "TMEM176A", "TMEM176B", "TMEM37", "TMEM51", "TMEM86A", "TNFAIP2", "TNFAIP8L2", "TNFSF13", "TNFSF13B", "TPP1", "TRAFD1", "TREM1", "TREM2", "TRPM2", "TTYH3", "TWF2", "TYMP", "TYROBP", "UBE2D1", "UBXN11", "UNC93B1", "VAMP8", "VMO1", "VSIG4", "WDFY2", "ZEB2", "ZNF267", "ZNF385A", "CCR5", "ENG", "FUT4", "ITGAL", "ITGAM", "LAMP2")
geneSets <- GeneSet(genes, setName="Macrophages_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Macrophages_JA<-getAUC(cells_AUC)
Macrophages_JA<-t(Macrophages_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Macrophages_JA)
FeaturePlot(GC_all_immune, features = "Macrophages_JA", label = T)
VlnPlot(GC_all_immune, features = "Macrophages_JA", pt.size = 0)

############################### "Cytotoxic_CD8_Tcells_JA" #################################################
genes<- c("CCL3", "CCL4", "CD2", "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CST7", "GZMA", "GZMB", "IFNG", "NKG7", "PRF1", "APOBEC3C", "B2M", "CCL5", "CSF1", "FASLG", "GZMH", "HLA-C")
geneSets <- GeneSet(genes, setName="Cytotoxic_CD8_Tcells_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Cytotoxic_CD8_Tcells_JA<-getAUC(cells_AUC)
Cytotoxic_CD8_Tcells_JA<-t(Cytotoxic_CD8_Tcells_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Cytotoxic_CD8_Tcells_JA)
FeaturePlot(GC_all_immune, features = "Cytotoxic_CD8_Tcells_JA", label = T)
VlnPlot(GC_all_immune, features = "Cytotoxic_CD8_Tcells_JA", pt.size = 0)

############################### "Tregs_JA" #################################################
genes<- c("CCR4", "CD4", "CNGB1", "CTLA4", "ENTPD1", "FOXP3", "IKZF2", "IL2RA", "ISG20", "ITGAE", "LAG3", "LRRC32", "NT5E", "SELL", "TNFRSF18", "TNFRSF4", "ARID3B", "CARD16", "CCR8", "CD177", "DDX60", "ERI1", "EZH2", "FAS", "GATA3", "GBP5", "GCHFR", "HNRNPC", "HPRT1", "IL10RA", "IL1R2", "IL2RG", "IL32", "LAIR2", "LAYN", "LOC541471", "LTB", "PARK7", "PIM2", "PMAIP1", "PPP2CA", "RASGRP1", "RORA", "RTKN2", "S100A4", "SAMD9", "STAM", "SUMO2", "USP15", "WDR1", "ZBTB32")
geneSets <- GeneSet(genes, setName="Tregs_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Tregs_JA<-getAUC(cells_AUC)
Tregs_JA<-t(Tregs_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Tregs_JA)
FeaturePlot(GC_all_immune, features = "Tregs_JA", label = T)
VlnPlot(GC_all_immune, features = "Tregs_JA", pt.size = 0)

############################### "CD4_Tcells_JA" #################################################
genes<- c("AIM1", "AQP3", "BCL11B", "CAMK4", "CCR4", "CCR8", "CD28", "CD4", "CD40LG", "CD5", "CD6", "DGKA", "DUSP16", "EML4", "F5", "FAAH2", "FAM102A", "FBLN7", "FLT3LG", "FOXP3", "FYB", "ICOS", "IL6R", "IL7R", "ITGB2-AS1", "ITK", "LAT", "LEPROTL1", "LOC100128420", "LOC285740", "MAF", "MAL", "PASK", "PBX4", "PBXIP1", "PIK3IP1", "PIM2", "SEPT6", "SLAMF1", "SPOCK2", "SUSD3", "TBC1D4", "TCF7", "TESPA1", "TIAM1", "TMEM66", "TNFAIP3", "TNFSF8", "TNIK", "TPT1")
geneSets <- GeneSet(genes, setName="CD4_Tcells_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
CD4_Tcells_JA<-getAUC(cells_AUC)
CD4_Tcells_JA<-t(CD4_Tcells_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, CD4_Tcells_JA)
FeaturePlot(GC_all_immune, features = "CD4_Tcells_JA", label = T)
VlnPlot(GC_all_immune, features = "CD4_Tcells_JA", pt.size = 0)

############################### "CD8_Tcells_JA" #################################################
genes<- c("ADORA2A", "ASB2", "ASXL2", "CBLB", "CCDC141", "CCDC64", "CD27", "CD3E", "CD84", "CD8A", "CD8B", "CLEC2D", "CRTAM", "CXCL13", "DNAJB1", "GPR171", "GZMK", "IFNG", "ITGA4", "ITGAE", "ITM2A", "JAKMIP1", "KLHL28", "LAG3", "LYST", "MAP4K1", "MCOLN2", "MIAT", "MIR155HG", "NELL2", "PAG1", "PCED1B", "PDCD1", "PRDM1", "PTPN7", "RAB27A", "RNF19A", "SIRPG", "SIT1", "SNAP47", "THEMIS", "TIGIT", "TIMD4", "TMEM155", "TNFRSF9", "TNIP3", "TOX", "TTC24", "TTN", "ZBED2")
geneSets <- GeneSet(genes, setName="CD8_Tcells_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
CD8_Tcells_JA<-getAUC(cells_AUC)
CD8_Tcells_JA<-t(CD8_Tcells_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, CD8_Tcells_JA)
FeaturePlot(GC_all_immune, features = "CD8_Tcells_JA", label = T)
VlnPlot(GC_all_immune, features = "CD8_Tcells_JA", pt.size = 0)

############################### "Tcells_JA" #################################################
genes<- c("ARHGEF1", "ASB2", "ATHL1", "BCL11B", "C16ORF54", "CASP8", "CCDC64", "CCND2", "CD2", "CD247", "CD27", "CD28", "CD3D", "CD3E", "CD3G", "CD5", "CD6", "CD7", "CD96", "CDC42SE2", "CELF2", "CNOT6L", "CORO1A", "CST7", "CTLA4", "CXCR3", "CXCR6", "CYTIP", "DEF6", "DENND2D", "EMB", "EVL", "FYB", "FYN", "GATA3", "GNG2", "GPR171", "GPR174", "GPRIN3", "GRAP2", "GZMA", "GZMM", "HNRNPA1P10", "ICOS", "IL12RB1", "IL21R", "IL2RB", "IL2RG", "IL32", "INPP4B", "IPCEF1", "ITGAL", "ITK", "JAK3", "KCNA3", "LAT", "LCK", "LIME1", "LOC100130231", "MBOAT1", "MIAT", "NLRC5", "PAG1", "PARP8", "PCED1B", "PCED1B-AS1", "PDCD1", "PIP4K2A", "PRDM1", "PRF1", "PRKCQ", "PTPN22", "PTPN7", "PTPRC", "PYHIN1", "RASAL3", "RASGRP1", "RGS1", "RHOF", "RNF213", "SCML4", "01-SEP", "SH2D1A", "SH2D2A", "SH3KBP1", "SIRPG", "SIT1", "SKAP1", "SLC9A3R1", "SPATA13", "SPN", "SPOCK2", "STAT4", "SYTL3", "TBC1D10C", "TC2N", "TESPA1", "THEMIS", "TIGIT", "TMEM66", "TNFAIP3", "TOX", "TRAF1", "TRAT1", "TTC39C", "TUBA4A", "UBASH3A", "WIPF1", "ZAP70", "ZC3HAV1")
geneSets <- GeneSet(genes, setName="Tcells_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Tcells_JA<-getAUC(cells_AUC)
Tcells_JA<-t(Tcells_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Tcells_JA)
FeaturePlot(GC_all_immune, features = "Tcells_JA", label = T)
VlnPlot(GC_all_immune, features = "Tcells_JA", pt.size = 0)

############################### "CD4_Tcells_Exhausted_JA" #################################################
genes<- c("AIM1", "AQP3", "BCL11B", "CAMK4", "CCR4", "CCR8", "CD28", "CD4", "CD40LG", "CD5", "CD6", "DGKA", "DUSP16", "EML4", "F5", "FAAH2", "FAM102A", "FBLN7", "FLT3LG", "FOXP3", "FYB", "ICOS", "IL6R", "IL7R", "ITGB2-AS1", "ITK", "LAT", "LEPROTL1", "LOC100128420", "LOC285740", "MAF", "MAL", "PASK", "PBX4", "PBXIP1", "PIK3IP1", "PIM2", "SEPT6", "SLAMF1", "SPOCK2", "SUSD3", "TBC1D4", "TCF7", "TESPA1", "TIAM1", "TMEM66", "TNFAIP3", "TNFSF8", "TNIK", "TPT1")
geneSets <- GeneSet(genes, setName="CD4_Tcells_Exhausted_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
CD4_Tcells_Exhausted_JA<-getAUC(cells_AUC)
CD4_Tcells_Exhausted_JA<-t(CD4_Tcells_Exhausted_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, CD4_Tcells_Exhausted_JA)
FeaturePlot(GC_all_immune, features = "CD4_Tcells_Exhausted_JA", label = T)
VlnPlot(GC_all_immune, features = "CD4_Tcells_Exhausted_JA", pt.size = 0)

############################### "CD4_Tcells_naive_JA" #################################################
genes<- c("ABLIM1", "ATM", "CAMK4", "CCR7", "EEF1A1", "EEF1B2", "EEF1G", "FAM65B", "FHIT", "GIMAP5", "IL7R", "LDHB", "LDLRAP1", "LEF1", "LOC100130231", "NAP1L1", "NOSIP", "PABPC1", "PIK3IP1", "SELL", "SERINC5", "SF1", "TCF7", "TMEM66", "TPT1", "TRABD2A", "TSC22D3", "TXNIP", "UBA52")
geneSets <- GeneSet(genes, setName="CD4_Tcells_naive_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
CD4_Tcells_naive_JA<-getAUC(cells_AUC)
CD4_Tcells_naive_JA<-t(CD4_Tcells_naive_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, CD4_Tcells_naive_JA)
FeaturePlot(GC_all_immune, features = "CD4_Tcells_naive_JA", label = T)
VlnPlot(GC_all_immune, features = "CD4_Tcells_naive_JA", pt.size = 0)

############################### "T_folicular_helper_cells_JA" #################################################
genes<- c("BCL6", "CD3D", "CD3E", "CD3G", "CD4", "CD40LG", "CD84", "CXCR5", "ICOS", "IL6R", "PDCD1", "SLAMF1", "STAT3", "TNFSF4")
geneSets <- GeneSet(genes, setName="T_folicular_helper_cells_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
T_folicular_helper_cells_JA<-getAUC(cells_AUC)
T_folicular_helper_cells_JA<-t(T_folicular_helper_cells_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, T_folicular_helper_cells_JA)
FeaturePlot(GC_all_immune, features = "T_folicular_helper_cells_JA", label = T)
VlnPlot(GC_all_immune, features = "T_folicular_helper_cells_JA", pt.size = 0)

############################### "Th1_cells_JA" #################################################
genes<- c("CCR1", "CCR5", "CD4", "CSF2", "CXCR3", "DPP4", "HAVCR2", "IFNA1", "IFNGR1", "IL2", "KLRD1", "TNF", "TNFSF11")
geneSets <- GeneSet(genes, setName="Th1_cells_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Th1_cells_JA<-getAUC(cells_AUC)
Th1_cells_JA<-t(Th1_cells_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Th1_cells_JA)
FeaturePlot(GC_all_immune, features = "Th1_cells_JA", label = T)
VlnPlot(GC_all_immune, features = "Th1_cells_JA", pt.size = 0)

############################### "Th2_cells_JA" #################################################
genes<- c("CCR3", "CCR4", "CCR7", "CCR8", "CD4", "CSF2", "CXCR4", "GATA3", "HAVCR1", "ICOS", "IL10", "IL13", "IL1R1", "IL4", "IL5", "IL6", "PTGDR2")
geneSets <- GeneSet(genes, setName="Th2_cells_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Th2_cells_JA<-getAUC(cells_AUC)
Th2_cells_JA<-t(Th2_cells_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Th2_cells_JA)
FeaturePlot(GC_all_immune, features = "Th2_cells_JA", label = T)
VlnPlot(GC_all_immune, features = "Th2_cells_JA", pt.size = 0)

############################### "Th9_cells_JA" #################################################
genes<- c("CD3D", "CD3E", "CD3G", "CD4", "GATA3", "IRF4", "STAT6")
geneSets <- GeneSet(genes, setName="Th9_cells_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Th9_cells_JA<-getAUC(cells_AUC)
Th9_cells_JA<-t(Th9_cells_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Th9_cells_JA)
FeaturePlot(GC_all_immune, features = "Th9_cells_JA", label = T)
VlnPlot(GC_all_immune, features = "Th9_cells_JA", pt.size = 0)

############################### "Th17_cells_JA" #################################################
genes<- c("CCR4", "CCR6", "CD38", "CD3D", "CD3E", "CD3G", "CD4", "IL17A", "IL17F", "IL1R1", "IL21", "IL22", "KLRB1", "LINC-ROR", "STAT3")
geneSets <- GeneSet(genes, setName="Th17_cells_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Th17_cells_JA<-getAUC(cells_AUC)
Th17_cells_JA<-t(Th17_cells_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Th17_cells_JA)
FeaturePlot(GC_all_immune, features = "Th17_cells_JA", label = T)
VlnPlot(GC_all_immune, features = "Th17_cells_JA", pt.size = 0)

############################### "Th22_cells_JA" #################################################
genes<- c("AHR", "CCR10", "CCR4", "CCR6", "CD3D", "CD3E", "CD3G", "CD4")
geneSets <- GeneSet(genes, setName="Th22_cells_JA")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Th22_cells_JA<-getAUC(cells_AUC)
Th22_cells_JA<-t(Th22_cells_JA)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Th22_cells_JA)
FeaturePlot(GC_all_immune, features = "Th22_cells_JA", label = T)
VlnPlot(GC_all_immune, features = "Th22_cells_JA", pt.size = 0)

###############################%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sade-Feldman (SF) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#############################
#######################################################################################################################################################
#Only top 100 genes (based on p-value) taken for signature
############################### "Bcells_SF" #################################################
genes<- c("IGHD", "PAX5", "FCRL1", "CR2", "VPREB3", "FCER2", "CD19", "EBF1", "CD22", "BANK1", "CLEC17A", "FCRLA", "FCRL2", "MS4A1", "BLK", "RALGPS2", "TCL1A", "TLR10", "FAM129C", "CNR2", "ARHGAP24", "HLA-DOB", "IGHM", "KIAA0226L", "CD79A", "BCL11A", "STAP1", "AFF3", "SPIB", "SWAP70", "ADAM28", "BLNK", "IGKC", "IGHG3", "WDFY4", "MEF2C", "IGHG2", "SMIM14", "POU2AF1", "CD79B", "CD40", "HVCN1", "TCF4", "CXCR5", "IGHG1", "NCF1", "FGD2", "IRF8", "CD83", "CIITA", "LY9", "HLA-DQB2", "HLA-DQA2", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "PKIG", "HLA-DOA", "HLA-DMB", "FCRL5", "SELL", "KIAA0125", "IGLC3", "CYBB", "RABEP2", "CCR6", "COBLL1", "BTK", "ADAM19", "PLCG2", "HLA-DPA1", "ALOX5", "AKAP2", "CD72", "HLA-DRB5", "SYK", "MICAL3", "SNX2", "CD55", "HLA-DRB1", "RASGRP3", "SIGLEC14", "HLA-DMA", "HLA-DPB1", "RAB30", "LAT2", "IGLL5", "PLAC8", "FAIM3", "STX7", "TPD52", "PHACTR1", "PDLIM1", "ST6GAL1", "SNX29", "POU2F2", "IGLC2", "PIK3C2B", "MGAT5", "LY86")
geneSets <- GeneSet(genes, setName="Bcells_SF")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Bcells_SF<-getAUC(cells_AUC)
Bcells_SF<-t(Bcells_SF)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Bcells_SF)
FeaturePlot(GC_all_immune, features = "Bcells_SF", label = T)
VlnPlot(GC_all_immune, features = "Bcells_SF", pt.size = 0)

############################### "Plasma_cells_SF" #################################################
genes<- c("SDC1", "DERL3", "IGLC1", "MZB1", "IGLL5", "TNFRSF17", "POU2AF1", "IGHA1", "IGHG4", "CD79A", "PRDX4", "IGLV3-1", "FKBP11", "IGHG2", "IGHG3", "HID1", "CPNE5", "XBP1", "ITM2C", "IGHG1", "FKBP2", "SPAG4", "SLC17A9", "TRAM2", "SDF2L1", "PDK1", "IGKC", "VIMP", "IGKV3OR2-268", "SEC11C", "IGHV3OR16-9", "RAB30", "P2RX1", "PPAPDC1B", "CRELD2", "SYVN1", "IGHV3-48", "ANKRD28", "COBLL1", "IGHA2", "TPD52", "IGLV6-57", "BLNK", "GAB1", "IGHV3-7", "RRBP1", "SELM", "PYCR1", "MEI1", "EAF2", "SLAMF7", "IGHV3-21", "SEL1L", "SPATS2", "C19ORF10", "IGLC7", "CHPF", "SEL1L3", "LMAN1", "CLIC4", "ST6GALNAC4", "UBE2J1", "TXNDC11", "FNDC3B", "MANEA", "IGHV3-11", "IGHV3-23", "HM13", "SIL1", "ERLEC1", "IGKV4-1", "KDELR1", "IGLC2", "PNOC", "DNAJB9", "KDELR2", "PLCG2", "PDIA4", "HDLBP", "CREB3L2", "IGKV3D-15", "ZBP1", "MANF", "CD38", "TMEM258", "ST6GAL1", "SSR3", "TXNDC15", "CASP10", "ARSA", "CLPTM1L", "SSR4", "SRPRB", "SPCS2", "IGLC3", "IGHV1-69", "FNDC3A", "C11orf24", "SPCS3")
geneSets <- GeneSet(genes, setName="Plasma_cells_SF")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Plasma_cells_SF<-getAUC(cells_AUC)
Plasma_cells_SF<-t(Plasma_cells_SF)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Plasma_cells_SF)
FeaturePlot(GC_all_immune, features = "Plasma_cells_SF", label = T)
VlnPlot(GC_all_immune, features = "Plasma_cells_SF", pt.size = 0)

############################### "Monocytes_Macrophages_SF" #################################################
genes<- c("MARCO", "FPR3", "CLEC5A", "AQP9", "CLEC10A", "HNMT", "OLR1", "CD300E", "ARHGEF10L", "HK3", "ANPEP", "LGALS2", "FCN1", "GPR84", "CCL2", "CLEC4E", "CSTA", "TREM1", "FCGR1B", "TMEM176B", "C19ORF59", "OSCAR", "CYP2S1", "STAB1", "VCAN", "FPR1", "SERPINA1", "TLR2", "VSIG4", "TMEM176A", "VMO1", "MAFB", "FCGR1A", "SLC37A2", "CD14", "TLR4", "CD300LF", "CXCL3", "CD33", "GPNMB", "FOLR2", "LILRA3", "LPCAT2", "TREM2", "PLXDC2", "CD163", "RIN2", "NFAM1", "CD300C", "IL1B", "MSR1", "LILRA6", "ALDH3B1", "CXCL2", "C5AR1", "ADAP2", "IL8", "C15ORF48", "MRC1L1", "BST1", "PLAU", "SIGLEC9", "MS4A4A", "SIRPA", "RNASE1", "CPVL", "PLBD1", "LILRA2", "SLC11A1", "PYGL", "HCK", "LRP1", "C1QC", "LYZ", "LRRC25", "ADM", "C1QA", "APOBEC3A", "ZNF385A", "TGM2", "TNFAIP2", "CSF3R", "SLCO2B1", "LILRB2", "SIGLEC7", "TNS1", "S100A9", "RAB20", "EPB41L3", "C1QB", "PTAFR", "CD93", "PRAM1", "ST3GAL6", "S100A8", "IL1RN", "SLC8A1", "FBP1", "LILRA1")
geneSets <- GeneSet(genes, setName="Monocytes_Macrophages_SF")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Monocytes_Macrophages_SF<-getAUC(cells_AUC)
Monocytes_Macrophages_SF<-t(Monocytes_Macrophages_SF)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Monocytes_Macrophages_SF)
FeaturePlot(GC_all_immune, features = "Monocytes_Macrophages_SF", label = T)
VlnPlot(GC_all_immune, features = "Monocytes_Macrophages_SF", pt.size = 0)

############################### "DCs_SF" #################################################
genes<- c("CLEC4C", "LILRA4", "PTPRS", "IL3RA", "SERPINF1", "PLD4", "C1ORF186", "DNASE1L3", "TSPAN13", "SPIB", "SMPD3", "EPHB1", "EGLN3", "MAP1A", "TNFRSF21", "PTCRA", "LILRB4", "BCL11A", "TCF4", "SCAMP5", "MPEG1", "TLR9", "ZFAT", "PFKFB2", "TGFBI", "DERL3", "CST3", "PACSIN1", "CSF2RB", "SLC15A4", "SLC12A3", "TPM2", "RNASE6", "APP", "IGKC", "CLIC3", "CD68", "TYROBP", "NRP1", "P2RY14", "FAM129C", "KRT5", "CCDC50", "PTGDS", "FCER1G", "SOX4", "WDFY4", "CYB561A3", "P2RY6", "GAS6", "PLAC8", "ITM2C", "TRAF4", "IRF8", "MYBL2", "SMIM5", "MZB1", "CHAF1A", "SCN9A", "LGMN", "GRN", "TLR7", "PLVAP", "CYBB", "UGCG", "CORO1C", "NPC2", "BLNK", "NOTCH4", "CBFA2T3", "SUSD1", "PPP1R14B", "THEMIS2", "RNF130", "POLB", "CIITA", "GAPT", "IRF7", "FGD2", "IRF4", "PLEKHD1", "SULF2", "CD36", "CD4", "FCHSD2", "IGHM", "CSF2RA", "GAB1", "GPR183", "SH2B3", "MS4A6A", "CCDC88A", "PTPRE", "RRBP1", "NEK8", "PLEK", "GNA15", "TCL1A", "AC023590.1")
geneSets <- GeneSet(genes, setName="DCs_SF")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
DCs_SF<-getAUC(cells_AUC)
DCs_SF<-t(DCs_SF)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, DCs_SF)
FeaturePlot(GC_all_immune, features = "DCs_SF", label = T)
VlnPlot(GC_all_immune, features = "DCs_SF", pt.size = 0)

############################### "CD8_Tcells_exhausted_SF" #################################################
genes<- c("PDCD1", "PRF1", "CD8A", "NKG7", "GZMA", "CCL4", "KLRK1", "CCL5", "CST7", "SIRPG", "GZMK", "CD8B", "HAVCR2", "CD3G", "CTSW", "CCL3", "IFNG", "CD3D", "CCL4L2", "TIGIT", "CD27", "CXCR6", "GZMH", "LYST", "LAG3", "CD2", "CD3E", "APOBEC3G", "TRAC", "IL32", "CD38", "KLRC4-KLRK1", "IL2RB", "IKZF3", "KLRD1", "GBP5", "GZMB", "FASLG", "RARRES3", "TOX", "SLAMF7", "PTPRCAP", "GIMAP4", "KLRC4", "CBLB", "PRKCH", "LCK", "VCAM1", "SNAP47", "GIMAP7", "ITM2A", "PYHIN1", "INPP4B", "CCR5", "IFITM1", "F2R", "GIMAP6", "DTHD1", "LCP2", "CHST12", "APOBEC3C", "CXCL13", "PVRIG", "OASL", "SIT1", "DENND2D", "SH2D1A", "TRGC2", "SLFN5", "RAB27A", "IFI6", "APOBEC3D", "GIMAP5", "LY6E", "HCST", "CD96", "CD84", "GPR174", "SLFN12L", "GBP2", "PSTPIP1", "PRDM1", "CD7", "FYN", "GPR171", "ZAP70", "TBC1D10C", "CXCR3", "EVL", "ID2", "GPR56", "PSMB9", "ITGAE", "S100PBP", "JAKMIP1", "THEMIS", "MAP4K1", "SH2D2A", "FCRL3")
geneSets <- GeneSet(genes, setName="CD8_Tcells_exhausted_SF")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
CD8_Tcells_exhausted_SF<-getAUC(cells_AUC)
CD8_Tcells_exhausted_SF<-t(CD8_Tcells_exhausted_SF)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, CD8_Tcells_exhausted_SF)
FeaturePlot(GC_all_immune, features = "CD8_Tcells_exhausted_SF", label = T)
VlnPlot(GC_all_immune, features = "CD8_Tcells_exhausted_SF", pt.size = 0)

############################### "Tregs_SF" #################################################
genes<- c("FOXP3", "TNFRSF4", "TNFRSF18", "MAF", "CD4", "SPOCK2", "CTLA4", "ICOS", "ICA1", "RTKN2", "BATF", "IL2RA", "TBC1D4", "TRAC", "ARID5B", "TIGIT", "FBLN7", "TMEM173", "CD28", "ZC3H12D", "PBXIP1", "CD2", "DUSP4", "MAGEH1", "TNFRSF25", "IL32", "LTB", "PHTF2", "GK", "IKZF2", "STAM", "CORO1B", "PIM2", "HS3ST3B1", "TIAM1", "CD5", "ETV7", "RORA", "PHACTR2", "MBOAT1", "DUSP16", "HNRNPLL", "CTSB", "PELI1", "TLK1", "CD3D", "THADA", "SIRPG", "FAS", "DNPH1", "RHBDD2", "SLAMF1", "KLRB1", "IL6ST", "LCK", "SDC4", "PIK3IP1", "ITM2A", "BTLA", "CD247", "ITK", "PHLDA1", "LAT", "CD82", "SOD1", "RAB11FIP1", "TNIK", "SKAP1", "ZC3H7A", "IL2RB", "LBH", "LIMS1", "MICAL2", "TP53INP1", "PBX4", "EPSTI1", "BIRC3", "CD7", "GOLGA8B", "CALM3", "CD27", "TFRC", "MSI2", "NR3C1", "PRDM1", "IL7R", "P2RY10", "LY75", "GBP2", "DGKA", "NABP1", "SH2D2A", "CNST", "FAIM3", "NDFIP1", "ADAM19", "TRAF3", "GPRIN3", "NCF4")
geneSets <- GeneSet(genes, setName="Tregs_SF")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Tregs_SF<-getAUC(cells_AUC)
Tregs_SF<-t(Tregs_SF)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Tregs_SF)
FeaturePlot(GC_all_immune, features = "Tregs_SF", label = T)
VlnPlot(GC_all_immune, features = "Tregs_SF", pt.size = 0)

############################### "Cytotoxic_lymphocytes_SF" #################################################
genes<- c("NKG7", "CCL5", "GZMA", "FGFBP2", "SAMD3", "GNLY", "GZMH", "TGFBR3", "KLRG1", "AOAH", "PRF1", "CST7", "KLRK1", "CTSW", "FCRL6", "A2M", "CCL4", "TRGC1", "ZAP70", "GZMM", "KLRD1", "TRGC2", "CCL4L2", "CCL4L1", "GZMB", "CD8A", "SPON2", "TRDC", "ANXA1", "TC2N", "SYNE1", "PXN", "GRAP2", "HOPX", "GIMAP7", "FCGR3A", "PYHIN1", "PTPN4", "SCML4", "MYO1F", "S1PR1", "SLAMF7", "PLAC8", "STAT4", "STOM", "KLRB1", "KLRC4-KLRK1", "GZMK", "SLC9A3R1", "C20ORF112", "SLFN12L", "TSEN54", "RORA", "BIN2", "GIMAP5", "ITM2C", "C12ORF75", "SLFN5", "PIM1", "SORL1", "NLRC3", "C5ORF56", "RASAL3", "GLIPR2", "MGAT4A", "SYTL1", "TBCD", "GPR56", "THEMIS", "TPST2", "TNF", "PATL2")
geneSets <- GeneSet(genes, setName="Cytotoxic_lymphocytes_SF")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Cytotoxic_lymphocytes_SF<-getAUC(cells_AUC)
Cytotoxic_lymphocytes_SF<-t(Cytotoxic_lymphocytes_SF)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Cytotoxic_lymphocytes_SF)
FeaturePlot(GC_all_immune, features = "Cytotoxic_lymphocytes_SF", label = T)
VlnPlot(GC_all_immune, features = "Cytotoxic_lymphocytes_SF", pt.size = 0)


############################### "CD8_Tcells_exhausted_HS_SF" #################################################
genes<- c("VCAM1", "DUSP4", "CCL5", "NKG7", "TNFRSF9", "PRF1", "CD8A", "NAB1", "CTLA4", "LYST", "HAVCR2", "CXCL13", "CST7", "GZMB", "KLRK1", "PRDM1", "CRTAM", "CBLB", "GPR56", "PDCD1", "PHLDA1", "DFNB31", "CREM", "CD8B", "CCL4", "TNFAIP3", "RGS1", "CCL4L1", "SYTL3", "GZMH", "RNF19A", "CXCR6", "HNRNPLL", "KIR2DL4", "TMEM2", "ATXN1", "DTHD1", "GZMA", "LAG3", "GEM", "ZNF331", "CTSW", "SAMSN1", "SH2D2A", "TOX", "TTN", "TIGIT", "CD7", "CD27", "FYN", "CCDC64", "PDE4D", "TNFSF9", "RAB27A", "SLC7A5", "PDE3B", "MYO7A", "RGS2", "NELL2", "CNOT6L", "CCL4L2", "SLA2", "ENTPD1", "CADM1", "KLRC4-KLRK1", "DUSP2", "VPS37B", "KLRC4", "PRKCH", "SNAP47", "PTPN22", "ETS1", "PAM", "AHI1", "METRNL", "GZMK", "SLA", "FAM46C", "ITM2A", "MCTP2", "GABARAPL1", "RUNX3", "GOLIM4", "PFKFB3", "ASXL2", "TSPYL2", "KLRD1", "CLEC2D", "TUBA4A", "TGIF1", "ITGA4", "CCND2", "DNAJA1", "APOBEC3G", "TNFRSF1B", "BTG3", "HSPH1", "BHLHE40", "SIRPG")
geneSets <- GeneSet(genes, setName="CD8_Tcells_exhausted_HS_SF")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
CD8_Tcells_exhausted_HS_SF<-getAUC(cells_AUC)
CD8_Tcells_exhausted_HS_SF<-t(CD8_Tcells_exhausted_HS_SF)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, CD8_Tcells_exhausted_HS_SF)
FeaturePlot(GC_all_immune, features = "CD8_Tcells_exhausted_HS_SF", label = T)
VlnPlot(GC_all_immune, features = "CD8_Tcells_exhausted_HS_SF", pt.size = 0)

############################### "Memory_Tcells_SF" #################################################
genes<- c("TCF7", "IL7R", "LEF1", "CCR7", "TNFRSF25", "DGKA", "SORL1", "CAMK4", "S1PR1", "LTB", "SERINC5", "GIMAP5", "FLT3LG", "PLAC8", "TC2N", "SELL", "FAM65B", "MGAT4A", "PASK", "ABLIM1", "C20ORF112", "PIK3IP1", "OXNAD1", "TESPA1", "FAM102A", "DENND2D", "FOXP1", "DHRS3", "GIMAP7", "CD5", "RCAN3", "TMEM63A", "SCML4", "ATM", "SATB1", "CD28", "CCDC109B", "AAK1", "ICAM2", "KLRB1", "EPB41", "GPR183", "NOSIP", "GOLGA8A", "TTC39C", "GOLGA8B", "TMEM123", "CHMP7", "NELL2", "RAPGEF6", "KIAA0922", "GIMAP2")
geneSets <- GeneSet(genes, setName="Memory_Tcells_SF")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Memory_Tcells_SF<-getAUC(cells_AUC)
Memory_Tcells_SF<-t(Memory_Tcells_SF)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Memory_Tcells_SF)
FeaturePlot(GC_all_immune, features = "Memory_Tcells_SF", label = T)
VlnPlot(GC_all_immune, features = "Memory_Tcells_SF", pt.size = 0)

############################### "Lymphocytes_exhausted_SF" #################################################
genes<- c("SPC25", "CDCA5", "KIF15", "CDC45", "DLGAP5", "HIST1H3G", "KIF18B", "RRM2", "UBE2C", "HJURP", "ESCO2", "SPC24", "BIRC5", "CDC6", "CDCA8", "AURKB", "ZWINT", "GTSE1", "DTL", "RAD51", "CDCA3", "MELK", "CKAP2L", "ANLN", "ASF1B", "TYMS", "NCAPG", "TK1", "PKMYT1", "KIFC1", "KIAA0101", "CCNB2", "CDC20", "TROAP", "CLSPN", "ASPM", "GINS2", "KIF23", "KIF2C", "RAD51AP1", "NUF2", "TOP2A", "CDK1", "MKI67", "MLF1IP", "DHFR", "KIF11", "CENPW", "TPX2", "CASC5", "CDKN3", "CCNA2", "BUB1B", "MCM2", "UBE2T", "MCM4", "TCF19", "BUB1", "FEN1", "WDR34", "NCAPG2", "CENPF", "MAD2L1", "NCAPH", "FANCI", "CENPM", "CDCA7", "RNASEH2A", "SMC2", "STMN1", "CKS1B", "WDR76", "NUSAP1", "CENPK", "TMEM106C", "NCAPD2", "PCNA", "LIG1", "MCM3", "MCM7", "RRM1", "DUT", "DNAJC9", "MCM6", "HELLS", "FANCD2", "KIF22", "MCM5", "SAE1", "RANBP1", "FABP5", "TMPO", "H2AFV", "SNRPD1", "CBX5", "HIRIP3", "NUDT1", "ANP32B", "TIMELESS")
geneSets <- GeneSet(genes, setName="Lymphocytes_exhausted_SF")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_all_immune@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Lymphocytes_exhausted_SF<-getAUC(cells_AUC)
Lymphocytes_exhausted_SF<-t(Lymphocytes_exhausted_SF)
GC_all_immune@meta.data<-cbind(GC_all_immune@meta.data, Lymphocytes_exhausted_SF)
FeaturePlot(GC_all_immune, features = "Lymphocytes_exhausted_SF", label = T)
VlnPlot(GC_all_immune, features = "Lymphocytes_exhausted_SF", pt.size = 0)
dev.off()

