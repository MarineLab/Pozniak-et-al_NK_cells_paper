#R4.4
library(Seurat)
library(ggplot2)
library(dplyr)
library(nichenetr)
library(colorRamp2)
library(AUCell)
library(tidyverse)
library(GSEABase)
library(data.table)
library(SeuratWrappers)
library(ggpubr)
library(RColorBrewer)
setwd("/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/YUMM_NRAS")
dir <- "/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/YUMM_NRAS"
list.files(dir) -> samples
seurat_list <- list()
#samp <- samples[1]

for(samp in samples){
  data_10x <- Read10X_h5(filename = paste(dir,samp,"/filtered_feature_bc_matrix.h5",sep = "/"))
  
  seu <- CreateSeuratObject(data_10x,min.cells = 10, min.features = 200,project = samp)
  rm(data_10x)
  
  seu[["percent.mt"]]  <- PercentageFeatureSet(seu, pattern = "mt-")
  seu <- SCTransform(seu, verbose = TRUE, vars.to.regress = c("percent.mt"))
  
  #Generate sce object
  as.SingleCellExperiment(seu,assay = "RNA") -> sce
  #Run doublet finder witn default parameters
  scDblFinder(sce) -> sce
  #Save classification in seurat object
  colData(sce) %>% data.frame() %>% .$scDblFinder.class -> seu@meta.data$Doublets
  seurat_list[[samp]] <- seu
  rm(seu,sce)
}

#Merge and add metadata
merged_seu <- merge(seurat_list$NRO001,y = c(seurat_list$NRO002,seurat_list$NRO003,seurat_list$NRO004,seurat_list$NRO005,seurat_list$NRO006,seurat_list$NRO007,seurat_list$NRO008),project="NK_mouse")
merged_seu <- subset(merged_seu, subset = Doublets == "singlet")

#Rename seurat clusters based on manual curation of clusters
merged_seu@meta.data$"Sample_type" <- plyr::revalue(as.character(merged_seu$orig.ident),
                                               c("NRO001" = "YUMM52aPD1",
                                                 "NRO002" = "YUMM52aNK11",
                                                 "NRO003" = "YUMM52control",
                                                 "NRO004" = "YUMM52aPD1aNK11",
                                                 "NRO005" = "NRAScontrol",
                                                 "NRO006" = "NRASaPD1",
                                                 "NRO007" = "NRASaNK11",
                                                 "NRO008" = "NRASaPD1aNK11")) 
Idents(merged_seu) <- merged_seu@meta.data$Sample_type

pdf("plots/QC.pdf", width = 12, height = 7)
VlnPlot(merged_seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) 
dev.off()

#Subset by quality
seu_sub <- subset(merged_seu, subset = nFeature_RNA > 500 &  nFeature_RNA < 6000 &  percent.mt < 10)
pdf("plots/QC_subsets.pdf", width = 12, height = 7)
VlnPlot(seu_sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) 
dev.off()

#saveRDS(merged_seu,"merged.rds")
#saveRDS(seu_sub,"subset.rds")

seu <- readRDS("/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/YUMM_NRAS/subset.rds")
mouse_cell_cycle_genes <- readRDS("/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/mouse_cell_cycle_genes.rds")

s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seu <- SCTransform(seu, verbose = TRUE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'),vst.flavor = "v2")
seu <- RunPCA(seu , verbose = TRUE,npcs = 100)
ElbowPlot(seu,ndims = 50)
seu  <- FindNeighbors(seu , dims = 1:30,reduction = "pca")
seu  <- RunUMAP(seu , dims=1:30,reduction.name = "umap.unintegrated")
DimPlot(seu,group.by = "Sample_type",reduction = "umap.unintegrated")

##################### SCVI
seu <- IntegrateLayers(object = seu, method = scVIIntegration, new.reduction = "integrated.scvi",normalization.method = "SCT",
                          conda_env = "/Users/micromamba/envs/scvi-env", verbose = FALSE,group.by = "orig.ident")
ElbowPlot(seu)
seu <- FindNeighbors(seu, reduction = "integrated.scvi", dims = 1:30)#25
seu <- RunUMAP(seu, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")
seu <- FindClusters(seu, resolution = 0.3, cluster.name = "seurat_clusters")

DimPlot(seu,reduction = "umap.scvi", label = TRUE)
seu <- PrepSCTFindMarkers(seu)
All_Markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
write.table(All_Markers, "/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/YUMM_NRAS/All_Markers.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)
pdf("/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/YUMM_NRAS/UMAP_TME.pdf", width = 8, height = 6)
DimPlot(seu,group.by = "Sample_type",reduction = "umap.scvi")+NoAxes()
dev.off()


##########################          Immune          ##########################
genes4 <- c("Acap1",  "Akna",  "Alox5ap",  "Ankrd44",  "Apobec3g",  "Arhgap15",  "Arhgap25",  "Arhgap30",  "Arhgap4",  "Arhgap9",  "Arhgdib",  "Atp2a3",  "Bin2",  "C16orf54",  "Ccdc88b",  "Cd37",  "Cd48",  "Cd52",  "Cd53",  "Cd69",  "Cd84",  "Cdc42se2",  "Celf2",  "Cntrl",  "Coro1a",  "Csk",  "Cxcr4",  "Cyth4",  "Cytip",  "Def6",  "Dennd1c",  "Dock2",  "Dock8",  "Dusp2",  "Evi2b",  "Fermt3",  "Fgd3",  "Fnbp1",  "Gbp5",  "Gpr65",  "Gpsm3",  "Hcls1",  "Hmha1",  "Ikzf1",  "Il10ra",  "Il16",  "Il2rg",  "Inpp5d",  "Itga4",  "Itgal",  "Itgb2",  "Lair1",  "Laptm5",  "Lcp1",  "Lilrb3",  "Limd2",  "Lpxn",  "Lsp1",  "Ly9",  "Map4k1",  "Myo1g",  "Nckap1l",  "Nr4a2",  "Parp8",  "Parvg",  "Pik3cd",  "Pim2",  "Plcb2",  "Plekha2",  "Prkcb",  "Psd4",  "Pstpip2",  "Ptk2b",  "Ptpn22",  "Ptpn6",  "Ptpn7",  "Ptprc",  "Rac2",  "Rassf5",  "Rcsd1",  "Rgs1",  "Rhoh",  "Rps6ka1",  "Samsn1",  "Sash3",  "Sla",  "Snx20",  "Sp140",  "Stk17b",  "Tagap",  "Tbc1d10c",  "Tmc6",  "Tmc8",  "Tmsb4x",  "Traf3ip3",  "Tsc22d3",  "Tstd1",  "Ucp2",  "Vav1",  "Wipf1")
geneSets <- GeneSet(genes4, setName="Immune")
geneSets
cells_rankings <- AUCell_buildRankings(seu@assays[["SCT"]]@counts,splitByBlocks=TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Immune<-getAUC(cells_AUC)
Immune<-t(Immune)
seu@meta.data<-cbind(seu@meta.data, Immune)
pdf("/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/YUMM_NRAS/UMAP_Immune_signature.pdf", width = 7, height = 5)
FeaturePlot(seu, features = c ("Immune"), reduction = "umap.scvi")+NoAxes()
dev.off()

##########################          Malignant          ##########################
genes4 <- c("Aasdhppt", "Aatf", "Acn9", "Acot7", "Acsl3", "Adipor1", "Adsl", "Ahcy", "Aif1l", "Ak2", "Aldoa", "Alx1", "Amz2", "Anapc11", "Ankrd54", "Anp32a", "Ap1s2", "Apeh", "Apoa1bp", "Apoc2", "Apod", "Apoo", "Arpc1a", "Atic", "Atp1a1", "Atp5c1", "Atp5g1", "Atp5g2", "Atp5g3", "Atp6v1c1", "Atp6v1e1", "Atp6v1g1", "Baiap2", "Bancr", "Bcan", "Bcas3", "Bcl2l13", "Birc7", "Bzw2", "C10orf90", "C11orf31", "C17orf89", "C1orf43", "C1orf85", "C1qbp", "C22orf32", "C4orf48", "Ca14", "Ca8", "Cacybp", "Capn3", "Cbx3", "Ccnd1", "Cct2", "Cct3", "Cct4", "Cct6a", "Cct8", "Cdh19", "Cdh3", "Cdk2", "Cdk4", "Cep170", "Chchd6", "Chd1l", "Clcn7", "Clns1a", "Cmc2", "Coa3", "Coa4", "Coa6", "Cox5b", "Cox6a1", "Cox7a2", "Cox7a2l", "Cox7c", "Cox8a", "Csag1", "Csag2", "Csag3", "Cspg4", "Cyc1", "Cyp27a1", "Daam2", "Dancr", "Dap3", "Dct", "Dcxr", "Ddit3", "Ddt", "Dll3", "Dnah14", "Dnaja4", "Drg1", "Ednrb", "Eif3c", "Eif3d", "Eif3e", "Eif3h", "Eif3l", "Eif3m", "Eno1", "Eno2", "Entpd6", "Epb41l4a-as1", "Erbb3", "Esrp1", "Etv4", "Etv5", "Exosc4", "Fah", "Fahd2b", "Fam103a1", "Fam162a", "Fam178b", "Farp2", "Fasn", "Fbxo32", "Fbxo7", "Fdft1", "Fkbp4", "Fmn1", "Fxyd3", "Gale", "Gapdh", "Gapdhs", "Gas2l3", "Gas5", "Gas7", "Gcsh", "Gdf15", "Gjb1", "Gmnn", "Gmpr", "Gpatch4", "Gpm6b", "Gpr137b", "Gpr143", "Gps1", "Gstp1", "Gtf2f2", "Gyg2", "H2afz", "Hax1", "Hddc2", "Hist1h2ac", "Hist1h2bd", "Hist3h2a", "Hmg20b", "Hmga1", "Hps4", "Hps5", "Hsbp1", "Hsp90aa1", "Hsp90ab1", "Hspa4", "Hspa9", "Hspd1", "Hspe1", "Hsph1", "Igsf11", "Igsf3", "Igsf8", "Ilf2", "Immp2l", "Inpp5f", "Irf4", "Isyna1", "Kcnj13", "Lage3", "Ldhb", "Lhfpl3-as1", "Linc00473", "Linc00518", "Linc00673", "Loc100126784", "Loc100127888", "Loc100130370", "Loc100133445", "Loc100505865", "Loc146481", "Loc340357", "Loxl4", "Lsm2", "Lzts1", "Mad2l1bp", "Magea12", "Magea2", "Magea2b", "Magea3", "Magea4", "Magea6", "Magec1", "Maged2", "Mdh1", "Mdh2", "Mettl23", "Mettl9", "Mfi2", "Mia", "Mif", "Mitf", "Mki67ip", "Mlana", "Mlph", "Mok", "Morn2", "Mrpl12", "Mrpl21", "Mrpl23", "Mrpl24", "Mrpl38", "Mrpl40", "Mrps21", "Mrps23", "Mrps25", "Mrps26", "Mrps6", "Msi2", "Mthfs", "Mxi1", "Myo10", "Nars2", "Nav2", "Ndufa4", "Ndufaf3", "Ndufb9", "Ndufs2", "Nedd4l", "Nelfcd", "Nfya", "Ngrn", "Nhp2", "Nme1", "Nop58", "Npm1", "Nsg1", "Nt5dc3", "Nup93", "Oca2", "Pacsin2", "Pafah1b3", "Page5", "Paics", "Pax3", "Pebp1", "Peg10", "Pex19", "Pfdn2", "Phactr1", "Phb", "Phf5a", "Phlda1", "Pigy", "Pir", "Plekhb1", "Plp1", "Pmel", "Pold2", "Polr2f", "Pomgnt1", "Ppil1", "Prame", "Prdx6", "Psmb4", "Psmb7", "Psmd4", "Pttg1", "Puf60", "Pygb", "Pyurf", "Qdpr", "Qpct", "Rab17", "Rab38", "Rabggtb", "Rad51c", "Rae1", "Ran", "Rap1gap", "Rgs20", "Ropn1", "Ropn1b", "Rrs1", "Rsl1d1", "Rtkn", "Ruvbl2", "S100a1", "S100b", "Samm50", "Scd", "Sdc3", "Sdcbp", "Sdhc", "Sec11c", "Sf3a1", "Sgcd", "Shc4", "Skp1", "Slc19a1", "Slc19a2", "Slc24a5", "Slc25a13", "Slc25a4", "Slc39a4", "Slc45a2", "Slc5a3", "Slc6a15", "Slc7a5", "Slmo2", "Sms", "Snca", "Snhg16", "Snhg6", "Snrpc", "Snrpd1", "Snrpe", "Sod1", "Sord", "Sort1", "Sox10", "Spcs1", "Srp9", "St13", "St3gal4", "St3gal6", "St6galnac2", "Stam", "Stip1", "Stk32a", "Stmn1", "Stoml2", "Stra13", "Stradb", "Stx7", "Stxbp1", "Supt4h1", "Syngr1", "Tbc1d10a", "Tbc1d16", "Tbc1d7", "Tbca", "Tbrg4", "Tex2", "Tfap2a", "Timm50", "Tmed10", "Tmem147", "Tmem177", "Tmem255a", "Tmx4", "Tom1l1", "Tomm20", "Tomm22", "Tomm6", "Tomm7", "Top1mt", "Trim2", "Trim63", "Trmt112", "Tsnax", "Ttll4", "Tubb2a", "Tubb2b", "Tubb4a", "Tyr", "Tyrp1", "Uba2", "Ubl3", "Uchl5", "Uqcrh", "Utp18", "Vat1", "Vdac1", "Wbp11", "Wbp2", "Wbscr22", "Wdfy1", "Wdr43", "Xage1a", "Xage1b", "Xage1c", "Xage1d", "Xage1e", "Xylb", "Zcchc17", "Zfas1", "Zfp106", "Znf280b")
geneSets <- GeneSet(genes4, setName="Malignant")
geneSets
#cells_rankings <- AUCell_buildRankings(seu@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Malignant<-getAUC(cells_AUC)
Malignant<-t(Malignant)
seu@meta.data<-cbind(seu@meta.data, Malignant)
pdf("/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/YUMM_NRAS/UMAP_Malignant_signature.pdf", width = 7, height = 5)
FeaturePlot(seu, features = c ("Malignant"), reduction = "umap.scvi")+NoAxes()
dev.off()


#proportion of immune cells
seu$ImmuneVSrest <- ifelse(seu$Immune >0.065, "IMMUNE", "REST")
table(seu$ImmuneVSrest)
seu@meta.data$"Sample_type" <- plyr::revalue(as.character(seu$orig.ident),
                                                c("NRO001" = "YUMM52aPD1",
                                                  "NRO002" = "YUMM52aNK11",
                                                  "NRO003" = "YUMM52control",
                                                  "NRO004" = "YUMM52aPD1aNK11",
                                                  "NRO005" = "NRAScontrol",
                                                  "NRO006" = "NRASaPD1",
                                                  "NRO007" = "NRASaNK11",
                                                  "NRO008" = "NRASaPD1aNK11")) 
TEST <- seu@meta.data
cell_num <- TEST %>%
  dplyr::mutate(sample_id = as.factor(paste(`Sample_type`, ImmuneVSrest, sep = "_"))) %>%
  mutate(Sample_type = as.factor(Sample_type)) %>%
  group_by(sample_id) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::separate(sample_id, c("Sample_type", "ImmuneVSrest"))%>%
  complete(Sample_type, ImmuneVSrest, fill = list(n = 0))

cell_num
total_cells<- TEST %>%
  group_by(Sample_type) %>%
  dplyr::summarise(total = n())

total_cells
cell_percentage<- left_join(cell_num, total_cells) %>%
  mutate(percentage = n/total*100)
cell_percentage
cell_percentage$Sample_type <- factor(cell_percentage$Sample_type, 
                                      levels = c("Ctrl", "aNK", "aPD1", "Combi"))

ggplot(cell_percentage, aes(Sample_type, y=percentage, fill=ImmuneVSrest)) + geom_bar(stat="identity",  position = "stack",) +RotatedAxis() +  labs(fill = "origin") + theme(axis.title.x=element_blank()) + theme_classic() +RotatedAxis()

##################### immune subset
Immune <- subset(seu, subset = Immune >0.065) #0.05
DimPlot(Immune, reduction = "umap.scvi")
Immune <- SCTransform(Immune, verbose = TRUE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'),vst.flavor = "v2")
Immune<- RunPCA(Immune, features = VariableFeatures(object = Immune))
Immune  <- FindNeighbors(Immune , dims = 1:20,reduction = "pca")
Immune  <- RunUMAP(Immune , dims=1:20,reduction.name = "umap.unintegrated")
Immune <- FindClusters(Immune, resolution = 1.5, cluster.name = "seurat_clusters")
DimPlot(Immune)

##################### SCVI
Immune <- IntegrateLayers(object = Immune, method = scVIIntegration, new.reduction = "integrated.scvi",normalization.method = "SCT",
                          conda_env = "/Users/micromamba/envs/scvi-env", verbose = FALSE,group.by = "orig.ident")
ElbowPlot(Immune)
Immune <- FindNeighbors(Immune, reduction = "integrated.scvi", dims = 1:30)
Immune <- RunUMAP(Immune, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")
Immune <- FindClusters(Immune, resolution = 0.3, cluster.name = "scvi_clusters")
DimPlot(Immune,group.by = "Sample_type",reduction = "umap.scvi")
DimPlot(Immune,group.by = "scvi_clusters",reduction = "umap.scvi", label = TRUE)

############## marker genes
Immune <- PrepSCTFindMarkers(Immune)
Immune_Markers <- FindAllMarkers(Immune, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.3)
write.table(Immune_Markers, "/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/YUMM_NRAS/Immune_Markers.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)
top5 <- Immune_Markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("/Users/Documents/PROJECTS/NK_NRAS_INK4/Markers_heatmap_Immune.pdf", width = 25, height = 30)
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
DoHeatmap(Immune, size = 4, draw.lines = FALSE, features = top5$gene) + scale_fill_gradientn(colours = rev(mapal))+ NoLegend()# + theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
dev.off()
#saveRDS(Immune, "/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/YUMM_NRAS/subset_integrated_Immune.rds")
Immune <- readRDS("/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/YUMM_NRAS/subset_integrated_Immune.rds")

##################### singleR
library(celldex)
library(SingleR)
ref <- fetchReference("immgen", "2024-02-26")
ref_filtered <- ref[, ref$label.main != "B cells, pro"]
ref_filtered <- ref_filtered[, ref_filtered$label.main != "Stromal_cells"]
ref_filtered <- ref_filtered[, ref_filtered$label.main != "Eosinophils"]
ref_filtered <- ref_filtered[, ref_filtered$label.main != "Epithelial cells"]
ref_filtered <- ref_filtered[, ref_filtered$label.main != "Endothelial cells"]
ref_filtered <- ref_filtered[, ref_filtered$label.main != "Stem cells"]
ref_filtered <- ref_filtered[, ref_filtered$label.main != "Fibroblasts"]
ref_filtered <- ref_filtered[, ref_filtered$label.main != "Tgd"]
ref_filtered <- ref_filtered[, ref_filtered$label.main != "Microglia"]
ref_filtered <- ref_filtered[, ref_filtered$label.main != "ILC"]
table(ref$label.main)
# Extract the data from the Seurat object
data_matrix <- GetAssayData(object = Immune, assay = "SCT", slot = "data")
# Perform SingleR annotation
pred.hesc <- SingleR(test = data_matrix, ref = ref_filtered, assay.type.test = 1, labels = ref_filtered$label.main)
# Check the results
head(pred.hesc)
Immune[["SingleR.labels"]] <- pred.hesc$labels
DimPlot(Immune, group.by = "SingleR.labels", reduction = "umap.scvi")+NoAxes()
DimPlot(Immune, group.by = "SingleR.labels", reduction = "umap.scvi", split.by = "Sample_type")+NoAxes()

#asign cluster identity
Immune@meta.data$"Cell_type" <- plyr::revalue(as.character(Immune$scvi_clusters),
                                              c("0" = "MacrophagesC1qc",
                                                "1" = "Monocytes",
                                                "2" = "Cd8TcellsActivated",
                                                "3" = "MacrophagesMrc1", 
                                                "4" = "Cd8TcellsTox",
                                                "5" = "cDC2s",
                                                "6" = "Tregs",
                                                "7" = "Neutrophils",
                                                "8" = "mDCs",
                                                "9" = "cDC1s",
                                                "10" = "NKcells",
                                                "11" = "Ribosomal",
                                                "12" = "pDCs",
                                                "13" = "Cd4Tcells",
                                                "14" = "Bcells",
                                                "15" = "Granulocytes",
                                                "16" = "Bcells"
                                                ))
Immune$Sample_type <- factor(Immune$Sample_type,levels = c("NRAScontrol",
                                                  "NRASaPD1",
                                                  "NRASaNK11",
                                                  "NRASaPD1aNK11",
                                                  "YUMM52control",
                                                  "YUMM52aPD1",
                                                  "YUMM52aNK11",
                                                  "YUMM52aPD1aNK11"))
colourCount = length(unique(Immune$seurat_clusters))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
pdf("/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/YUMM_NRAS/UMAP_only_immune_cells.pdf", width = 8, height = 6)
DimPlot(Immune, group.by = "Cell_type", reduction = "umap.scvi", cols = getPalette(colourCount))+NoAxes()
dev.off()

#####################################################Proportions
Immune@meta.data$"Model" <- plyr::revalue(as.character(Immune$Sample_type),
                                             c("YUMM52aPD1" = "YUMM52",
                                               "YUMM52aNK11" = "YUMM52",
                                               "YUMM52aPD1aNK11" = "YUMM52",
                                               "YUMM52control"= "YUMM52",
                                               "NRAScontrol"= "NRAS",
                                               "NRASaPD1" = "NRAS",
                                               "NRASaNK11" = "NRAS",
                                               "NRASaPD1aNK11" = "NRAS")) 
TEST <- Immune@meta.data
cell_num <- TEST %>%
  dplyr::mutate(sample_id = as.factor(paste(`Sample_type`, Cell_type, Model, sep = "_"))) %>%
  mutate(Sample_type = as.factor(Sample_type)) %>%
  group_by(sample_id) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::separate(sample_id, c("Sample_type", "Cell_type", "Model"))%>%
  complete(Sample_type, Cell_type, fill = list(n = 0))

cell_num
total_cells<- TEST %>%
  group_by(Sample_type) %>%
  dplyr::summarise(total = n())

total_cells
cell_percentage<- left_join(cell_num, total_cells) %>%
  mutate(percentage = n/total*100)
cell_percentage
cell_percentage$Sample_type <- factor(cell_percentage$Sample_type 
                                      ,levels = c("NRAScontrol",
                                      "NRASaPD1",
                                      "NRASaNK11",
                                      "NRASaPD1aNK11",
                                      "YUMM52control",
                                      "YUMM52aPD1",
                                      "YUMM52aNK11",
                                      "YUMM52aPD1aNK11"))

pdf("/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/YUMM_NRAS/Immune_stacked.pdf", width = 5, height = 5)
ggplot(cell_percentage, aes(Sample_type, y=percentage, fill=Cell_type)) +
  geom_bar(stat="identity", position = "stack", color ="black", linewidth = 0.1) +  # Add black borders
  scale_fill_manual(values = getPalette(colourCount)) +
  theme_classic() +
  labs(fill = "origin") +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+NoLegend()+ theme(
    plot.title = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black")
  )
dev.off()




