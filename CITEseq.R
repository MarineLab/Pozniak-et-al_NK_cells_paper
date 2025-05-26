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

setwd("/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/NRAS_INK4A_CD45_RNA_CITE_seq/Data")
dir <- "/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/NRAS_INK4A_CD45_RNA_CITE_seq/Data"
list.files(dir) -> samples

seurat_list <- list()
# samp <- samples[1]

for(samp in samples){
  data_10x <- Read10X_h5(filename = paste(dir,samp,"/filtered_feature_bc_matrix.h5",sep = "/"))
  rownames(x = data_10x[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
                                                   x = rownames(x = data_10x[["Antibody Capture"]]))
  
  seu1 <- CreateSeuratObject(data_10x[["Gene Expression"]],min.cells = 10, min.features = 200,project = samp)
  if("Antibody Capture" %in% names(data_10x)){
    protein_data <- data_10x$`Antibody Capture`
    
    # Ensure consistent cell barcodes between RNA and protein data
    common_cells <- intersect(colnames(seu1), colnames(protein_data))
    protein_data <- protein_data[, common_cells]
    seu1 <- subset(seu1, cells = common_cells)
    
    seu1[["ADT"]] <- CreateAssayObject(counts = protein_data)
  }
  rm(data_10x, protein_data)
  
  seu1[["percent.mt"]] <- PercentageFeatureSet(seu1, pattern = "^mt-")
  seu1<- SCTransform(seu1, verbose = TRUE, vars.to.regress = c("percent.mt"))
  
  #Generate sce object
  as.SingleCellExperiment(seu1,assay = "RNA") -> sce
  #Run doublet finder witn default parameters
  scDblFinder(sce) -> sce
  #Save classification in seurat object
  colData(sce) %>% data.frame() %>% .$scDblFinder.class -> seu1@meta.data$Doublets
  seurat_list[[samp]] <- seu1
  rm(seu1,sce)
}

#Merge and add metadata
merged_seu1 <- merge(seurat_list$NRO009,y = c(seurat_list$NRO010),project="NK_mouse_CITEseq")
merged_seu1 <- subset(merged_seu1, subset = Doublets == "singlet")

#Rename seurat clusters based on manual curation of clusters
merged_seu1@meta.data$"Sample_type" <- plyr::revalue(as.character(merged_seu1$orig.ident),
                                                    c("NRO009" = "NRAS_control",
                                                      "NRO010" = "NRAS_aPD1")) 
Idents(merged_seu1) <- merged_seu1@meta.data$Sample_type

#pdf("plots/QC.pdf", width = 12, height = 7)
VlnPlot(merged_seu1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) 
#dev.off()


#Subset by quality
seu_sub1 <- subset(merged_seu1, subset = nFeature_RNA > 500 &  nFeature_RNA < 6000 &  percent.mt < 10)
#pdf("plots/QC_subsets.pdf", width = 12, height = 7)
VlnPlot(seu_sub1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) 
#dev.off()

saveRDS(merged_seu1,"merged.rds")
saveRDS(seu_sub1,"subset.rds")
seu = seu_sub1
mouse_cell_cycle_genes <- readRDS("/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/mouse_cell_cycle_genes.rds")
s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seu <- SCTransform(seu, verbose = TRUE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'),vst.flavor = "v2")
seu <- RunPCA(seu , verbose = TRUE,npcs = 100)
ElbowPlot(seu,ndims = 50)
seu  <- FindNeighbors(seu , dims = 1:20,reduction = "pca")
seu  <- RunUMAP(seu , dims=1:20,reduction.name = "umap.unintegrated")
seu <- FindClusters(seu, resolution = 1, cluster.name = "seurat_clusters")
DimPlot(seu,reduction = "umap.unintegrated", group.by = "seurat_clusters")
VlnPlot(seu, features = c("Mlana"))
seu <- NormalizeData(seu, assay = "ADT", normalization.method = "CLR", margin = 2)
seu <- ScaleData(seu, assay = "ADT")

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
FeaturePlot(seu, features = c ("Immune"))

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
FeaturePlot(seu, features = c ("Malignant"), reduction = "umap.unintegrated")


##################### immune subset
Immune <- subset(seu, subset = Immune >0.056) #0.05
DimPlot(Immune, reduction = "umap.unintegrated", group.by = "Sample_type")
Immune <- SCTransform(Immune, verbose = TRUE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'),vst.flavor = "v2")
Immune<- RunPCA(Immune, features = VariableFeatures(object = Immune))
Immune  <- FindNeighbors(Immune , dims = 1:10,reduction = "pca")
Immune  <- RunUMAP(Immune , dims=1:10,reduction.name = "umap.unintegrated")
Immune <- FindClusters(Immune, resolution = 1, cluster.name = "seurat_clusters")
DimPlot(Immune, label= TRUE)
DimPlot(Immune, split.by = "Sample_type", label = TRUE)

DefaultAssay(Immune) <- "SCT"
Immune <- PrepSCTFindMarkers(Immune)
Immune_Markers <- FindAllMarkers(Immune, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
write.table(Immune_Markers, "/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/NRAS_INK4A_CD45_RNA_CITE_seq/Immune_Markers_CITE-seq.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)
#saveRDS(Immune,"Immune.rds")
Immune <- readRDS("/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/NRAS_INK4A_CD45_RNA_CITE_seq/Immune.rds")
colourCount = length(unique(Immune$seurat_clusters))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
pdf("/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/NRAS_INK4A_CD45_RNA_CITE_seq/UMAP_only_immune_cells.pdf", width = 5, height = 5)
DimPlot(Immune, cols = getPalette(colourCount), label = TRUE)+NoAxes()+NoLegend()
dev.off()

###################### visualize proteins
Immune <- NormalizeData(Immune, assay = "ADT", normalization.method = "CLR", margin = 2)
Immune <- ScaleData(Immune, assay = "ADT")
VlnPlot(Immune, features = "adt_CX3CR1", pt.size = 0)  + stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95) + NoLegend()

################### SUBSET NK CELLS and test DIFFERTIAL EXPRESSION
NK_only <- Immune %>% subset(seurat_clusters == "11")
DimPlot(NK_only)
NK_only <- SCTransform(NK_only, verbose = TRUE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'),vst.flavor = "v2")
NK_only<- RunPCA(NK_only, features = VariableFeatures(object = NK_only))
ElbowPlot(NK_only)
NK_only  <- FindNeighbors(NK_only , dims = 1:10,reduction = "pca")
NK_only  <- RunUMAP(NK_only , dims=1:10,reduction.name = "umap.unintegrated")
NK_only <- FindClusters(NK_only, resolution = 0.5, cluster.name = "seurat_clusters")
DimPlot(NK_only)
DimPlot(NK_only, group.by = "Sample_type")
NK_only <- SetIdent(NK_only, value =NK_only$Sample_type)
adt_markers <- FindMarkers(NK_only, ident.1 ="NRAS_control", assay = "ADT")

library(EnhancedVolcano)
library(ggrepel)
adt_markers$gene_name <- rownames(adt_markers)
# Ensure custom_colors vector matches the number of rows in adt_markers and is named
custom_colors <- ifelse(adt_markers$avg_log2FC > 0 & adt_markers$p_val_adj<10e-3, '#7CAE00',
                        ifelse(adt_markers$avg_log2FC < 0 & adt_markers$p_val_adj<10e-3, '#F8766D', 'grey'))
names(custom_colors)[custom_colors == '#7CAE00'] <- "NRASCotrol"
names(custom_colors)[custom_colors == '#F8766D'] <- 'NRAS_aPD1'

# Plot using EnhancedVolcano with custom colors
pdf("/Users/Library/CloudStorage/OneDrive-KULeuven/Documents/PROJECTS/NK_NRAS_INK4/NRAS_INK4A_CD45_RNA_CITE_seq/Volcano_plot.pdf", width = 6, height = 6)
EnhancedVolcano(
  adt_markers,
  lab = adt_markers$gene_name,
  selectLab = rownames(adt_markers)[which(names(custom_colors) %in% c('NRASCotrol', 'NRAS_aPD1'))],
  x = 'avg_log2FC',
  y = 'p_val_adj',
  pCutoff = 10e-3,
  FCcutoff = 0,
  pointSize = 2.0,            # Adjust point size
  labSize = 3.0,              # Adjust label size
  drawConnectors = TRUE,      # Draw connectors
  widthConnectors = 0.4,      # Adjust connector width
  max.overlaps = Inf,         # Allow all overlaps
  title = NULL,               # No title
  colCustom = custom_colors,
  colAlpha = 1
)
dev.off()


