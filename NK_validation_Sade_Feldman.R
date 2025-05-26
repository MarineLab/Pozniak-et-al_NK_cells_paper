library(harmony)
library(data.table)
library(circlize)
library(Seurat)
library(AUCell)
library(GSEABase)
library(dplyr)
#read count data
df<-read.table("/Users/u0128760/Documents/PUBLIC_DATA/scRNAseq Sade-Feldman/tpm_all_scp.txt", header = TRUE, row.names = 1, fill = TRUE)
dim(df)
df[1:15, 16290:16291]
SF <- CreateSeuratObject(
  df,
  project = "SeuratProject",
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = NULL,
  min.cells = 10,
  min.features = 1000
)
#read meta data
metadata_SF <- fread("/Users/u0128760/Documents/PUBLIC_DATA/scRNAseq Sade-Feldman/cells_all_scp.txt")
metadata_SF <- metadata_SF[-1,]

SF@meta.data$x_merge <- rownames(SF@meta.data)
metadata_SF$x_merge <- metadata_SF$NAME
SF@meta.data <- SF@meta.data %>% inner_join(metadata_SF, by="x_merge")
rownames(SF@meta.data) <- SF@meta.data$x_merge
SF$patient_new <- gsub("_", "", SF@meta.data$patient)
SF$patient_new
table(SF$therapy)

SF@meta.data$"Therapy" <- plyr::revalue(as.character(SF$therapy),
                                        c(`CTLA4 (baseline) ; PD1 (post I)` = "CTLAbasePD1posti",
                                          `CTLA4 (baseline); PD1 (post I and II)` = "CTLAbasePD1postiii",
                                          "DNDP" = "DNDP",
                                          "PD1" = "PD1", 
                                          "CTLA4+PD1" = "CTLA4PD1"
                                        ))

SF <- SF %>% subset(Therapy %in% c("PD1", "CTLA4PD1"))
table(SF$patient, SF$response)

# from DGE (marker genes) from CD8T and NK clusters
genes<- c("TYROBP","FCER1G","TRDC","KLRB1","KLRF1","NCAM1","GNLY","XCL1","PLAC8","SH2D1B","PTGDR","XCL2","CD300A","MATK","FGR","CXXC5","C1orf162","LAT2","CEBPD","HOPX","TXK","KLRC1","AREG","FOS","CTSW","KLRD1","IFITM3","IER2","MYBL1","GADD45B","IFITM2","TNFRSF18","SELL","NFKBIA","PTPN12","EFHD2","XBP1","RPL3","NR4A1","MAP3K8","FCGR3A","CMC1","IL2RB","EEF1B2","RPL10","FOSB","PLEK","JAK1","GSTP1","IFITM1","KLF2","CD247","NFKBID","TPT1","TXNIP","RPS23","BTG2","CD7","EEF1A1","REL","NKG7","IL7R","BTG1","DUSP1","PPP1R15A","CD69","TAGLN2","CLIC3","GZMB","AC020916.1","LTB","JUNB","EIF3G","PLCG2","ZFP36","PRF1")
geneSets <- GeneSet(genes, setName="NK_1")
geneSets
cells_rankings <- AUCell_buildRankings(SF@assays$RNA@counts,splitByBlocks=TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=5, assign=TRUE)
NK_1<-getAUC(cells_AUC)
NK_1<-t(NK_1)
SF@meta.data<-cbind(SF@meta.data, NK_1)

SF$NK_ON_OFF <- ifelse(SF$NK_1 >0.2 & SF@assays$RNA@data["CD8A",]<0.00001 & SF$cluster_all =="G08" , "NK", "NON-NK") #chosen on the first bump on the histogram 0.2
table(SF$NK_ON_OFF, SF$cluster_all)

VlnPlot(SF, features = c("CD3E","CD8A", "XCL1", "XCL2", "PRF1", "GNLY", "GZMB"), group.by  = "NK_ON_OFF")

#Calculate percentage of NK cells from all immune cells
TEST <- SF@meta.data
cell_num <- TEST %>%
  mutate(sample_id = as.factor(paste(patient_new, response, prepost, Therapy, sep="_"))) %>%
  mutate(NK_ON_OFF = as.factor(NK_ON_OFF)) %>%
  group_by(sample_id, NK_ON_OFF, .drop=FALSE) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::separate(sample_id, c("patient_new", "response", "prepost", "therapy"))
cell_num
TEST1 <-  SF@meta.data
total_cells<- TEST1 %>%
  group_by(patient_new) %>%
  dplyr::summarise(total = n())
total_cells
cell_percentage<- left_join(cell_num, total_cells) %>%
  mutate(percentage = n/total*100)
cell_percentage
cell_percentage$response <- factor(cell_percentage$response, level=c("R", "NR")) 

################### overall
cell_percentage <- subset(cell_percentage, response != "DN")
cell_percentage <- subset(cell_percentage, patient_new != "PostP28")
cell_percentage <- subset(cell_percentage, NK_ON_OFF == "NK")

ggboxplot(cell_percentage, x = "response", y = "percentage",
          color = "response",
          shape = "response",
          add = "jitter", level=c("R", "NR"))+
  stat_compare_means(aes(group = response, label = paste0("p = ", ..p.format..), method =  "wilcox.test"), label.x = c(1.5)) +
  facet_wrap(~NK_ON_OFF, ncol = 13) + theme_classic()  + RotatedAxis()
#################### treatment dynamics
cell_percentage$var <- paste(cell_percentage$prepost,cell_percentage$response, sep = "_" )
cell_percentage$var <- factor(cell_percentage$var, level=c("Pre_NR", "Post_NR", "Pre_R", "Post_R")) 
table(cell_percentage$var)
my_comp <- list(c("Post_NR", "Post_R"), c("Pre_NR", "Pre_R"), c("Pre_R", "Post_R"), c("Pre_NR", "Post_NR"))
pdf("/Users/u0128760/Documents/PROJECTS/Grand_Challenge/Public_data/SADE_FELDMAN_data_analysis/NK_score_RvsNRES_timepoint.pdf", width = 3, height = 4)
ggboxplot(cell_percentage, x = "var", y = "percentage",
          fill = "prepost",
          shape = "response",
          palette = c("#E7B800","#0969A2"),
          add = "jitter")+
  stat_compare_means(comparisons =  my_comp, label = "p.format", method =  "wilcox.test") + theme_classic() + RotatedAxis() +NoLegend()
dev.off()


