library("DESeq2")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("ggplot2")
library(annotate)
library(org.Hs.eg.db)
library("genefilter")
library("Biobase")
library("matrixStats")
library("IHW")
library(data.table)
library(org.Hs.eg.db)
library(GSVA)
library(dplyr)
library(ggplot2)
library(ggpubr)
fpkm_riaz <- read.csv("/Users/Documents/PROJECTS/RIAZ_transcriptome/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv", stringsAsFactors = F)
fpkm_riaz[1:5, 1:5]
fpkm_riaz[1, ] <- colnames(fpkm_riaz)
fpkm_riaz[1:5, 1:5]
fpkm_riaz_hugo <- mapIds(org.Hs.eg.db,
                         keys = as.character(fpkm_riaz$X),
                         keytype="ENTREZID",
                         column="SYMBOL", multiVals = "first")
fpkm_riaz$HUGO <- fpkm_riaz_hugo
class(fpkm_riaz$HUGO)
table(duplicated(fpkm_riaz$HUGO))
test <- fpkm_riaz[!duplicated(fpkm_riaz$HUGO),] ### removed duplicated gene names 
test[1:5, 1:5]
test <- na.omit(test)
row.names(test) <- test$HUGO
test$X <- NULL
test$HUGO <- NULL
dim(test)
test[1:5, 1:5]
test1 <- as.matrix(test)
dim(test1)
test1[1:5, 1:5]
# Convert to numeric matrix while preserving row and column names
test1 <- as.matrix(test1)
test1 <- matrix(as.numeric(test1), ncol = ncol(test1))
dimnames(test1) <- list(rownames(test), colnames(test))
gs1 <-list("Antigen_presentation",c("HLA-DRA","CD74","HLA-DPA1","HLA-DRB1","HLA-DPB1","GBP1","HLA-B","GBP2","HLA-DRB5","HLA-C","HLA-DQA1","IRF1","GBP4","HLA-DMA","WARS","TAP1","HLA-E","HLA-A","B2M","PSMB9","STAT1","HLA-DQB1","RARRES3","APOL6","HLA-F","TAPBP","PSMB8","SERPING1","CIITA","IL18BP","APOL1","HLA-DMB","TNFSF13B","UBE2L6","PARP14","XAF1","IFIH1","TAP2","TRIM22","PSME1","IFI35","LRP2","PSME2","EPSTI1","PARP9","SAMD9L","IFI6","ISG15","HAPLN3","LAP3","CTSS","IFI44L","DTX3L","GBP3","NLRC5","TNC","VAMP5","DDX60","BST2","BTN3A1","IFITM3","NMI","CASP1","BTN3A2","IFITM1","PARP12","OAS1","SAMHD1","ZNFX1","IRF7","IFIT3","OAS3","BTN3A3","IFIT2","S100A10","CTSO","GSDMD","TRIM69","CD47","IFI16","RTP4","C1R","APOL2","SP100","OAS2","STAT2","C1S","IFI44","CST3","PLSCR1","HERC6","MIA","TRIM56","PDLIM4","HELZ2","LGALS3BP","EIF2AK2","ITIH6","CDH19","MX1"))
test1<- log10(test1+1)
gsva.APC <- gsva(test1, gs1)
dim(gsva.APC)
row.names(gsva.APC) <-c("Antigen_presentation")
test_with_gsva <-rbind(test1,gsva.APC) 
test_with_gsva[1:5, 1:5]
test_with_gsva <- as.data.frame(test_with_gsva)
test_with_gsva[1:5, 1:5]
fpkm_riaz_t <- t(test_with_gsva)
fpkm_riaz_t[1:5, 1:5]
class(fpkm_riaz_t)
fpkm_riaz_t <- as.data.frame(fpkm_riaz_t)
class(fpkm_riaz_t$NAT2)
class(fpkm_riaz_t$NAT2)
Sample <- substr(rownames(fpkm_riaz_t),1,nchar(rownames(fpkm_riaz_t))-11)
fpkm_riaz_t[1:5, 1:5]
i <- colnames(fpkm_riaz_t)
fpkm_riaz_t[ , i] <- apply(fpkm_riaz_t[ , i], 2,            # Specify own function within apply
                           function(x) as.numeric(as.character(x)))
fpkm_riaz_t$Sample <- Sample

#write.table(fpkm_riaz_t,"/Users/Documents/PROJECTS/RIAZ_transcriptome/Riaz.txt", row.names = T,col.names = T,sep = "\t", quote = F)

clinical_data <- readxl::read_excel("/Users/Documents/PROJECTS/RIAZ_transcriptome/Riaz_clinical.xlsx")
dim(clinical_data)

merged <- merge(fpkm_riaz_t, clinical_data, by = "Sample") 
merged[1:5, 1:5]

merged$RESPONSE <- ifelse(merged$BOR == "CR", "Responders",
                          ifelse(merged$BOR == "PR", "Responders", "NonResponders"))
table(merged$RESPONSE)
merged$RESPONSE_1 <- ifelse(merged$BOR == "CR", "Responders",
                            ifelse(merged$BOR == "PR", "Responders",
                                   ifelse(merged$BOR == "SD", "StableDisease", "NonResponders")))
table(merged$BOR)
table(merged$PopCateg)
table(merged$RESPONSE, merged$PreOn)
merged <- merged %>% subset(Response...11 !="NE") 
#merged <- merged %>% subset(PopCateg !="ProgOn") 
#merged <- merged %>% subset(PopCateg !="ProgPre") 
table(merged$RESPONSE, merged$PreOn)
merged$var <- paste(merged$PreOn,merged$RESPONSE, sep = "_" )
table(merged$var)
my_comp <- list(c("On_NonResponders", "On_Responders"), c("Pre_NonResponders", "Pre_Responders"), c("Pre_NonResponders", "On_NonResponders"), c("Pre_Responders", "On_Responders"))

merged$var1 <- paste(merged$PreOn,merged$RESPONSE_1, sep = "_" )
table(merged$var1)
my_comp_1 <- list(c("Pre_StableDisease", "On_StableDisease"), c("Pre_NonResponders", "On_NonResponders"), c("Pre_Responders", "On_Responders"))


########################## CIBERSORT
library(data.table)
ciber_output <- fread("/Users/Documents/PROJECTS/Grand_Challenge/CibersortX/CIBERSORTx_Job2_output/CIBERSORTxGEP_Job2_Fractions-Adjusted.txt")
ciber_output <- as.data.frame(ciber_output)
ciber_output$Mixture<- substr(ciber_output$Mixture, 1 ,nchar(ciber_output$Mixture)-11)
ciber_output$Sample <- ciber_output$Mixture
clinical_data <- readxl::read_excel("/Users/Documents/PROJECTS/RIAZ_transcriptome/Riaz_clinical.xlsx")
dim(clinical_data)
class(clinical_data)
merged_1 <- merge(ciber_output, clinical_data, by = "Sample") 
merged_1[1:5, 1:10]
merged_1$RESPONSE <- ifelse(merged_1$BOR == "CR", "Responders",
                            ifelse(merged_1$BOR == "PR", "Responders","NonResponders"))
table(merged_1$RESPONSE)
merged_1$RESPONSE_1 <- ifelse(merged_1$BOR == "CR", "Responders",
                              ifelse(merged_1$BOR == "PR", "Responders",
                                     ifelse(merged_1$BOR == "SD", "StableDisease", "NonResponders")))

merged_1 <- merged_1 %>% subset(Response...11 !="NE")
table(merged_1$RESPONSE, merged_1$PreOn)
merged_1$var<- paste(merged_1$PreOn,merged_1$RESPONSE, sep = "_" )
table(merged_1$var)
merged_1$var <- factor(merged_1$var, levels = c("Pre_NonResponders", "On_NonResponders","Pre_Responders","On_Responders"))
merged_1 <- merged_1 %>%
  mutate(var = recode(var , "Pre_NonResponders" = "Pre_NR", "On_NonResponders"="On_NR","Pre_Responders"="Pre_R","On_Responders" = "On_R"))
my_comp <- list(c("On_NR", "On_R"), c("Pre_NR", "Pre_R"), c("Pre_NR", "On_NR"), c("Pre_R", "On_R"))
pdf("/Users/Documents/PROJECTS/RIAZ_transcriptome/Results/NK_RvsNR_Timepoint.pdf", height = 4, width = 3)
ggboxplot(merged_1, x = "var", y = "NK",
          fill = "PreOn",
          shape = "RESPONSE",
          palette = c("tan","#0969A2"),
          add = "jitter")+
  stat_compare_means(comparisons =  my_comp, label = "p.format", method =  "wilcox.test")  + theme_classic() + RotatedAxis()+NoLegend()
dev.off()


