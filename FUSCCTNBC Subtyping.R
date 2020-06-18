#####################################################################################################
## FUSCCTNBC Subtyping V11
## MaDing_MD2 20180828 23:59
#####################################################################################################

rm(list = ls())
FUSCCdir <- "C:/FUSCCTNBC/" # Where you put the FUSCCTNBC files, mainly for data loading process
wkdir    <- "C:/FUSCCTNBC/Analysis/Finalizing_CC/Subtyping/" # Set your working working directory
setwd(wkdir) 

library(stringr)
library(pheatmap)

load(paste(FUSCCdir, "/Data/FUSCCTNBC_MergedData_V11_181120.Rdata", sep = ""))
FUSCCTNBC_ClinicalData <- read.table("C:/FUSCCTNBC/Data/ClinicalData/FUSCCTNBC_MergedClinicalOmicsData_V11_181120.txt",
                                     header = T, stringsAsFactors = F, sep = "\t")
row.names(FUSCCTNBC_ClinicalData) <- FUSCCTNBC_ClinicalData$Project_ID 

ls()
dir("C:/FUSCCTNBC/Data/Expression")


######################################################################################################
## 1. Data Loading and preparing
######################################################################################################
Added_RNAseq   <- read.table(paste(FUSCCdir, "Data/Expression/20181226_TNBC_RNAseq_add/tnbc_fpkm_gsum_45308c127_withtype_20181226.txt", sep = ""), 
                             header = T, stringsAsFactors = F, row.names = 1, sep = "\t")
FUSCCTNBC_Added.RNAseq  <- Added_RNAseq[, -ncol(Added_RNAseq)]
colnames(FUSCCTNBC_Added.RNAseq) <- gsub("_FPKM", "", colnames(FUSCCTNBC_Added.RNAseq), fixed = T)
Protein.Coding_Genes    <- row.names(Added_RNAseq)[which(Added_RNAseq$genetype == "protein_coding")]

## QC filtering --------------------------------------------------------------------------------------
QC <- read.table(paste(FUSCCdir, "Data/Expression/20181226_TNBC_RNAseq_add/TNBC_QC_part_20181226.txt", sep = ""), 
                 header = T, stringsAsFactors = F, row.names = 1, sep = "\t")
QC.passed.Samples  <- na.omit(row.names(QC)[QC$Raw_GCContent < 60 & QC$OverallMappingRatio >= 50]) 
QC.failed.Samples  <- na.omit(row.names(QC)[QC$Raw_GCContent >= 60 | QC$OverallMappingRatio < 50])
QC.passed.Samples  <- as.character(QC.passed.Samples) ; QC.failed.Samples <- as.character(QC.failed.Samples)
Low.purity.Samples <- c("FUSCCTNBC346") # "FUSCCTNBC430", "HTA3205"

Additional.Filter  <- TRUE
if(Additional.Filter)
{
  Customed.defined.Samples <- 0.01*QC$OverallMappingRatio * QC$Trimmed_TotalRead
  names(Customed.defined.Samples) <- row.names(QC)
  sort(Customed.defined.Samples)
  Customed.defined.Samples <- 
   c("HTA3310", "FUSCCTNBC401", "FUSCCTNBC111", "FUSCCTNBC402", "FUSCCTNBC029", "FUSCCTNBC386", "FUSCCTNBC454",
     "FUSCCTNBC440" # 10M
      #,"FUSCCTNBC325", "FUSCCTNBC436", "FUSCCTNBC318", "FUSCCTNBC199" # 12M
      #, "FUSCCTNBC345" #, "FUSCCTNBC406" # 13M
      #, "FUSCCTNBC393" # 14M
      #, "FUSCCTNBC450", "FUSCCTNBC465" # 15M
      #, "FUSCCTNBC413", "FUSCCTNBC346", "FUSCCTNBC055", "FUSCCTNBC075"
     ) ## we should not keep : FUSCCTNBC345 FUSCCTNBC199
  FUSCCTNBC_Added.RNAseq <- FUSCCTNBC_Added.RNAseq[, !colnames(FUSCCTNBC_Added.RNAseq) %in% Customed.defined.Samples]
}
FUSCCTNBC_Added.RNAseq <- FUSCCTNBC_Added.RNAseq[, !colnames(FUSCCTNBC_Added.RNAseq) %in% QC.failed.Samples]
FUSCCTNBC_Added.RNAseq <- FUSCCTNBC_Added.RNAseq[, !colnames(FUSCCTNBC_Added.RNAseq) %in% Low.purity.Samples]

## Log2 transform and data alignment -------------------------------------------------------------
FUSCCTNBC_Added.RNAseq[1:5,1:5]
FUSCCTNBC_Added.RNAseq_log2 <- apply(FUSCCTNBC_Added.RNAseq, 2, function(x){log2(x+1)})
dim(FUSCCTNBC_Added.RNAseq_log2)
dim(FUSCCTNBC_RNAseqShi_log2)

Tested_Genes <- intersect( row.names(FUSCCTNBC_Added.RNAseq_log2), row.names(FUSCCTNBC_RNAseqShi_log2)) 
length(Tested_Genes)

######################################################################################################
## 2. Batch effect
######################################################################################################
##  Controled Random Gene ----------------------------------------------------------------------------
Customed_RNAseq_Batch1 <- FUSCCTNBC_RNAseqShi_log2[Tested_Genes, ]
Customed_RNAseq_Batch1 <- Customed_RNAseq_Batch1[, !str_detect(colnames(Customed_RNAseq_Batch1), fixed("PT"))]
Customed_RNAseq_Batch2 <- FUSCCTNBC_Added.RNAseq_log2[Tested_Genes, ]
GeneFilter_Data <- cbind(Customed_RNAseq_Batch1, Customed_RNAseq_Batch2)
dim(GeneFilter_Data)

## FPKM over 1000
Mean_FPKM.By.Genes <- apply(t(GeneFilter_Data), 2, function(x){mean((2^x) - 1)})
FPKM_Too.High      <- row.names(GeneFilter_Data)[which(Mean_FPKM.By.Genes > 2000)]
## FPKM = 0 in over 30% samples
Zero_FPKM.By.Genes <- apply(t(GeneFilter_Data), 2, function(x){length(which(x==0))})
FPKM_Too.Low       <- row.names(GeneFilter_Data)[which(Zero_FPKM.By.Genes > ncol(GeneFilter_Data) * 0.3)]

TEMP_Genes <- Tested_Genes[!Tested_Genes %in% FPKM_Too.Low]
#TEMP_Genes <- TEMP_Genes[!TEMP_Genes %in% FPKM_Too.High]
TEMP_Genes <- TEMP_Genes[TEMP_Genes %in% Protein.Coding_Genes] 

set.seed(0)
Random_Gene <- sample(TEMP_Genes, 2000)
BE_Data <- GeneFilter_Data[Random_Gene, ]

## PCA ------------------------------------------------------------------------------------------------
library(psych) ; library(ggbiplot) ; library(devtools)
mydata<- t(BE_Data)
wine.pca <- prcomp(mydata, scale. = FALSE)
Platform <- ifelse(str_detect( colnames(BE_Data), fixed("PT") ), "RNAseq_Previous_Normal", 
                   ifelse(colnames(BE_Data) %in% colnames(Tested_RNAseq_Batch1), "RNAseq_Previous_Tumor", "RNAseq_2018"))
pdf("BE PCA Random 2000 Coding FPKM Tumor Only 181230.pdf", width = 10, height = 5)
ggbiplot(wine.pca, obs.scale = 1, var.scale = 1, size = 2,
         groups = Platform , ellipse = F, circle = F, var.axes=F,varname.size=1) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal',
        panel.background =element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill='transparent', color='black'))
sort(wine.pca$x[,1])[names(sort(wine.pca$x[,1])) %in% colnames(Tested_RNAseq_Batch2)]
dev.off()

## Clustering-------------------------------------------------------------------------------------------
Platform <- ifelse(colnames(BE_Data) %in% colnames(Tested_RNAseq_Batch1), "RNAseq_Previous", "RNAseq_2018")
Platform <- data.frame(Platform, stringsAsFactors = F)
row.names(Platform)  <- colnames(BE_Data)
#Platform$Sample_Type <- ifelse(str_detect(row.names(Platform), fixed(".PT")), "Normal", "Tumor")
PlatformCol <- c(list(Platform = c(RNAseq_Previous="#00BFC4",RNAseq_2018 ="#F8766D") ), 
                 list(Sample_Type = c(Normal="green",Tumor ="blue")) )
library(pheatmap)
pdf("BE Heatmap Random 2000 Coding FPKM Tumor Only 181230.pdf", width = 8, height = 8) # Tumor Only  Coding Genes
pheatmap(BE_Data, fontsize_col = 0.5, fontsize_row = 0.1, annotation_col = Platform, treeheight_col = 150,
         clustering_method = "ward.D", annotation_colors = PlatformCol)
dev.off()

#####################################################################################################
## 3. Independent Clustering
#####################################################################################################
Customed_RNAseq_Batch1 <- FUSCCTNBC_RNAseqShi_log2[Tested_Genes, ]
Customed_RNAseq_Batch1 <- Customed_RNAseq_Batch1[, !str_detect(colnames(Customed_RNAseq_Batch1), fixed("PT"))]
Customed_RNAseq_Batch2 <- FUSCCTNBC_Added.RNAseq_log2[Tested_Genes, ]
GeneFilter_Data <- cbind(Customed_RNAseq_Batch1, Customed_RNAseq_Batch2)

## FPKM over 1000
Mean_FPKM.By.Genes <- apply(t(GeneFilter_Data), 2, function(x){mean((2^x) - 1)})
FPKM_Too.High      <- row.names(GeneFilter_Data)[which(Mean_FPKM.By.Genes > 2500)]
## FPKM = 0 in over 30% samples
Zero_FPKM.By.Genes <- apply(t(GeneFilter_Data), 2, function(x){length(which(x==0))})
FPKM_Too.Low       <- row.names(GeneFilter_Data)[which(Zero_FPKM.By.Genes > ncol(GeneFilter_Data) * 0.3)]

TEMP_Genes <- Tested_Genes[Tested_Genes %in% Protein.Coding_Genes]
#TEMP_Genes <- Tested_Genes[!Tested_Genes %in% FPKM_Too.High]
TEMP_Genes <- TEMP_Genes[!TEMP_Genes %in% FPKM_Too.Low]

Tested_RNAseq_Tumor  <- GeneFilter_Data[TEMP_Genes, ]
FUSCCTNBC_RNAseq.Tumor_log2 <- Tested_RNAseq_Tumor
#save(FUSCCTNBC_RNAseq.Tumor_log2, file = "FUSCCTNBC_RNAseq.Tumor_log2_Filtered.Tumor_V12_181230.Rdata")

## Filtering SD for Clustering Genes
SDs <- NULL
for(i in 1:nrow(Tested_RNAseq_Tumor))
{
  SDs <- c(SDs, sd(as.numeric(Tested_RNAseq_Tumor[i,])))
  if(i/1000 == round(i/1000)) {print(i/nrow(Tested_RNAseq_Tumor))}
}
names(SDs) <- row.names(Tested_RNAseq_Tumor)
SDs <- sort(SDs, decreasing = T)
Clustering_Genes <- names(SDs[1:2000])

## K means clustering -------------------------------------------------------------
Scale.Index <- FALSE
Tested_Exp <- Tested_RNAseq_Tumor[Clustering_Genes, ]
if(Scale.Index) {
  Tested_Exp <- t(apply(t(Tested_Exp), 2, scale))
  colnames(Tested_Exp) <- colnames(Tested_RNAseq_Tumor) } else
{
  Tested_Exp <- t(apply(t(Tested_Exp), 2, function(x){x - median(x)}))
  colnames(Tested_Exp) <- colnames(Tested_RNAseq_Tumor) }

FUSCCTNBC_ExpforClustering <- Tested_Exp
save(FUSCCTNBC_ExpforClustering, file = "D:/FUSCCTNBC_ExpforClustering_log2FPKM_181231.Rdata")

set.seed(1234)
km1 <- kmeans(t(Tested_Exp),4,iter.max = 1000,nstart = 100)
table(km1$cluster)
FUSCC_Subtype <- as.character(km1$cluster) ; names(FUSCC_Subtype) <- names(km1$cluster)
FUSCC_Subtype[FUSCC_Subtype == "3"] <- "BLIS" ; FUSCC_Subtype[FUSCC_Subtype == "2"] <- "IM"
FUSCC_Subtype[FUSCC_Subtype == "4"] <- "MES" ; FUSCC_Subtype[FUSCC_Subtype == "1"] <- "LAR"
Patient_order <- names(sort(FUSCC_Subtype))
write.csv(sort(FUSCC_Subtype), "Patient_order_RNAseq_190108.csv")

set.seed(1234)
km2 <- kmeans(Tested_Exp, 4,iter.max = 1000,nstart = 100)
table(km2$cluster)
FUSCC_Subtype_Gene <- as.character(km2$cluster) ; names(FUSCC_Subtype_Gene) <- names(km2$cluster)
FUSCC_Subtype_Gene[FUSCC_Subtype_Gene == "1"] <- "BLIS" ; FUSCC_Subtype_Gene[FUSCC_Subtype_Gene == "4"] <- "IM"
FUSCC_Subtype_Gene[FUSCC_Subtype_Gene == "2"] <- "LAR"  ; FUSCC_Subtype_Gene[FUSCC_Subtype_Gene == "3"] <- "MES"
Gene_order <- names(sort(FUSCC_Subtype_Gene))
write.csv(sort(FUSCC_Subtype_Gene), "Gene_order_RNAseq_190108.csv")

##
FUSCC_Subtype <- matrix(sort(FUSCC_Subtype), ncol=1)
rownames(FUSCC_Subtype) <- Patient_order
Tested_Res <- cbind(FUSCC_Subtype, FUSCCTNBC_ClinicalData[Patient_order , c("mRNA_Subtype", "LiuYR_Subtype", "Lehmann_Subtype")])
xtabs(~FUSCC_Subtype + mRNA_Subtype, data = Tested_Res)
xtabs(~FUSCC_Subtype + LiuYR_Subtype, data = Tested_Res) 
xtabs(~FUSCC_Subtype + Lehmann_Subtype, data = Tested_Res)
Tested_Res[c("FUSCCTNBC116","FUSCCTNBC170","FUSCCTNBC202","FUSCCTNBC313","FUSCCTNBC370"), ]
write.csv(Tested_Res,"FUSCCTNBC_Clustering.woScale_181230.csv")
"FUSCCTNBC289    FUSCCTNBC013    FUSCCTNBC113"

## Heatmap ----------------------------------------------------------------------
Plot.Index <- FALSE
if(Plot.Index) {
  library(pheatmap)
  Tested_Exp_reorder <- Tested_Exp[Gene_order, Patient_order]
  Tested_Exp_reorder[1:5,1:5]
  
  Plotted_Anno <- Tested_Res
  Plotted_Anno[is.na(Plotted_Anno)] <- "Unknown"
  AR <- unlist(Tested_Exp["AR",Patient_order])
  CD8A <-  unlist(Tested_Exp["CD8A",Patient_order])
  Plotted_Anno <- cbind(Plotted_Anno, AR, CD8A)
  
  Lehmann_Subtype <- list(
    Lehmann_Subtype = c(BL1="#E52B2C",BL2="#F1902F", IM= "#96CA89",  MSL = "#C1DFF4",
                M= "#C6B09A", LAR = "#1D4B99", UNS="light grey", Unknown = "white") )
  LiuYR_Subtype  <- list(
    LiuYR_Subtype = c(BLIS ="#E52B2C", IM= "#96CA89", LAR = "#1D4B99", MES= "#C1DFF4", Unknown = "white"))
  
  mRNA_Subtype <- list(mRNA_Subtype =c(BLIS="#EE4923", IM="#7CC243", LAR="#9180BA", MES="#3EAADF"))
  FUSCC_Subtype<- list(FUSCC_Subtype =c(BLIS="#EE4923", IM="#7CC243", LAR="#9180BA", MES="#3EAADF"))
  Plotted_Anno.Col <- c(Lehmann_Subtype, LiuYR_Subtype, mRNA_Subtype, FUSCC_Subtype)
  
  pos <- max(Tested_Exp_reorder)
  neg <- min(Tested_Exp_reorder)
  poscut <- 100
  negcut <- 100
  mybreaks1 <- -(seq(-sqrt(-neg), 0, sqrt(-neg)/negcut)) ^ 2
  mybreaks2 <- (seq(0, sqrt(pos), sqrt(pos)/poscut)[-1]) ^ 2
  mybreaks <- c(mybreaks1, mybreaks2)
  mycolors <- c(colorRampPalette(c("Green", "Black"))(negcut), colorRampPalette(c("Black", "Red"))(poscut))
  png(paste("Heatmap_woScale.", Scale.Index, "_FPKMthre500", ".png", sep = ""), width = 860, height = 640) #pdf(file = "heatmap.pdf")
  pheatmap(Tested_Exp_reorder, color = mycolors, breaks=mybreaks , cluster_rows = F, cluster_cols=F,
           show_colnames = F, show_rownames = T, fontsize_row=0.5, 
           annotation_col = Plotted_Anno, annotation_colors = Plotted_Anno.Col,
           treeheight_col = 50)
  dev.off() 
}

#####################################################################################################
## 4. Further testing in FUSCC
#####################################################################################################
Demo_Clustering <- read.csv("Demo.Clustering.woScale.csv", header = T, stringsAsFactors = F, row.names = 1) 
#Demo.Clustering.woScale.csv  #Demo.Clustering.Scaled.csv

Demo_ClinicalData <- FUSCCTNBC_ClinicalData[row.names(Demo_Clustering), ]
Demo_ClinicalData <- cbind(Demo_ClinicalData, Demo_Clustering)
Demo_ClinicalData <- Demo_ClinicalData[Demo_ClinicalData$FUSCC_Subtype  %in% c("BLIS", "IM"),]

########### New FUSCC ----------------------------------------
SurvT <-  Demo_ClinicalData$RFS_time_Months
SurvS <-  Demo_ClinicalData$RFS_Status
Tested_Var  <- Demo_ClinicalData$FUSCC_Subtype
Color_Demo  <- rainbow(length(table(Demo_ClinicalData$FUSCC_Subtype)))

#SurvS[SurvT >= 60] <- 0
#SurvT[SurvT >= 60] <- 60

library(survival)
MySurv <- Surv(SurvT, SurvS)
survdiff(MySurv ~ Tested_Var)
MySurvFit <- survfit(MySurv ~ Tested_Var)
plot(MySurvFit, col = Color_mRNA_Subtype, lwd = 4, xaxt = "n", ylim = c(0.2, 1))
abline(v = 60, col = "#CACACA", lty = 3, lwd = 2)
Customed_Labels <- 12*c(1:10)
axis(1, at = Customed_Labels, labels = as.character(Customed_Labels), cex = 2)

Demo_Clustering[Demo_ClinicalData$mRNA_Subtype == "LAR" 
                  & Demo_ClinicalData$RFS_Status == "1",]

Customed_LAR <- c("FUSCCTNBC098","FUSCCTNBC403","HTA4160","HTA3431")
Customed_Group <- Demo_Clustering$FUSCC_Subtype
names(Customed_Group) <- row.names(Demo_Clustering)
Customed_Group[names(Customed_Group) %in% Customed_LAR] <- 5

boxplot(Tested_Exp["AR", names(Customed_Group)] ~ Customed_Group)
stripchart( Tested_Exp["AR", names(Customed_Group)] ~ Customed_Group ,
            vertical = TRUE, log = "y", method = 'jitter', jitter = 0.2, cex = 1,
            pch = c(16, 16), add = T )

########### MRNA -------------------------------------------
SurvT <-  FUSCCTNBC_ClinicalData$RFS_time_Months
SurvS <-  FUSCCTNBC_ClinicalData$RFS_Status
Tested_Var  <- FUSCCTNBC_ClinicalData$mRNA_Subtype

#SurvS[SurvT >= 60] <- 0
#SurvT[SurvT >= 60] <- 60

library(survival)
MySurv <- Surv(SurvT, SurvS)
survdiff(MySurv ~ Tested_Var)
MySurvFit <- survfit(MySurv ~ Tested_Var)
plot(MySurvFit, col = Color_mRNA_Subtype, lwd = 4, xaxt = "n", ylim = c(0.2, 1))
abline(v = 60, col = "#CACACA", lty = 3, lwd = 2)
Customed_Labels <- 12*c(1:10)
axis(1, at = Customed_Labels, labels = as.character(Customed_Labels), cex = 2)

#method 11: consensus clusting
date()
#source("https://bioconductor.org/biocLite.R")
#biocLite("ConsensusClusterPlus")
library(ConsensusClusterPlus)
rcc = ConsensusClusterPlus(Tested_Exp,maxK=10,reps=1000,pItem=0.8,pFeature=1,title="consensus",distance="euclidean",clusterAlg="km",plot="png", verbose = T)
date()

###########################################################################################
BLIS_Patients <-Demo_ClinicalData$Project_ID[Demo_ClinicalData$FUSCC_Subtype == "BLIS"
                                                 #& FUSCCTNBC_ClinicalData$Chemotherapy == "TRUE"
                                                 #& Demo_ClinicalData$Lost == FALSE
                                                 #& Demo_ClinicalData$Note != "NAC"
                                                 & Demo_ClinicalData$Ascat_ACF != 1
                                                 ]

PatientFilter <- na.omit(BLIS_Patients)
BLIS_ClinicalData <- Demo_ClinicalData[BLIS_Patients, ]
##

SurvT <-  BLIS_ClinicalData$RFS_time_Months
SurvS <-  BLIS_ClinicalData$RFS_Status
Tested_Var <- BLIS_ClinicalData$HRD
#Tested_Var <- ifelse(Tested_Var >= quantile(Tested_Var, 0.5), "High", "Low")
Tested_Var <- ifelse(Tested_Var >= 42, "High", "Low")

SurvS[SurvT >= 120] <- 0
SurvT[SurvT >= 120] <- 120

library(survival)
MySurv <- Surv(SurvT, SurvS) 
survdiff(MySurv ~ Tested_Var)
MySurvFit <- survfit(MySurv ~ Tested_Var)
plot(MySurvFit, col = Color_JCO[c("Yellow", "Blue")], lwd = 4, xaxt = "n", ylim = c(0.3, 1))
Customed_Labels <- 12*c(1:10)
axis(1, at = Customed_Labels, labels = as.character(Customed_Labels), cex = 2)










###### 1.5 Heatmap --------------------------------------------------------------------------------------
if(Type == "RSEM")
{Plotted_Exp <- Tested_Exp[GenesforClustering[GenesforClustering %in% row.names(Tested_Exp)], row.names(Tested_ClinicalData_RO)]} else
{Plotted_Exp <- Tested_Exp[GenesforClustering_multiIDs$ENSEMBL[GenesforClustering_multiIDs$ENSEMBL %in% 
                                                                 row.names(Tested_Exp)],row.names(Tested_ClinicalData_RO)]}
#BackupRownames <- colnames(Plotted_Exp)
#Plotted_Exp    <- t(apply(t(Plotted_Exp), 2, scale))
#colnames(Plotted_Exp) <-BackupRownames

pos <- 8.7
neg <- -8
poscut <- 50
negcut <- 50
mybreaks1 <- -(seq(-sqrt(-neg), 0, sqrt(-neg)/negcut)) ^2
mybreaks2 <- (seq(0, sqrt(pos), sqrt(pos)/poscut)[-1]) ^2
mybreaks <- c(mybreaks1, mybreaks2)
mycolors <- c(colorRampPalette(c("green", "black"))(negcut), colorRampPalette(c( "black", "red"))(poscut))
pdf("TN.pamr.Subtype.TCGA.RSEM.pdf", width = 16, height = 16)
pheatmap(Plotted_Exp, annotation_col = Plotted_Anno, annotation_colors = Plotted_AnnoCol,
         cluster_rows = F, cluster_cols = F,  breaks = mybreaks,
         color = mycolors, fontsize_row = 5, fontsize_col = 1)
dev.off()

write.csv(Tested_ClinicalData_RO, "Tested_TCGA.ClinicalData_with.TNBCsubtype4RSEM_181224.csv")


