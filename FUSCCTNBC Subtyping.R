#################################################################################################################################################
# FUSCCTNBC Clustering
#################################################################################################################################################

setwd("C:/FUSCCTNBC/Analysis/Finalizing_CC/Subtyping")
load("C:/FUSCCTNBC/Analysis/Finalizing_CC/Finalizing_Data/FUSCCTNBC_RNAseqShi.Tumor.Coding_log2_Subtyping.Version_360x15766_V12_190101.Rdata")

## Filtering SD for Clustering Genes
SDs <- NULL
for(i in 1:nrow(FUSCCTNBC_RNAseqShi.Tumor.Coding_log2))
{
  SDs <- c(SDs, sd(as.numeric(FUSCCTNBC_RNAseqShi.Tumor.Coding_log2[i,])))
  if(i/1000 == round(i/1000)) {print(i/nrow(FUSCCTNBC_RNAseqShi.Tumor.Coding_log2))}
}

names(SDs)       <- row.names(FUSCCTNBC_RNAseqShi.Tumor.Coding_log2)
SDs              <- sort(SDs, decreasing = T)
Clustering_Genes <- names(SDs[1:2000])

Tes_Exp <- FUSCCTNBC_RNAseqShi.Tumor.Coding_log2[Clustering_Genes, ]


## ConsensusClusterPlus
date()
#source("https://bioconductor.org/biocLite.R")
#biocLite("ConsensusClusterPlus")
library(ConsensusClusterPlus)
rcc = ConsensusClusterPlus(Tes_Exp, maxK=10,reps=1000,pItem=0.8,pFeature=1, title="FUSCCTNBC_Consensus",distance="euclidean",
                           clusterAlg="km",plot="png", verbose = T)
date()

## Sample Clustering
set.seed(1234)
km1 <- kmeans(t(Tes_Exp),4,iter.max = 1000,nstart = 100)
table(km1$cluster)
FUSCC_Subtype <- as.character(km1$cluster) ; names(FUSCC_Subtype) <- names(km1$cluster)
FUSCC_Subtype[FUSCC_Subtype == "1"] <- "BLIS" ; FUSCC_Subtype[FUSCC_Subtype == "2"] <- "IM"
FUSCC_Subtype[FUSCC_Subtype == "3"] <- "MES"  ; FUSCC_Subtype[FUSCC_Subtype == "4"] <- "LAR"
Patient_order <- names(sort(FUSCC_Subtype))
write.csv(sort(FUSCC_Subtype), "Patient_order_RNAseq_190208.csv")




