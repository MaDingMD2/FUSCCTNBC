setwd("C:/00_wkdir/JNCI/Github Files")

library(stringr)

Clinvar_RS <- read.csv("Clinvar_RS.csv", header = T, stringsAsFactors = F)

## Load CharGer results files
CharGer.GATK <- read.table("TNBC34_GATK.BamReadcount.MAFfiltered.VEP.CharGer.tsv" ,sep = "\t", 
                           header = T , stringsAsFactors = F, fill = T, quote = "")
# CharGer.Del <- read.table("TNBC34_Pindel.del_BamReadcount.MAFfiltered.VEP.CharGer.tsv" ,sep = "\t", 
#                            header = T , stringsAsFactors = F, fill = T, quote = "")
# CharGer.Ins <- read.table("TNBC34_Pindel.ins_BamReadcount.MAFfiltered.VEP.CharGer.tsv" ,sep = "\t", 
#                           header = T , stringsAsFactors = F, fill = T, quote = "")

## Load vcf files
GATK_vcf <- read.table("TNBC34_GATK.BamReadcount.MAFfiltered.VEP.vcf" ,sep = "\t", 
                           header = F , stringsAsFactors = F, fill = T, quote = "")
# Del_vcf  <- read.table("TNBC34_Pindel.del_BamReadcount.MAFfiltered.VEP.vcf" ,sep = "\t", 
#                          header = F , stringsAsFactors = F, fill = T, quote = "")
# Ins_vcf  <- read.table("TNBC34_Pindel.ins_BamReadcount.MAFfiltered.VEP.vcf" ,sep = "\t", 
#                          header = F , stringsAsFactors = F, fill = T, quote = "")

## Add Clinvar to vcf
GATK_vcf_RS <- GATK_vcf[GATK_vcf$V3 != ".", ]
GATK_vcf_RS$CLNSIG <- rep("", nrow(GATK_vcf_RS))
GATK_vcf_RS$CLNDN<-  rep("",nrow(GATK_vcf_RS) )

for(i in 1: nrow(GATK_vcf_RS))
{
  if(GATK_vcf_RS[i, 3] %in% Clinvar_RS$RS)
  {
    TEMP_Index  <- which(Clinvar_RS$RS == GATK_vcf_RS[i, 3])
    GATK_vcf_RS$CLNSIG[i] <- Clinvar_RS$CLNSIG[TEMP_Index]
    GATK_vcf_RS$CLNDN[i] <- Clinvar_RS$CLNDN[TEMP_Index] }
  }
  
## Add Clinvar to CharGer output
GATK_vcf_RS$Var_ID  <- paste(GATK_vcf_RS$V1, ":", GATK_vcf_RS$V2, sep = "")
row.names(GATK_vcf_RS) <- GATK_vcf_RS$Var_ID 
CharGer.GATK$Var_ID <- paste(CharGer.GATK$Chromosome, ":", CharGer.GATK$Start, ":",
                             CharGer.GATK$Reference , ":", CharGer.GATK$Alternate, sep = "")
row.names(CharGer.GATK) <- CharGer.GATK$Var_ID 

CharGer.GATK$RS_MD     <- rep("",nrow(CharGer.GATK) )
CharGer.GATK$CLNSIG_MD <- rep("",nrow(CharGer.GATK) )
CharGer.GATK$CLNDN_MD  <- rep("",nrow(CharGer.GATK) )
  
CharGer.GATK$RS        <- GATK_vcf_RS[row.names(CharGer.GATK), "RS"]
CharGer.GATK$CLNSIG    <- GATK_vcf_RS[row.names(CharGer.GATK), "CLNSIG"]
CharGer.GATK$CLNDN     <- GATK_vcf_RS[row.names(CharGer.GATK), "CLNDN"]


for(i in 1: nrow(GATK_vcf_RS))
{
  if(GATK_vcf_RS[i, 3] %in% Clinvar_RS$RS)
  {
    TEMP_Index  <- which(Clinvar_RS$RS == GATK_vcf_RS[i, 3])
    GATK_vcf_RS$CLNSIG[i] <- Clinvar_RS$CLNSIG[TEMP_Index]
    GATK_vcf_RS$CLNDN[i] <- Clinvar_RS$CLNDN[TEMP_Index] }
}


write.csv(CharGer.GATK, "CharGer_GATK_UpdateClinvar.csv", row.names = F)

