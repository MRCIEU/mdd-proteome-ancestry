## 
library(data.table)
library(readxl)
library(dplyr)

mdd_targets <- read_xlsx("C:/Users/Mutema/Downloads/41591_2021_1310_MOESM3_ESM (1).xlsx",sheet= 'ST4', skip = 2 )

#select the 18 genes from the list
mdd <-c("GABRA1","GRIN1","CACNA1H","ESR1","ADRA2A","VDR","DRD2","SRD5A1",
        "CHRM3","COMT","HTR1D","HTR2A","HRH1","PPARG","GRIA1","CACNA2D1",
        "CACNA1C","PDE10A")
fileterd_mdd_targets <-mdd_targets %>% filter(`gene name`%in% mdd)

#read in the db SNP reference
ref <- fread("./data/dbsnp.v153.b37.vcf")
