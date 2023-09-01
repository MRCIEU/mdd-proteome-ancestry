library(data.table)
library(dplyr)
#library(gwasvcf)
library(here)
library(tidyr)


# Reading and subsetting for only the 11 proteins 
prot1 <- readRDS(file.path("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon/data", "combined_pqtl.rds"))
MDDprot <- read.table(file.path("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon/data", "Allprotlist_chrposcode.txt"))
names(MDDprot) <- c("chr", "pos", "rsid", "prot")
prot <- subset(prot1, prot %in% MDDprot$prot)
chrpos <- paste0(prot$chr, ":", prot$pos)

# AFR
afr <- fread(file.path("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon/scratch", "FE_afr_exclu_whi_jhs_23andMe_qced_rsid.txt.gz"))
afr <- subset(afr, rsid %in% prot$rsid)
dim(afr)

# SAS
sas <- fread(file.path("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon/scratch", "FE_sas_exclu_whi_jhs_23andMe_qced_rsid.txt.gz"))
sas <- subset(sas, rsid %in% prot$rsid)
dim(sas)


# EAS
eas <- fread(file.path("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon/scratch", "FE_eas_exclu_whi_jhs_23andMe_qced_rsid.txt.gz"))
eas <- subset(eas, rsid %in% prot$rsid)
dim(eas)

#HIS
his <- fread(file.path("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon/scratch", "FE_his_exclu_whi_jhs_23andMe_qced_rsid.txt.gz"))
his <- subset(eas, rsid %in% prot$rsid)
dim(his)

# EUR
eur <- fread(file.path("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon/scratch", "PGC_UKB_depression_genome-wide.txt"))
eur <- subset(eur, MarkerName %in% prot$rsid)
dim(eur)

save(afr, sas, eas, his, eur, file=file.path("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon/data", "mdd_extract.rdata"))

