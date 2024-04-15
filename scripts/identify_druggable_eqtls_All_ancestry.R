library(glue)
library(dplyr)
library(here)
library(readxl)
library(data.table)
library(tidyr)

setwd("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon/")
act_pqtl <- read_xlsx("scratch/actionable_protein.xlsx", sheet="ST4", skip=2) 
mdd <- fread("scratch/FE_all_ancestry_exclu_whi_jhs_23andMe_qced_rsid.txt.gz")

head(mdd)
str(mdd)

head(act_pqtl)
str(act_pqtl)

####druggable genes
library(tibble)
# Create a vector of targets
targets <- c(
  "CACNA2D1", "CACNA1H", "CACNA1C", "ESR1", "ESR2", "DRD2", "HTR1D", "CHRM3", "GABRA1", "GABRG3",
  "GABRA6", "HRH1", "HRH3", "HRH4", "GRIA1", "GRIN1", "GRIN2B", "HDAC", "HTR2A", "HTR2C",
  "COX1", "COX2", "ADRA2A", "ADRB2", "SRD5A1", "SRD5A2", "CFTR", "COMT", "CPT1C", "PDE10A",
  "MKP1", "DUSP1", "PPARG", "PPARA", "NOS2", "MPO", "MAPK8", "MAPK14", "MMP1", "MMP7", "MMP8",
  "MMP13", "TH", "ABCB1", "VDR", "KCNA10", "CYP3A43", "ABCG2", "F2R")

# Create a data frame
druggables <- tibble(targets = targets)
print(head(druggables))


##subseting for druggables targets from eqtl data
act_pqtl <- subset(act_pqtl, `gene name` %in% druggables$targets )
table(act_pqtl$`gene name`)
#33 genes
unique(act_pqtl$`gene name`)
#"CFTR"     "ABCB1"    "ESR1"     "COMT"     "HTR2A"    "MAPK8"    "VDR"      "MAPK14"   "PDE10A"   "PPARG"    "CHRM3"    "HRH4"  "MMP7"     "ESR2"     "SRD5A1"   "GABRA6"   "HTR2C"    "DRD2"     "ADRA2A"   "CACNA1C"  "CACNA2D1" "GRIA1"    "ADRB2"    "HTR1D" "TH"       "F2R"      "GABRG3"   "PPARA"    "CACNA1H"  "MMP1"     "HRH1"     "GRIN2B"   "SRD5A2"  


##creating a code column in mdd and act_pqtl data 
mddcode<- mdd %>% mutate(code=paste0(Chromosome, ":", Position))

act_pqtlcode<-act_pqtl %>% mutate(code=paste0(chr, ":", pos))


##Subsetting by code column
mdd_pqtlcode <- subset(mddcode, code %in% act_pqtlcode$code)
str(mdd_pqtlcode)


##Filtering and writing output from Subsetting by code column
mdd_pqtlcode$fdr <- p.adjust(mdd_pqtlcode$P, "fdr")

table(mdd_pqtlcode$fdr < 0.05)

mdd_pqtlcode_sig <- subset(mdd_pqtlcode, fdr < 0.05)

pqtl_cis_sigcode <- subset(act_pqtlcode, code %in% mdd_pqtlcode_sig$code)

#pqtl_cis_sigcode <- subset(pqtl_cis, code %in% mdd_pqtlcode_sig$code)

#get unique gene names and write them out

write.table(pqtl_cis_sigcode, file="C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon-Final/data/Allancestry_eQTL.txt", row.names = F, sep="\t")


write.table(unique(pqtl_cis_sigcode$`gene name`), file="C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon-Final/data/Allancestry_Cis_eQTL.txt", row=F, col=F, qu=F)
protcode <- subset(act_pqtlcode, code %in% mdd_pqtlcode_sig$code)%>% dplyr::select(chr, pos, gene=`gene name`, EA=effect_allele, NEA=other_allele, tissue)

write.table(protcode, file="C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon-Final/data/Allprotlist_chrposcode.txt", row=F, col=F, qu=F)
