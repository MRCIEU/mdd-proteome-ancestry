library(data.table)
library(dplyr)
library(gwasvcf)
library(here)
library(tidyr)

set_bcftools()

prot <- readRDS(here("data", "combined_pqtl.rds"))
chrpos <- paste0(prot$chr, ":", prot$pos)

# AFR
afr <- query_gwas(here("scratch", "ukb-e-20126_p1_AFR.vcf.gz"), chrompos=chrpos)
afr <- afr %>% vcf_to_tibble() %>%
    mutate(pop="AFR")

# CSA
csa <- query_gwas(here("scratch", "ukb-e-20126_p5_CSA.vcf.gz"), chrompos=chrpos)
csa <- csa %>% vcf_to_tibble() %>%
    mutate(pop="AFR")
csa

# EAS
eas <- fread(here("scratch", "mdd_eas.txt.gz"))
eas <- subset(eas, MarkerName %in% prot$rsid)
dim(eas)

# EUR
eur <- fread(here("scratch", "mdd_eur.txt.gz"))
eur <- subset(eur, MarkerName %in% prot$rsid)
dim(eur)

save(afr, csa, eas, eur, file=here("data", "mdd_extract.rdata"))

