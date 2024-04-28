library(TwoSampleMR)
library(dplyr)
library(here)
library(ggplot2)

# Read in exposure data

pqtl <- readRDS(file.path("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon-2023-6/data", "ancestry_pqtl.rds"))
class(pqtl)


lapply(pqtl, class)

lapply(pqtl, dim)

# only keep european
pqtl <- pqtl$EUR
head(pqtl)
unique(pqtl$prot)
#"AKT3"   "BTN2A1" "BTN3A2" "CCS"    "CD40"   "CDH6"   "DPT"    "DSG3"   "FES"    "HLA-E"  "ISLR2"  "ITIH3"  "LTBP2"  "SPRY2" "SSC4D"  "TIMP4"  "TNR"   
# Get the p-value
pqtl$pval <- 10^-pqtl$LOG10P
head(pqtl)
pqtl_sig<-pqtl[pqtl$pval < 5e-8,]

# How many SNPs for one protein?
subset(pqtl, prot == "AKT3")

# How many snps for all proteins?
table(pqtl$prot)

# pqtl data doesn't have rsid - get that from combined_pqtl

pqtl_combined <- readRDS(file.path("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon-2023-6/data", "combined_pqtl.rds"))

head(pqtl_combined)

# add rsid to pqtl
pqtl <- left_join(
  pqtl,
  subset(pqtl_combined, select=c(CHROM, GENPOS, rsid)),
  by=c("CHROM", "GENPOS")
)


head(pqtl)

# Get the outcome data

#load(file.path("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon-2023-6/data", "mdd_extract.rdata"))

load(file.path("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon/data", "mdd_extract.rdata"))
# Look at eur MDD GWAS
head(eur_pos)

dim(eur_pos)

# TwoSampleMR
# Format data
# Harmonise
# MR


# Format data
exp_dat <- format_data(pqtl,
                       type="exposure",
                       snp_col="rsid",
                       phenotype_col = "prot",
                       beta_col = "BETA",
                       se_col = "SE",
                       eaf_col = "A1FREQ",
                       effect_allele_col = "ALLELE1",
                       other_allele_col = "ALLELE0",
                       pval_col = "pval",
                       chr_col = "CHROM",
                       pos_col = "GENPOS"
)


head(exp_dat)


# Format outcome data
out_dat <- format_data(eur_pos,
    type="outcome",
    snp_col="MarkerName",
    effect_allele_col="A1",
    other_allele_col="A2",
    eaf_col="Freq",
    beta_col="LogOR",
    se_col="StdErrLogOR",
    pval_col="P"
)

out_dat$outcome <- "MDD"
head(out_dat)

 
# Harmonise
# Assume forward strand
dat <- harmonise_data(exp_dat, out_dat, action=1)

head(dat)

# Keep just best SNP for each protein in the harmonised data

dim(dat)

dat <- dat %>%
    group_by(id.exposure) %>%
    arrange(pval.exposure, desc(abs(beta.exposure))) %>%
    slice_head(n=1)
dim(dat)


write.csv(dat, "C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon-Final/results/EUR_SNP.csv", row.names = FALSE)


# Perform MR

res <- mr(dat)
res

ggplot(res, aes(y=exposure, x=b)) +
geom_point() +
geom_errorbarh(aes(xmin=b-se*1.96, xmax=b+se*1.96)) +
geom_vline(xintercept=0)

ggsave(file=file.path("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon-Final/results", "eur.pdf"))
save(res, file=file.path("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon-Final/results", "eur.rdata"))
