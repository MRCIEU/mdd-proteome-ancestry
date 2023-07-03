library(TwoSampleMR)
library(dplyr)
#install.packages("here")
library(here)
library(ggplot2)
library(png)

# Read in exposure data

pqtl <- readRDS(here("data", "ancestry_pqtl.rds"))
#pqtl <- readRDS(file.path("data", "ancestry_pqtl.rds"))
class(pqtl)


lapply(pqtl, class)

lapply(pqtl, dim)

# only keep eas

pqtl <- pqtl$EAS
head(pqtl)

# Get the p-value
pqtl$pval <- 10^-pqtl$LOG10P
head(pqtl)

# How many SNPs for one protein?
subset(pqtl, prot == "AKT3")

# How many snps for all proteins?
table(pqtl$prot)

# pqtl data doesn't have rsid - get that from combined_pqtl

#pqtl_combined <- readRDS(file.path("data", "combined_pqtl.rds"))
pqtl_combined <- readRDS(here("data", "combined_pqtl.rds"))

head(pqtl_combined)

# add rsid to pqtl
pqtl <- left_join(
  pqtl,
  subset(pqtl_combined, select=c(CHROM, GENPOS, rsid)),
  by=c("CHROM", "GENPOS")
)

head(pqtl)
dim(pqtl)

# Get the outcome data

file.exists(file.path("data", "mdd_extract.rdata"))
load(file.path("data", "mdd_extract.rdata"))

# Look at eur MDD GWAS
head(eas)
dim(eas)

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

# check
head(exp_dat)


# Format outcome data
head(eas)
out_dat <- format_data(eas,
                       type="outcome",
                       snp_col="MarkerName",
                       effect_allele_col="Allele1",
                       other_allele_col="Allele2",
                       eaf_col="Freq1",
                       beta_col="Effect",
                       se_col="StdErr",
                       pval_col="P.SE"
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

# Perform MR

res <- mr(dat)
res

ggplot(res, aes(y=exposure, x=b)) +
  geom_point() +
  geom_errorbarh(aes(xmin=b-se*1.96, xmax=b+se*1.96)) +
  geom_vline(xintercept=0)

save(res, file=file.path("results", "eas.rdata"))
ggsave(file=file.path("results", "easplot.png"))

