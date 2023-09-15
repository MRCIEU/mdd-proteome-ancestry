library(TwoSampleMR)

library(dplyr)
library(here)
library(ggplot2)

# Read in exposure data

pqtl <- readRDS(here("data", "ancestry_pqtl.rds"))
pqtl <- readRDS(file.path("data", "ancestry_pqtl.rds"))
class(pqtl)


lapply(pqtl, class)

lapply(pqtl, dim)

# only keep central south asian
pqtl <- pqtl$CSA
head(pqtl)

# Get the p-value
pqtl$pval <- 10^-pqtl$LOG10P
head(pqtl)

# How many SNPs for one protein?
subset(pqtl, prot == "AKT3")

# How many snps for all proteins?
table(pqtl$prot)

# pqtl data doesn't have rsid - get that from combined_pqtl

pqtl_combined <- readRDS(file.path("data", "combined_pqtl.rds"))

head(pqtl_combined)

# add rsid to pqtl
pqtl <- left_join(
  pqtl,
  subset(pqtl_combined, select=c(CHROM, GENPOS, rsid)),
  by=c("CHROM", "GENPOS")
)

head(pqtl)


# Get the outcome data

load(file.path("data", "mdd_extract.rdata"))

# Look at CSA MDD GWAS
head(eur)
dim(eur)

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
head(csa)
out_dat <- format_data(csa,
                       type="outcome",
                       snp_col="rsid",
                       effect_allele_col="ALT",
                       other_allele_col="REF",
                       eaf_col="AF",
                       beta_col="ES",
                       se_col="SE",
                       pval_col="LP"
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

save(res, file=file.path("results", "csa.rdata"))
ggsave(file=file.path("results", "csaplot.png"))
