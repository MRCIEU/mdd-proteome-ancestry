library(glue)
library(dplyr)
library(here)
library(readxl)
library(data.table)
library(tidyr)

bcftools <- "bcftools"
vcffile <- "/mnt/storage/private/mrcieu/research/mr-eve/vcf-reference-datasets/dbsnp/dbsnp.v153.b37.vcf.gz"

generate_vid <- function(d, ea = "ea", nea = "nea", eaf = "eaf", beta = "beta", rsid = "rsid", chr = "chr", position = "position") {
  toflip <- d[[ea]] > d[[nea]]
  d[[eaf]][toflip] <- 1 - d[[eaf]][toflip]
  d[[beta]][toflip] <- d[[beta]][toflip] * -1
  temp <- d[[nea]][toflip]
  d[[nea]][toflip] <- d[[ea]][toflip]
  d[[ea]][toflip] <- temp
  d[["rsido"]] <- d[[rsid]]
  d[[rsid]] <- paste0(d[[chr]], ":", d[[position]], "_", toupper(strtrim(d[[ea]], 5)), "_", toupper(strtrim(d[[nea]], 5)))
  d
}

a <- read_xlsx(here("scratch", "media-2.xlsx"), sheet="ST6", skip=2) %>%
    tidyr::separate(`Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`, sep=":", into=c("chr37","pos37","oa", "ea", "imp", "v")) %>%
    mutate(code=paste0(chr37, ":", pos37))

a <- generate_vid(a, "ea", "oa", "A1FREQ (discovery)", "BETA (discovery, wrt. A1)", "code", "chr37", "pos37")

pqtl_cis <- subset(a, `cis/trans` == "cis") %>%
    mutate(code=paste0(chr37, ":", pos37))


# b37
mdd <- fread(here("scratch", "PGC_UKB_depression_genome-wide.txt.gz"))

rsid <- mdd$MarkerName
tmp <- tempfile()
write.table(unique(rsid), file=paste0(tmp, ".snplist"), row.names = FALSE, col.names = FALSE, quote = FALSE)
cmd <- glue("{bcftools} query -f '%CHROM %POS %ID\\n' -i'ID=@{tmp}.snplist' {vcffile} > {tmp}.txt")
system(cmd)

map <- fread(glue("{tmp}.txt"))
names(map) <- c("chr", "pos", "rsid")
mdd <- left_join(mdd, map, by=c("MarkerName" = "rsid"))
str(mdd)
mdd <- generate_vid(mdd, "A1", "A2", "Freq", "LogOR", "MarkerName", "chr", "pos")

mdd_pqtl <- subset(mdd, code %in% pqtl_cis$code)
str(mdd_pqtl)

mdd_pqtl$fdr <- p.adjust(mdd_pqtl$P, "fdr")
table(mdd_pqtl$fdr < 0.05)

mdd_pqtl_sig <- subset(mdd_pqtl, fdr < 0.05)

pqtl_cis_sig <- subset(pqtl_cis, code %in% mdd_pqtl_sig$code)

write.table(pqtl_cis_sig$`Assay Target`, file=here("data", "protlist_ukb.txt"), row=F, col=F, qu=F)

prot <- subset(pqtl_cis, code %in% mdd_pqtl_sig$code) %>%
  dplyr::select(chr37, pos37, rsID, prot=`Assay Target`)

write.table(prot, file=here("data", "protlist_ukb_chrpos.txt"), row=F, col=F, qu=F)

