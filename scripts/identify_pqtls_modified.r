##Loading Libraries
library(glue)
library(dplyr)
library(here)
library(readxl)
library(data.table)
library(tidyr)


#Setting File Paths:
##bcftools <- "/usr/bin/bcftools"
vcffile <- "scratch/dbsnp.v153.b37.vcf.gz"


#Defining the generate_vid function that:Checks if a variant should be flipped based on the effect allele (ea) and the non-effect allele (nea), Adjusts allele frequencies and effect size (beta) values based on the above flipping and Generates a new ID for the variant using the chromosome, position, and allele information.
generate_vid <- function(d, ea = "ea", nea = "nea", eaf = "eaf", beta = "beta", rsid = "rsid", chr = "chr", position = "position") {
  toflip <- d[[ea]] > d[[nea]]
  d[[eaf]][toflip] <- 1 - d[[eaf]][toflip]
  d[[beta]][toflip] <- d[[beta]][toflip] * -1
  temp <- d[[nea]][toflip]
  # d[[nea]][toflip] <- d[[ea]][toflip]
  d[[ea]][toflip] <- temp
  d[["rsido"]] <- d[[rsid]]
  d[[rsid]] <- paste0(d[[chr]], ":", d[[position]], "_", toupper(strtrim(d[[ea]], 5)), "_", toupper(strtrim(d[[nea]], 5)))
  d
}


#Reading and Transforming the ukpQTL data: The column Variant ID is split into several columns and A new code column is created from the chromosome and position.

a <- read_xlsx(here("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon/scratch", "media-2.xlsx"), sheet="ST6", skip=2) %>%
  tidyr::separate(`Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`, sep=":", into=c("chr37","pos37","oa", "ea", "imp", "v")) %>%
  mutate(code=paste0(chr37, ":", pos37))


a <- generate_vid(a, "ea", "oa", "A1FREQ (discovery)", "BETA (discovery, wrt. A1)", "code", "chr37", "pos37")

#Filtering pQTL data to only retain variants that are in cis

pqtl_cis <- subset(a, `cis/trans` == "cis") %>%
    mutate(code=paste0(chr37, ":", pos37))


#Reading in the MDD GWAS summary statistics for all ancestries and extracting rsids 
mdd <- fread(here("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon/scratch", "FE_all_ancestry_exclu_whi_jhs_23andMe_qced_rsid.txt.gz"))

rsid <- mdd$rsid

output_path <- "C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon/scratch" 
snplist_filename <- file.path(output_path, "mysnplist.snplist")
txt_filename <- file.path(output_path, "newoutput.txt")


write.table(unique(rsid), file = snplist_filename, row.names = FALSE, col.names = FALSE, quote = FALSE)

##querying the VCF file using bcftools 
#cmd <- "bcftools query -f '%CHROM %POS %ID\n' -i'ID=@{snplist_filename}' {vcffile} > {txt_filename}"
#system(cmd)


map <- fread(txt_filename)

names(map) <- c("chr", "pos", "rsid")
mdd <- left_join(mdd, map, by=c("rsid" = "rsid"))
str(mdd)
mdd_new <- generate_vid(mdd, "A1", "A2", "Freq", "LogOR", "MarkerName", "chr", "pos")
str(mdd_new)
##creating a code column in new-mdd call it mdd_code
mddcode<- mdd_new %>% mutate(code=paste0(Chromosome, ":", Position))
str(mddcode)


#Subsetting by rsid column
mdd_pqtlold <- subset(mddcode, rsid %in% pqtl_cis$rsID)
str(mdd_pqtlold)



##Subsetting by code column
mdd_pqtlcode <- subset(mddcode, code %in% pqtl_cis$code)
str(mdd_pqtlcode)





##Filtering and writing out put from Subsetting by code column
mdd_pqtlcode$fdr <- p.adjust(mdd_pqtlcode$P, "fdr")
table(mdd_pqtlcode$fdr < 0.05)

mdd_pqtlcode_sig <- subset(mdd_pqtlcode, fdr < 0.05)

pqtl_cis_sigcode <- subset(pqtl_cis, code %in% mdd_pqtlcode_sig$code)

write.table(pqtl_cis_sigcode$`Assay Target`, file=here("data", "Allprotlistcode.txt"), row=F, col=F, qu=F)

protcode <- subset(pqtl_cis, code %in% mdd_pqtlcode_sig$code) %>%
  dplyr::select(chr37, pos37, rsID, prot=`Assay Target`)

write.table(protcode, file=here("data", "Allprotlist_chrposcode.txt"), row=F, col=F, qu=F)



##Filtering and writing out put from Subsetting by rsid column
mdd_pqtlold$fdr <- p.adjust(mdd_pqtlold$P, "fdr")
str(mdd_pqtlold)
table(mdd_pqtlold$fdr < 0.05)

mdd_pqtl_sig <- subset(mdd_pqtlold, fdr < 0.05)

pqtl_cis_sig <- subset(pqtl_cis, rsID %in% mdd_pqtl_sig$rsid)

write.table(pqtl_cis_sig$`Assay Target`, file=here("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon/data", "Allrsidprotlist.txt"), row=F, col=F, qu=F)

prot <- subset(pqtl_cis, rsID %in% mdd_pqtl_sig$rsid) %>%
  dplyr::select(chr37, pos37, rsID, prot=`Assay Target`)

write.table(prot, file=here("C:/Users/HP OMEN GAMING/Desktop/mr-uganda-hackathon/data", "Allrsidprotlist_chrpos.txt"), row=F, col=F, qu=F)


