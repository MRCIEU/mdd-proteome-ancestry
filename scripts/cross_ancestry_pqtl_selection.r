library(glue)
library(here)
library(dplyr)
library(data.table)
library(parallel)
options(mc.cores=5)
library(tidyr)

# Get the EUR pQTLs that have an assoc with EUR MDD
prot <- read.table(here("data", "protlist_ukb_chrpos.txt"))
names(prot) <- c("chr", "pos", "rsid", "prot")

# This is where the combined pQTL data are
combined_path <- "/projects/MRC-IEU/research/projects/icep2/wp1/028/working/data/ukbb_pqtl/gib/combined"
a <- list.files(combined_path)
a <- a[!grepl(".tar", a)]
a

# Create a data frame that has the path to each of the chromosome files that need to be used to identify the best SNP
prot$fullname <- a
prot
prot$filename <- paste0("combined_chr", prot$chr, "_", gsub("_", ":", prot$fullname), ".gz")
prot$path <- file.path(combined_path, prot$fullname, prot$filename)

# Check all the necessary files are there
stopifnot(all(file.exists(file.path(combined_path, prot$fullname, prot$filename))))


# Function to identify best SNP in region
# Just in case some SNPs are missing from some datasets, take the top 10 SNPs for each region
get_best_snp <- function(prot, i, radius) {
    x <- fread(prot$path[i]) %>%
        dplyr::filter(LOG10P > -log10(0.01)) %>%
        tidyr::separate(ID, sep=":", into=c("chr", "pos", "a0", "a1", "imp", "v1")) %>%
        mutate(pos = as.numeric(pos)) %>%
        subset(pos > (prot$pos[i] - radius) & pos < (prot$pos[i] + radius)) %>%
        arrange(desc(LOG10P)) %>%
        slice_head(n=20) %>%
        mutate(prot = prot$prot[i])
    x
}

# check
get_best_snp(prot, 1, 250000)

# Loop through all proteins and get best SNP
prot2 <- mclapply(1:nrow(prot), function(i) {
    get_best_snp(prot, i, 250000)
}, mc.cores=5) %>% bind_rows()

# Extract cross-ancestry pQTLs from each specific ancestry


paths <- tibble(
    pop=c("EUR", "EAS", "CSA", "AFR"),
    fp = c(
        "/projects/MRC-IEU/research/projects/icep2/wp1/028/working/data/ukbb_pqtl/",
        "/projects/MRC-IEU/research/projects/icep2/wp1/028/working/data/ukbb_pqtl/gib/east_asian",
        "/projects/MRC-IEU/research/projects/icep2/wp1/028/working/data/ukbb_pqtl/gib/central_south_asian",
        "/projects/MRC-IEU/research/projects/icep2/wp1/028/working/data/ukbb_pqtl/gib/african"
    )
)

o <- lapply(1:nrow(paths), function(i) {
    message(paths$pop[i])
    path <- file.path(paths$fp[i], prot$fullname, prot$filename)
    if(paths$pop[i] == "EUR")
    {
        path <- gsub("combined", "discovery", path)
        # TNR not available in EUR? just used combined
        tnrloc <- grep("TNR", path)
        path[tnrloc] <- prot$path[tnrloc]        
    }
    if(paths$pop[i] == "EAS")
    message(paste(file.exists(path), collapse=" "))
    path <- path[file.exists(path)]
    l <- mclapply(1:length(path), function(j) {
        message(j)
        protj <- prot$prot[j]
        ids <- subset(prot2, prot == protj)$GENPOS
        x <- fread(path[j]) %>% mutate(prot=protj)
        subset(x, GENPOS %in% ids)
    }, mc.cores=5)
    o <- bind_rows(l)
    o$pop <- paths$pop[i]
    return(o)
})


lapply(o, function(x) length(unique(x$prot)))

names(o) <- paths$pop
saveRDS(o, file=here("data", "ancestry_pqtl.rds"))

# get rsids for prot2

bcftools <- "bcftools"
vcffile <- "/mnt/storage/private/mrcieu/research/mr-eve/vcf-reference-datasets/dbsnp/dbsnp.v153.b37.vcf.gz"

tmp <- tempfile()
write.table(tibble(prot2$chr, a=prot2$pos, b=prot2$pos), file=paste0(tmp, ".snplist"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")
cmd <- glue("{bcftools} query -f '%CHROM %POS %ID\\n' -R {tmp}.snplist {vcffile} > {tmp}.txt")
system(cmd)

map <- fread(glue("{tmp}.txt"))
names(map) <- c("chr", "pos", "rsid")
map <- subset(map, !duplicated(rsid))
map <- subset(map, !duplicated(paste(chr, pos)))
dim(map)
map$chr <- as.character(map$chr)
prot3 <- left_join(prot2, map, by=c("chr", "pos"))
dim(prot3)
prot3
table(is.na(prot3$rsid))

saveRDS(prot3, file=here("data", "combined_pqtl.rds"))


