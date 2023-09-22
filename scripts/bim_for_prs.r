library(data.table)
library(dplyr)
library(tidyr)

# Strategy
## Find common SNPs across 1000G super populations
## Find which of those are in the combined pQTL sum stats
## Write out bim

p <- "/projects/MRC-IEU/research/projects/icep2/wp1/028/working/data/ukbb_pqtl/gib/combined/BTN3A2_P78410_OID20540_v1_Inflammation"

f <- file.path(p, list.files(p))
file.exists(f)

readchr <- function(f) {
    message(f)
    fread(f) %>%
        tidyr::separate(ID, sep=":", into=c("chr", "pos", "a1", "a2", "imp", "v1")) %>%
        dplyr::select(chr, pos, a1, a2)

}

a <- lapply(f, readchr) %>% bind_rows()


library(glue)
pop <- 
bims <- lapply(c("EUR", "AFR", "SAS", "EAS"), \(x) {
        glue("~/scratch/data/ldmat/{x}.bim") %>%
        fread(.) %>% 
        mutate(pop=x)
    })

str(bims)

rsids <- lapply(bims, \(x) x$V2) %>%
    Reduce(intersect, .)

length(rsids)

bim <- lapply(bims, \(x) subset(x, V2 %in% rsids))
str(bim)

head(readchr)

a$chrpos <- paste(a$chr, a$pos)
head(bim)

b <- bim[[1]]
b$chrpos <- paste(b$V1, b$V4)
head(b)

b <- subset(b, chrpos %in% a$chrpos)
write.table(b %>% select(V1, V2, V3, V4, V5, V6), file="data/prs.bim", row=F, col=F, qu=F)
