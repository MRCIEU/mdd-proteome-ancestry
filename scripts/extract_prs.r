## Notes
## Trying to make a performant way to extract data from ukb to perform the PRS assoc using summary level information
## However the MDD PRS is LD adjusted whereas the PRS is not.
## A couple of tests below suggest that we'll be getting overly precise estimates from the analysis here.
## Better to generate PRS in the individual level data and then perform the assoc in a standard way


library(queryukbppp)
library(here)
library(data.table)
library(dplyr)

get_snpid_from_chrpos_r <- function(chrpos, build=c("19", "38")[1], map_files=get_mapfiles()) {
    td <- tempdir(check=TRUE)
    chrpos <- strsplit(chrpos, ":")
    chrpos <- dplyr::tibble(chr=sapply(chrpos, \(x) x[1]), pos=sapply(chrpos, \(x) x[2]))
    l <- list()
    for(ch in unique(chrpos$chr)) {
        message(ch)
        # utils::write.table(subset(chrpos, chr == ch)$pos, file=file.path(td, "pos"), row.names=FALSE, col.names=FALSE, quote=FALSE)
        # cmd <- glue::glue("zgrep -wf {file.path(td, 'pos')} '{grep(paste0('chr', ch), map_files, value=TRUE)}' > {file.path(td, 'out')}")
        # system(cmd)
        # l[[ch]] <- data.table::fread(file.path(td, 'out'))
        fn <- grep(paste0('chr', ch, "_"), map_files, value=TRUE)
        message(fn)
        stopifnot(file.exists(fn))
        x <- fread(fn)
        x$chr <- ch
        l[[ch]] <- subset(x, POS19 %in% subset(chrpos, chr == ch)$pos)
        
    }
    a <- bind_rows(l)
    return(a)    
}

standardise <- function(d, ea="ea", oa="oa", eaf="eaf", beta="beta", chr="chr", pos="pos") {
    toflip <- d[[ea]] > d[[oa]]
    if(!is.null(d[[eaf]])) {
        d[[eaf]][toflip] <- 1 - d[[eaf]][toflip]
    }
    if(!is.null(d[[beta]])) {
        d[[beta]][toflip] <- d[[beta]][toflip] * -1
    }
    temp <- d[[oa]][toflip]
    d[[oa]][toflip] <- d[[ea]][toflip]
    d[[ea]][toflip] <- temp
    d[["varid"]] <- paste0(d[[chr]], ":", d[[pos]], "_", d[[ea]], "_", d[[oa]])
    d
}

query_pgwas_r <- function(tarfile, variants) {
    td <- tempdir()
    cmd <- glue::glue("tar xvf '{tarfile}' -C '{td}'")
    system(cmd)
    fd <- file.path(td, gsub(".tar", "", basename(tarfile)))
    fnut <- list.files(fd, full.names = TRUE)
    l <- list()
    for (ch in unique(variants$chr)) {
        message(ch)
        fn <- grep(paste0("chr", ch, "_"), fnut, value = T)
        if (!file.exists(fn)) {
            message("missing")
            next
        }
        x <- fread(fn)
        l[[ch]] <- subset(x, ID %in% variants$ID)
    }
    l <- bind_rows(l)
    con <- gzfile(fn)
    names(l) <- scan(con, nlines = 1, what = character())
    close(con)
    system(glue::glue("rm -r {fd}"))
    l <- inner_join(l, v, by="ID") %>%
        filter(!duplicated(ID)) %>%
        select(chr, pos=POS19, oa=ALLELE0, ea=ALLELE1, eaf=A1FREQ, beta=BETA, se=SE) %>%
        standardise()
    return(l)
}

prot <- protein_info()

a <- file.path(options()$ukbpppdir, "UKB-PPP pGWAS summary statistics", "European (discovery)") %>% list.files(full.names=TRUE)
length(a)

prs <- fread(here("data", "PRS", "combined_PRS_EUR.txt")) %>% 
    select(chr=V1, rsid=V2, pos=V3, ea=V4, oa=V5, beta=V6) %>%
    standardise()
prs

chrpos <- paste0(prs$V1, ":", prs$V3)
v <- get_snpid_from_chrpos_r(paste0(prs$V1, ":", prs$V3))

prs

t1 <- Sys.time()
o <- query_pgwas_r(a[1], v)
t1 - Sys.time()

head(prs)
head(o)


standardise

x <- inner_join(prs, o, by="varid")
summary(lm(beta.y ~ 0 + beta.x, data=x, weight=1/x$se^2))

nvar <- seq(10000, 900000, by=10000)
res1 <- lapply(nvar, \(n) {
    index <- sample(1:nrow(x), n)
    summary(lm(beta.y ~ 0 + beta.x, data=x[index], weight=1/x$se[index]^2))$coef %>% as_tibble() %>% slice_tail(n=1) %>% mutate(nvar=n)
}) %>% bind_rows()

res1
plot(Estimate ~ nvar, res1)
summary(lm(Estimate ~ nvar, res1))


o100 <- query_pgwas_r(a[100], v)
x <- inner_join(prs, o100, by="varid")
summary(lm(beta.y ~ 0 + beta.x, data=x, weight=1/x$se^2))
