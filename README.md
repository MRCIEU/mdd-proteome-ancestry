# Hackathon group 6 - MR course, Uganda, 2023

## Objectives

Estimate the influence of pharmacological targets for incidence of major depressive disorder (MDD) across multiple ancestries using Mendelian randomisation

## Analysis

Download data - MDD GWASs plus the 

```
./scripts/download_mdd.sh
```

Identify proteins relevant to MDD in Europeans. Slight challenge - pQTL data is chr:pos, but MDD data is rsid. Need to map rsid to chr:pos using dbSNP, and then lookup pQTLs in MDD GWAS.

Need to have dbSNP VCF downloaded, and [bcftools](https://samtools.github.io/bcftools/bcftools.html) on the PATH.

1. Map MDD data to chr:position
2. Lookup cis-pQTLs in MDD data
3. Identify pQTL effects in MDD with FDR < 0.05
4. Write out list of proteins that have a cis effect in MDD

```
Rscript scripts/identify_pqtls.r
```

This identifies 17 proteins to use in the multi-ancestry comparison for effects on MDD.

The Eur pQTLs are ascertained for high LD tagging in Europeans only, so it will disadvantage power in non-Europeans. Use the cross-ancestry pQTL meta analysis to analyse the region around a pQTL. Choose the best SNP in the region for each of the 17 pQTLs.

```
Rscript scripts/cross_ancestry_pqtl_selection.r
```

Extract the ancestry-agnostic pQTL SNPs from each of the MDD GWASs

```
Rscript scripts/extract_mdd.r
```


TODO:

- PGC multi-ancestry MDD data - waiting for publication for summary stats to become available
- Linda is running GWAS on MDD in UGR
- Catherine running GWAS on MDD in ED6 and Neurogap (Uganda cohorts)
- Sandra to be included in future conversations / work
- Examine drugability of proteins being analysed - Anthony

