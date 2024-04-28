# Hackathon group 6 - MR course, Uganda, 2023

## Objectives

Estimate the influence of pharmacological targets for incidence of major depressive disorder (MDD) across multiple ancestries using Mendelian randomisation

## Forward MR Analysis

1. Download data - MDD GWASs 

```
./scripts/download_mdd.sh
```

2. Identify proteins relevant to MDD in Europeans. Slight challenge - pQTL data is chr:pos, but MDD data is rsid. Need to map rsid to chr:pos using dbSNP, and then lookup pQTLs in MDD GWAS.

Need to have dbSNP VCF downloaded, and [bcftools](https://samtools.github.io/bcftools/bcftools.html) on the PATH.

    2.1. Map MDD data to chr:position
    2.2. Lookup cis-pQTLs in MDD data
    2.3. Identify pQTL effects in MDD with FDR < 0.05
    2.4. Write out list of proteins that have a cis effect in MDD

```
Rscript scripts/identify_pqtls.r
```

This identifies 17 proteins to use in the multi-ancestry comparison for effects on MDD.

3. The Eur pQTLs are ascertained for high LD tagging in Europeans only, so it will disadvantage power in non-Europeans. Use the cross-ancestry pQTL meta analysis to analyse the region around a pQTL. Choose the best SNP in the region for each of the 17 pQTLs.

```
Rscript scripts/cross_ancestry_pqtl_selection.r
```

4. Extract the ancestry-agnostic pQTL SNPs from each of the MDD GWASs

```
Rscript scripts/extract_mdd.r
```

5. Perform  per ancestry MR analysis
    EAS: `scripts/EAS_analysis_modified_Forward_MR.r`
	AFR: `scripts/afr_analysis_modified_Forward_MR.R`
	EUR: `scripts/eur_analysis_modified_Forward_MR.r`
	SAS: `scripts/csa_sas_modified_Foward_MR.r`


6. Perform combined MR and Heterogeneity analysis
    Use `scripts/combined_analysis.qmd`

## Reverse MR

1. Generated bim file with `scripts/bim_for_prs` script.
2. Calculate per ancestry PRS using PRScsx.py (refer to PRScsx), MDD per ancestry GWAS summary stats and bim file using `scripts/PRScsx.sh` .
3. Use `scripts/PRS/Association_Analysis.Rmd` to perform MDD PRS-protein association and heterogenity analysis.

## Foward MR for know druggable targets

1. Use script `scripts/identify_druggable_eqtls_Eur_ancestry` to obtain IVs for EUR and `scripts/identify_druggable_eqtls_All_ancestry` script for all ancestry GWAS data.
AFR, SAS and EAS have no significant Ivs

2. Perform per ancestry MR using rmd files 
	EAS:`scripts/Foward_mr_druggable_eas.Rmd`
	AFR: `scripts/Foward_mr_druggable_afr.Rmd`
	EUR: `scripts/Foward_mr_druggable_eur.Rmd`
	SAS: `scripts/Foward_mr_druggable_sas.Rmd`

3. Perform combined MR and heterogenity analysis using `scripts/combined_analysis _druggable_MR`
