# Analysis of MDD PRS in UK Biobank proteomics

We need to group individuals in UKB by ancestry, and then generate PRS scores for each individual based on ancestry. With these we can run associations against proteomic data. 

UKB scripting guidance taken from: https://github.com/pjgreer/ukb-rap-tools/tree/main/prs-calc


1. Identify ancestries in UKB individuals
2. Download proteomic data
3. Generate PRS for each ancestry
4. Perform analysis



1. Upload PRS scores
2. Get proteomic data and the list of IDs to include
3. Extract PRS SNPs and IDs from bgen files
4. Combine bgen files, filter for quality and combine into a single plink file
5. Generate scores



plink --bfile ~/repo/opengwas-api-internal/opengwas-api/app/ld_files/EUR --indep 2000 500 1.05 --out 1kgeur

