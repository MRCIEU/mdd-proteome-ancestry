# Analysis of MDD PRS in UK Biobank proteomics

We need to generate PRS scores for each individual in UKB using the RAP. With these we can run associations against proteomic data.

Guidance taken from: https://github.com/pjgreer/ukb-rap-tools/tree/main/prs-calc



1. Upload PRS scores
2. Get proteomic data and the list of IDs to include
3. Extract PRS SNPs and IDs from bgen files
4. Combine bgen files, filter for quality and combine into a single plink file
5. Generate scores
