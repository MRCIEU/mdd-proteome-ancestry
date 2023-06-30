# Datasets

- `protlist_liu.txt` - List of proteins for MDD identified in Liu et al (these are ignored as they were not available for analysis in UKB)
- `protlist_ukb.txt` - List of 17 proteins identified in EUR UKB participants that associate with EUR MDD - used for the rest of the analysis
- `protlist_ukb_chrpos.txt` - Info about the pQTLs for 17 proteins identified for EUR
- `combined_pqtl.rds` - R data file that has a single data frame containing the top 20 SNPs for each of the 17 protein levels, based on the combined analysis of all ancestries. This file also contains the chr:pos for builds 37 and 38, as well as the rsids for each of the SNPs. It also contains the overall effect estimates for the SNP-exposure associations (exposure = protein level)
- `ancestry_pqtl.rds` - R data file containing a list of data frames, one for each of the 4 ancestries that is being analysed. Each data frame is equivalent to the `combined_pqtl.rds`, but the effect estimates are based on using only one ancestry. So these are the ancestry-specific SNP-exposure estimates. Note that the rsid is missing from this file, so you could merge it with `combined_pqtl.rds` to match rsid to chr:pos
- `mdd_extract.rdata` - R data file that contains 4 data frames, one for each ancestry. The data frames are the ancestry-specific SNP-outcome associations for each of the pQTLs extracted in `combined_pqtl.rds`.

