#!/bin/bash

mkdir -p data

# East Asian MDD
# Giannakopoulu et al 2021 - https://pubmed.ncbi.nlm.nih.gov/34586374/
wget -O ../scratch/mdd_eas.txt.gz https://figshare.com/ndownloader/files/31424374
# wget https://figshare.com/ndownloader/files/34437842

# CSA MDD
# UKB Pan GWAS
wget -O ../scratch/ukb-e-20126_p5_CSA.vcf.gz https://gwas.mrcieu.ac.uk/files/ukb-e-20126_p5_CSA/ukb-e-20126_p5_CSA.vcf.gz
wget -O ../scratch/ukb-e-20126_p5_CSA.vcf.gz.tbi https://gwas.mrcieu.ac.uk/files/ukb-e-20126_p5_CSA/ukb-e-20126_p5_CSA.vcf.gz.tbi

# AFR MDD
# UKB Pan GWAS
wget -O ../scratch/ukb-e-20126_p1_AFR.vcf.gz https://gwas.mrcieu.ac.uk/files/ukb-e-20126_p1_AFR/ukb-e-20126_p1_AFR.vcf.gz
wget -O ../scratch/ukb-e-20126_p1_AFR.vcf.gz.tbi https://gwas.mrcieu.ac.uk/files/ukb-e-20126_p1_AFR/ukb-e-20126_p1_AFR.vcf.gz.tbi

# EUR MDD
# Howard et al 2019 - https://datashare.ed.ac.uk/handle/10283/3203
wget -O ../scratch/mdd_eur.txt https://datashare.ed.ac.uk/bitstream/handle/10283/3203/PGC_UKB_depression_genome-wide.txt?sequence=3&isAllowed=y
gzip ../scratch/mdd_eur.txt

# UKB pQTL EUR
# Sun et al 2023 - https://europepmc.org/article/ppr/ppr508031#supplementary-material
wget -O ../scratch/media-2.xlsx https://www.biorxiv.org/content/biorxiv/early/2022/06/18/2022.06.17.496443/DC2/embed/media-2.xlsx?download=true

