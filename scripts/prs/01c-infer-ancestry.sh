#!/bin/bash

# get paths etc
source config.sh

# 3. use KING to project samples to ancestries

cmd="
# Download king
wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar xzvf Linux-king.tar.gz
chmod 755 king
rm Linux-king.tar.gz

# Download reference datasets
wget https://www.kingrelatedness.com/ancestry/KGref.bed.xz
unxz KGref.bed.xz

wget https://www.kingrelatedness.com/ancestry/KGref.bim.xz
unxz KGref.bim.xz

wget https://www.kingrelatedness.com/ancestry/KGref.fam.xz
unxz KGref.fam.xz

# Install e1071 and data.table libraries
Rscript -e \"install.packages(c('e1071'))\"

# Run KING
./king -b KGref.bed,ukb_pruned.bed --pca --projection --pngplot --prefix ukb

# Remove unnecessary files
rm KGref.*
rm king
rm ukb_pruned*

# Get ancestry lists
awk '{ if (\$9 == \"EUR\") print \$1 }' ukb_InferredAncestry.txt > ids_EUR.txt
awk '{ if (\$9 == \"AFR\") print \$1 }' ukb_InferredAncestry.txt > ids_AFR.txt
awk '{ if (\$9 == \"SAS\") print \$1 }' ukb_InferredAncestry.txt > ids_SAS.txt
awk '{ if (\$9 == \"EAS\") print \$1 }' ukb_InferredAncestry.txt > ids_EAS.txt
awk '{ if (\$9 == \"AMR\") print \$1 }' ukb_InferredAncestry.txt > ids_AMR.txt
"

dx run swiss-army-knife \
    -iin="/data/ancestry/ukb_pruned.bed" \
    -iin="/data/ancestry/ukb_pruned.bim" \
    -iin="/data/ancestry/ukb_pruned.fam" \
    -icmd="${cmd}" \
    --destination="${project}:${datadir}/" \
    --brief \
    --yes
