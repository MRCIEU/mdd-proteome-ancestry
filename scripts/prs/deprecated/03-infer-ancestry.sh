#!/bin/bash

# get paths etc
source config.sh

# Aim
# 1. get a list of markers in linkage equilibtrium
# 2. extract those markers from bgen to a single plink format dataset
# 3. use KING to project samples to ancestries


# 1. get a list of markers in linkage equilibtrium
cmd="
wget http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz
tar xzvf 1kg.v3.tgz
plink --bfile ld_files/EUR --indep 2000 500 1.05 --out 1kgeur
rm -r ld_files
rm 1kg.v3.tgz
"

dx run swiss-army-knife \
    -iin="/data/ancestry/KGref.fam.xz" \
    -icmd="${cmd}" \
    --destination="${project}:${datadir}/" \
    --brief \
    --yes

# 3. use KING to project samples to ancestsries


# download king and its binaries
# is it possible 
wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar xzvf Linux-king.tar.gz
dx upload king --destination /data/ancestry/

dx run swiss-army-knife \
    -iin="/data/ancestry/KGref.bed.xz" \
    -icmd= \
    --destination="${project}:${datadir}/" \
    --brief \
    --yes

dx run swiss-army-knife \
    -iin="/data/ancestry/KGref.bim.xz" \
    -icmd="wget https://www.kingrelatedness.com/ancestry/KGref.bim.xz; unxz KGref.bim.xz" \
    --destination="${project}:${datadir}/" \
    --brief \
    --yes

dx run swiss-army-knife \
    -iin="/data/ancestry/KGref.fam.xz" \
    -icmd="wget https://www.kingrelatedness.com/ancestry/KGref.fam.xz; unxz KGref.fam.xz" \
    --destination="${project}:${datadir}/" \
    --brief \
    --yes



# Download
wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar xzvf Linux-king.tar.gz
chmod 755 king
rm Linux-king.tar.gz

wget https://www.kingrelatedness.com/ancestry/KGref.bed.xz
unxz KGref.bed.xz

wget https://www.kingrelatedness.com/ancestry/KGref.bim.xz
unxz KGref.bim.xz

wget https://www.kingrelatedness.com/ancestry/KGref.fam.xz
unxz KGref.fam.xz

# Install

./king -b KGref.bed,ukb_pruned.bed --pca --projection --pngplot --prefix ukb

rm KGref.*
rm king
rm ukb_pruned*
"

dx run swiss-army-knife \
    -iin="/data/ancestry/KGref.fam.xz" \
    -icmd="wget https://www.kingrelatedness.com/ancestry/KGref.fam.xz; unxz KGref.fam.xz" \
    --destination="${project}:${datadir}/" \
    --brief \
    --yes


