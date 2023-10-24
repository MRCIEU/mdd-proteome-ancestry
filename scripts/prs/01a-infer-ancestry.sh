#!/bin/bash

# get paths etc
source config.sh

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

