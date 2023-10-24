#!/bin/bash

# get paths etc
source config.sh

# 2. extract those markers from bgen to a single plink format dataset

# Extraction of snps from chromosome bgen -> vcf -> plink
for i in {1..22}; do
    run_snps="
    # extract to vcf
    bgenix -g ${data_field}_c${i}_b0_v3.bgen \
    -incl-rsids ${rsidlist} \
    -vcf | \
    gzip -c > chr_${i}.vcf.gz
    
    # convert to plink
    plink --vcf chr_${i}.vcf.gz --make-bed --out chr_${i} --double-id

    # remove vcf
    rm chr_${i}.vcf.gz
    "
    dx run swiss-army-knife \
    -iin="${imp_file_dir}/${data_field}_c${i}_b0_v3.bgen" \
    -iin="${imp_file_dir}/${data_field}_c${i}_b0_v3.sample" \
    -iin="${imp_file_dir}/${data_field}_c${i}_b0_v3.bgen.bgi" \
    -iin="${datadir}/${rsidlist}" \
    -icmd="${run_snps}" \
    --tag="SelectSNPs${i}" \
    --instance-type "mem2_ssd2_v2_x16" \
    --destination="${project}:${datadir}" \
    --brief \
    --yes \
    --priority "normal"
done


mergecmd="
rm -f mergelist.txt
for i in {2..22}; do echo \"chr_\${i}\" >> mergelist.txt; done

plink --bfile chr_1 --merge-list mergelist.txt --make-bed --out ukb_pruned

# update variant ids to be just rsid

# first remove SNPs that are duplicate rsids
cut -f 2 ukb_pruned.bim | cut -d ';' -f 1 | sort | uniq -c | awk '{ if(\$1 > 1) print \$2\";\" }' > duprsids.txt
grep -f duprsids.txt ukb_pruned.bim | cut -f 2 > remsnps.txt
plink --bfile ukb_pruned --exclude remsnps.txt --make-bed --out ukb_pruned

# second update bim file to just have rsids
cut -f 2 ukb_pruned.bim | cut -d ';' -f 1 > ukb_pruned
paste ukb_pruned.bim ukb_pruned | \
awk '{ print \$1, \$7, \$3, \$4, \$5, \$6 }' > temp_\${i}
cp ukb_pruned.bim ukb_pruned.bim.orig
mv temp_\${i} ukb_pruned.bim
rm ukb_pruned duprsids.txt remsnps.txt

# Update IDs to be UKB ids
cat ${data_field}_c1_b0_v3.sample | sed 1d | sed 1d | awk '{ print \$1, \$2, \"0\", \"0\", \$4, \"-9\"}' > temp
cp ukb_pruned.fam ukb_pruned.fam.orig
mv temp ukb_pruned.fam
"

dx run swiss-army-knife \
    -iin="${imp_file_dir}/${data_field}_c1_b0_v3.sample" \
    -iin="${datadir}/chr_1.bed" \
    -iin="${datadir}/chr_1.bim" \
    -iin="${datadir}/chr_1.fam" \
    -iin="${datadir}/chr_2.bed" \
    -iin="${datadir}/chr_2.bim" \
    -iin="${datadir}/chr_2.fam" \
    -iin="${datadir}/chr_3.bed" \
    -iin="${datadir}/chr_3.bim" \
    -iin="${datadir}/chr_3.fam" \
    -iin="${datadir}/chr_4.bed" \
    -iin="${datadir}/chr_4.bim" \
    -iin="${datadir}/chr_4.fam" \
    -iin="${datadir}/chr_5.bed" \
    -iin="${datadir}/chr_5.bim" \
    -iin="${datadir}/chr_5.fam" \
    -iin="${datadir}/chr_6.bed" \
    -iin="${datadir}/chr_6.bim" \
    -iin="${datadir}/chr_6.fam" \
    -iin="${datadir}/chr_7.bed" \
    -iin="${datadir}/chr_7.bim" \
    -iin="${datadir}/chr_7.fam" \
    -iin="${datadir}/chr_8.bed" \
    -iin="${datadir}/chr_8.bim" \
    -iin="${datadir}/chr_8.fam" \
    -iin="${datadir}/chr_9.bed" \
    -iin="${datadir}/chr_9.bim" \
    -iin="${datadir}/chr_9.fam" \
    -iin="${datadir}/chr_10.bed" \
    -iin="${datadir}/chr_10.bim" \
    -iin="${datadir}/chr_10.fam" \
    -iin="${datadir}/chr_11.bed" \
    -iin="${datadir}/chr_11.bim" \
    -iin="${datadir}/chr_11.fam" \
    -iin="${datadir}/chr_12.bed" \
    -iin="${datadir}/chr_12.bim" \
    -iin="${datadir}/chr_12.fam" \
    -iin="${datadir}/chr_13.bed" \
    -iin="${datadir}/chr_13.bim" \
    -iin="${datadir}/chr_13.fam" \
    -iin="${datadir}/chr_14.bed" \
    -iin="${datadir}/chr_14.bim" \
    -iin="${datadir}/chr_14.fam" \
    -iin="${datadir}/chr_15.bed" \
    -iin="${datadir}/chr_15.bim" \
    -iin="${datadir}/chr_15.fam" \
    -iin="${datadir}/chr_16.bed" \
    -iin="${datadir}/chr_16.bim" \
    -iin="${datadir}/chr_16.fam" \
    -iin="${datadir}/chr_17.bed" \
    -iin="${datadir}/chr_17.bim" \
    -iin="${datadir}/chr_17.fam" \
    -iin="${datadir}/chr_18.bed" \
    -iin="${datadir}/chr_18.bim" \
    -iin="${datadir}/chr_18.fam" \
    -iin="${datadir}/chr_19.bed" \
    -iin="${datadir}/chr_19.bim" \
    -iin="${datadir}/chr_19.fam" \
    -iin="${datadir}/chr_20.bed" \
    -iin="${datadir}/chr_20.bim" \
    -iin="${datadir}/chr_20.fam" \
    -iin="${datadir}/chr_21.bed" \
    -iin="${datadir}/chr_21.bim" \
    -iin="${datadir}/chr_21.fam" \
    -iin="${datadir}/chr_22.bed" \
    -iin="${datadir}/chr_22.bim" \
    -iin="${datadir}/chr_22.fam" \
    -icmd="${mergecmd}" \
    --tag="mergeSNPs" \
    --instance-type "mem2_ssd2_v2_x16" \
    --destination="${project}:${datadir}" \
    --brief \
    --yes \
    --priority "normal"

for i in {1..22}
do
    dx rm /data/ancestry/chr_${i}.bed
    dx rm /data/ancestry/chr_${i}.bim
    dx rm /data/ancestry/chr_${i}.fam
    dx rm /data/ancestry/chr_${i}.log
    dx rm /data/ancestry/chr_${i}.nosex
done

