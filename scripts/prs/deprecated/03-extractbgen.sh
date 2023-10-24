#!/bin/bash

project="project-GX2v3q0JYBkyK8P8X8VgFgBv"
imp_file_dir="/Bulk/Imputation/UKB imputation from genotype/"
data_field="ukb22828"
datadir="/data/ancestry"
rsidlist="1kgeur.prune.in"

echo $rsidlist

# Get list of independent SNPs
wget http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz
tar xzvf 1kg.v3.tgz
plink --bfile ld_files/EUR --indep 2000 500 1.05 --out 1kgeur
rm -r ld_files
rm 1kg.v3.tgz

dx upload 1kgeur.prune.in --destination ${datadir}/

# King
wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar xzvf Linux-king.tar.gz
dx upload king --destination /data/ancestry/

wget https://www.kingrelatedness.com/ancestry/KGref.bed.xz
# unxz KGref.bed.xz
dx upload KGref.bed --destination /data/ancestry/
dx run swiss-army-knife \
    -iin="/data/ancestry/KGref.bed.xz" \
    -icmd="unxz KGref.bed.xz" \
    --destination="${project}:${datadir}/" \
    --brief \
    --yes

wget https://www.kingrelatedness.com/ancestry/KGref.bim.xz
dx upload KGref.bim --destination /data/ancestry/
dx run swiss-army-knife \
    -iin="/data/ancestry/KGref.bim.xz" \
    -icmd="unxz KGref.bim.xz" \
    --destination="${project}:${datadir}/" \
    --brief \
    --yes

wget https://www.kingrelatedness.com/ancestry/KGref.fam.xz
dx upload KGref.fam --destination /data/ancestry/
dx run swiss-army-knife \
    -iin="/data/ancestry/KGref.fam.xz" \
    -icmd="unxz KGref.fam.xz" \
    --destination="${project}:${datadir}/" \
    --brief \
    --yes

rm KGref*

# Extraction of snps from bgen -> vcf -> plink 
for i in {1..22}; do
    run_snps="
    bgenix -g ${data_field}_c${i}_b0_v3.bgen \
    -incl-rsids ${rsidlist} \
    -vcf | \
    gzip -c > chr_${i}.vcf.gz
    
    plink --vcf chr_${i}.vcf.gz --make-bed --out chr_${i} --double-id

    cut -f 2 chr_${i}.bim | cut -d ';' -f 1 > chr_${i}
    
    paste chr_${i}.bim chr_${i} | \
    awk '{ print \$1, \$7, \$3, \$4, \$5, \$6 }' > temp_${i}
    cp chr_${i}.bim chr_${i}.bim.orig
    mv temp_${i} chr_${i}.bim

    rm chr_${i} chr_${i}.vcf.gz
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

combine_snps='
touch mergelist.txt

for i in {2..22}; do echo "/mnt/project/data/ancestry/chr_${i}" >> mergelist.txt; done

plink --bfile /mnt/project/data/ancestry/chr_1 --merge-list mergelist.txt --make-bed --out ukb_pruned
'

dx run swiss-army-knife \
    -icmd="${combine_snps}" \
    --tag="CombSNPs" \
    --instance-type "mem2_ssd2_v2_x16" \
    --destination="${project}:${data_file_dir}" \
    --brief \
    --yes

# this gives a bgen error


combine_snps='cat-bgen -g /Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c21_b0_v3.bgen /Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c22_b0_v3.bgen -og c21_22.bgen -clobber'

# try with the unprocessed bgen files
dx run swiss-army-knife \
    -iin="${imp_file_dir}/${data_field}_c1_b0_v3.sample" \
    -icmd="${combine_snps}" \
    --tag="CombSNPs" \
    --instance-type "mem2_ssd2_v2_x16" \
    --destination="${project}:${data_file_dir}" \
    --brief \
    --yes




# Calculate ancestry
cmd='${datadir}/king -b ${datadir}/KGref.bed,${datadir}/ukb_pruned.bed --pca --projection --pngplot --prefix ${datadir}/ukb'

dx run swiss-army-knife \
    -iin="$datadir/ukb_pruned" \
    -iin="$datadir/KGref.bed" \
    -iin="$datadir/KGref.bim" \
    -iin="$datadir/KGref.fam" \
    -iin="${imp_file_dir}/${data_field}_c1_b0_v3.sample" \
    -iin="${txt_file_dir}/${scorefile}" \
    -iin="${txt_file_dir}/${prsfile}" \
    -icmd="${combine_snps}" \
    --tag="CombSNPs" \
    --instance-type "mem2_ssd2_v2_x16" \
    --destination="${project}:${data_file_dir}" \
    --brief \
    --yes






wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar xzvf Linux-king.tar.gz
ls -l

wget https://www.kingrelatedness.com/ancestry/KGref.fam.xz




wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220402.zip && \
    mkdir plinktemp && mv plink_linux_x86_64_20220402.zip plinktemp && cd plinktemp && \
    unzip plink_linux_x86_64_20220402.zip && \
    mv plink /bin/ && \
    cd ../ && rm -r plinktemp
    
touch mergelist.txt

for i in {2..22}; do echo "/mnt/project/data/ancestry/chr_${i}.bed /mnt/project/data/ancestry/chr_${i}.bim.orig /mnt/project/data/ancestry/chr_${i}.fam" >> mergelist.txt; done

plink --bfile /mnt/project/data/ancestry/chr_1 --bim /mnt/project/data/ancestry/chr_1.bim.orig --merge-list mergelist.txt --make-bed --out ukb_pruned



# need to remove SNPs that are duplicate positions





docker pull lifebitai/bgenix

docker run -v /mntlifebitai/bgenix

combine_snps='
touch mergelist.txt

for i in {2..22}; do echo "/data/ancestry/chr_${i} >> mergelist.txt"; done

plink --bfile /data/ancestry/chr_1 --merge-list mergelist.txt --make-bed --out ukb_pruned
'


## remove duplicate snps

wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220402.zip && \
    mkdir plinktemp && mv plink_linux_x86_64_20220402.zip plinktemp && cd plinktemp && \
    unzip plink_linux_x86_64_20220402.zip && \
    mv plink /bin/ && \
    cd ../ && rm -r plinktemp

bfile="/mnt/project/data/ancestry/ukb_pruned"
cp $bfile.* .

cut -f 2 $bfile.bim | cut -d ";" -f 1 | sort | uniq -c | awk '{ if($1 > 1) print $2";" }' > duprsids.txt
grep -f duprsids.txt $bfile.bim | cut -f 2 > remsnps.txt
plink --bfile $bfile --exclude remsnps.txt --make-bed --out ukb_pruned

# update bim file to just have rsids
cut -f 2 ukb_pruned.bim | cut -d ';' -f 1 > ukb_pruned
paste ukb_pruned.bim ukb_pruned | \
awk '{ print $1, $7, $3, $4, $5, $6 }' > temp_${i}
cp ukb_pruned.bim ukb_pruned.bim.orig
mv temp_${i} ukb_pruned.bim
rm ukb_pruned


# king
cp /mnt/project/data/ancestry/king .
chmod 755 king
cp /mnt/project/data/ancestry/KGref* .

./king -b KGref.bed,ukb_pruned.bed --pca --projection --pngplot --prefix ukb
