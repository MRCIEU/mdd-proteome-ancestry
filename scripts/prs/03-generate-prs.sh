#!/bin/bash

# Upload PRS rsults
dx upload -r ../../data/PRS/ --destination /mdd-prs/

# Generate rsid lists

cmd="
cut -f 2 combined_PRS_AFR.txt > combined_PRS_AFR.rsids
cut -f 2 combined_PRS_EUR.txt > combined_PRS_EUR.rsids
cut -f 2 combined_PRS_EAS.txt > combined_PRS_EAS.rsids
cut -f 2 combined_PRS_SAS.txt > combined_PRS_SAS.rsids
"

dx run swiss-army-knife \
    -iin="/mdd-prs/PRS/combined_PRS_AFR.txt" \
    -iin="/mdd-prs/PRS/combined_PRS_EUR.txt" \
    -iin="/mdd-prs/PRS/combined_PRS_EAS.txt" \
    -iin="/mdd-prs/PRS/combined_PRS_SAS.txt" \
    -icmd="${cmd}" \
    --tag="rsids" \
    --destination="${project}:/mdd-prs/PRS/" \
    --brief \
    --yes \
    --priority "high"

# Extract rsids and individuals from each chromosome

for i in {1..22}; do
    for anc in EUR EAS SAS AFR
    do
        cmd="
        # Extract from bgen
        qctool -g ${data_field}_c${i}_b0_v3.bgen -s ${data_field}_c${i}_b0_v3.sample -incl-rsids combined_PRS_${anc}.rsids -incl-samples ids_${anc}.txt -og chr_${anc}_${i}.gen -os chr_${anc}_${i}.sample

        # convert to plink
        plink --gen chr_${anc}_${i}.gen --sample chr_${anc}_${i}.sample --make-bed --out chr_${anc}_${i}
        rm chr_${anc}_${i}.gen chr_${anc}_${i}.sample

        # remove duplicates
        cut -f 2 chr_${anc}_${i}.bim | sort | uniq -c | awk '{ if(\$1 > 1) print \$2 }' > dupsnps.txt
        plink --bfile chr_${anc}_${i} --exclude dupsnps.txt --make-bed --out chr_${anc}_${i}

        rm dupsnps.txt chr*~
        "
        dx run swiss-army-knife \
            -iin="${imp_file_dir}/${data_field}_c${i}_b0_v3.bgen" \
            -iin="${imp_file_dir}/${data_field}_c${i}_b0_v3.sample" \
            -iin="${imp_file_dir}/${data_field}_c${i}_b0_v3.bgen.bgi" \
            -iin="/mdd-prs/PRS/combined_PRS_${anc}.rsids" \
            -iin="/data/ancestry/ids_${anc}.txt" \
            -icmd="${cmd}" \
            --tag="SelectSNPs${i}" \
            --instance-type="mem2_ssd2_v2_x16" \
            --destination="${project}:/mdd-prs/scratch/" \
            --brief \
            --yes \
            --priority "normal"
    done
done

for anc in EUR EAS SAS AFR
do
    mergecmd="
    rm -f mergelist.txt
    for i in {2..22}; do echo \"chr_${anc}_\${i}\" >> mergelist.txt; done

    plink --bfile chr_${anc}_1 --merge-list mergelist.txt --make-bed --out ukb_${anc}
    "

    dx run swiss-army-knife \
        -iin="${datadir}/chr_${anc}_1.bed" \
        -iin="${datadir}/chr_${anc}_1.bim" \
        -iin="${datadir}/chr_${anc}_1.fam" \
        -iin="${datadir}/chr_${anc}_2.bed" \
        -iin="${datadir}/chr_${anc}_2.bim" \
        -iin="${datadir}/chr_${anc}_2.fam" \
        -iin="${datadir}/chr_${anc}_3.bed" \
        -iin="${datadir}/chr_${anc}_3.bim" \
        -iin="${datadir}/chr_${anc}_3.fam" \
        -iin="${datadir}/chr_${anc}_4.bed" \
        -iin="${datadir}/chr_${anc}_4.bim" \
        -iin="${datadir}/chr_${anc}_4.fam" \
        -iin="${datadir}/chr_${anc}_5.bed" \
        -iin="${datadir}/chr_${anc}_5.bim" \
        -iin="${datadir}/chr_${anc}_5.fam" \
        -iin="${datadir}/chr_${anc}_6.bed" \
        -iin="${datadir}/chr_${anc}_6.bim" \
        -iin="${datadir}/chr_${anc}_6.fam" \
        -iin="${datadir}/chr_${anc}_7.bed" \
        -iin="${datadir}/chr_${anc}_7.bim" \
        -iin="${datadir}/chr_${anc}_7.fam" \
        -iin="${datadir}/chr_${anc}_8.bed" \
        -iin="${datadir}/chr_${anc}_8.bim" \
        -iin="${datadir}/chr_${anc}_8.fam" \
        -iin="${datadir}/chr_${anc}_9.bed" \
        -iin="${datadir}/chr_${anc}_9.bim" \
        -iin="${datadir}/chr_${anc}_9.fam" \
        -iin="${datadir}/chr_${anc}_10.bed" \
        -iin="${datadir}/chr_${anc}_10.bim" \
        -iin="${datadir}/chr_${anc}_10.fam" \
        -iin="${datadir}/chr_${anc}_11.bed" \
        -iin="${datadir}/chr_${anc}_11.bim" \
        -iin="${datadir}/chr_${anc}_11.fam" \
        -iin="${datadir}/chr_${anc}_12.bed" \
        -iin="${datadir}/chr_${anc}_12.bim" \
        -iin="${datadir}/chr_${anc}_12.fam" \
        -iin="${datadir}/chr_${anc}_13.bed" \
        -iin="${datadir}/chr_${anc}_13.bim" \
        -iin="${datadir}/chr_${anc}_13.fam" \
        -iin="${datadir}/chr_${anc}_14.bed" \
        -iin="${datadir}/chr_${anc}_14.bim" \
        -iin="${datadir}/chr_${anc}_14.fam" \
        -iin="${datadir}/chr_${anc}_15.bed" \
        -iin="${datadir}/chr_${anc}_15.bim" \
        -iin="${datadir}/chr_${anc}_15.fam" \
        -iin="${datadir}/chr_${anc}_16.bed" \
        -iin="${datadir}/chr_${anc}_16.bim" \
        -iin="${datadir}/chr_${anc}_16.fam" \
        -iin="${datadir}/chr_${anc}_17.bed" \
        -iin="${datadir}/chr_${anc}_17.bim" \
        -iin="${datadir}/chr_${anc}_17.fam" \
        -iin="${datadir}/chr_${anc}_18.bed" \
        -iin="${datadir}/chr_${anc}_18.bim" \
        -iin="${datadir}/chr_${anc}_18.fam" \
        -iin="${datadir}/chr_${anc}_19.bed" \
        -iin="${datadir}/chr_${anc}_19.bim" \
        -iin="${datadir}/chr_${anc}_19.fam" \
        -iin="${datadir}/chr_${anc}_20.bed" \
        -iin="${datadir}/chr_${anc}_20.bim" \
        -iin="${datadir}/chr_${anc}_20.fam" \
        -iin="${datadir}/chr_${anc}_21.bed" \
        -iin="${datadir}/chr_${anc}_21.bim" \
        -iin="${datadir}/chr_${anc}_21.fam" \
        -iin="${datadir}/chr_${anc}_22.bed" \
        -iin="${datadir}/chr_${anc}_22.bim" \
        -iin="${datadir}/chr_${anc}_22.fam" \
        -icmd="${mergecmd}" \
        --tag="mergeSNPs" \
        --instance-type "mem2_ssd2_v2_x16" \
        --destination="${project}:/mdd-prs/PRS/" \
        --brief \
        --yes \
        --priority "normal"
done

# remove intermediate files
dx rm -r /mdd-prs/scratch

# generate PRS
for anc in EUR AFR SAS EAS
do
    cmd="plink --bfile ukb_${anc} --score combined_PRS_${anc}.txt 2 4 6 --out mdd_score_${anc}"
    dx run swiss-army-knife \
        -iin="/mdd-prs/PRS/ukb_${anc}.bed" \
        -iin="/mdd-prs/PRS/ukb_${anc}.bim" \
        -iin="/mdd-prs/PRS/ukb_${anc}.fam" \
        -iin="/mdd-prs/PRS/combined_PRS_${anc}.txt" \
        -icmd="${cmd}" \
        --tag="scores_${anc}" \
        --instance-type="mem2_ssd2_v2_x16" \
        --destination="${project}:/mdd-prs/PRS/" \
        --brief \
        --yes \
        --priority "normal"
done
