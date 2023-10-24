#!/bin/bash

# Unfinished!
#
# Figure out how to extract data to csv
# dx extract_dataset downloads data to local
# table-exporter extracts to remote
# can just do this interactively but trying to determine how it's done programmatically here

sdir="../../data/prs"

mkdir -p ${sdir}

dx extract_dataset project-GX2v3q0JYBkyK8P8X8VgFgBv:app15825_20231011044137.dataset --list-entities

dx extract_dataset project-GX2v3q0JYBkyK8P8X8VgFgBv:app15825_20231011044137.dataset --entities olink_instance_0 --list-fields > ${sdir}/olink0_fields.txt
dx extract_dataset project-GX2v3q0JYBkyK8P8X8VgFgBv:app15825_20231011044137.dataset --entities olink_instance_2 --list-fields > ${sdir}/olink2_fields.txt
dx extract_dataset project-GX2v3q0JYBkyK8P8X8VgFgBv:app15825_20231011044137.dataset --entities olink_instance_3 --list-fields > ${sdir}/olink3_fields.txt

awk '{print $1}' ${sdir}/olink0_fields.txt > ${sdir}/olink0_fields_codes.txt
dx upload -r ${sdir}/olink0_fields_codes.txt --destination mdd-prs/

fieldlist=$(tr '\n' ',' < ${sdir}/olink0_fields_codes.txt | sed s/,$//1)
head -n 1 ${sdir}/olink0_fields_codes.txt

dx extract_dataset project-GX2v3q0JYBkyK8P8X8VgFgBv:app15825_20231011044137.dataset --fields eid,$(head -n 1 ${sdir}/olink0_fields_codes.txt) --output ${sdir}/olink0.csv

wc -l ${sdir}/olink3_fields.txt
wc -l ${sdir}/olink0_fields.txt






dx run table-exporter \
    -idataset_or_cohort_or_dashboard="project-GX2v3q0JYBkyK8P8X8VgFgBv:app15825_20231011044137.dataset" \
    -field_names_file_txt="project-GX2v3q0JYBkyK8P8X8VgFgBv:mdd-prs/olink0_fields_codes.txt" \
    -