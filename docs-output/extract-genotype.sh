#!/bin/bash
# INPUT="/Users/hung/Data/DGV4VN/Impute_protocol/run_1/array/test100.chr20_QC3.vcf.gz"
INPUT=$1
input=$(basename ${INPUT})
TMP=${input%.vcf.gz}_ID.vcf.gz
OUTPUT="${input%_*}_recodeAt"
# subset based on list of sample and annotate ID with bcftools
bcftools query -l ../Impute/APMRA19.hg38.chr20_impute.1KGP3.dose.vcf.gz > sample_test.txt
bcftools view -S sample_test.txt ${INPUT} -Ou |\
  bcftools annotate --set-id +'%CHROM:%POS:%REF:%ALT' -Oz -o ${TMP}
# recode with plink
plink --vcf ${TMP} \
    --const-fid --keep-allele-order --recode A-transpose \
    --out ${OUTPUT}
rm ${TMP}