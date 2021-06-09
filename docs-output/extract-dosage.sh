#!/bin/bash
# extract imputed snp dosage
INPUT=$1
_input=$(basename $INPUT)
_input2=${_input#*impute.}
_input3=${_input2%.dose*}
OUTPUT="stat_impute_${_input3}.txt"
printf "CHROM\tPOS\tID\tREF\tALT\t" > ${OUTPUT}
bcftools query -l ${INPUT} | tr "\n" "\t" >> ${OUTPUT}
gsed -E 's/\t$/\n/g' ${OUTPUT} > ${OUTPUT}.a
mv ${OUTPUT}.a ${OUTPUT}
bcftools view -i 'INFO/IMPUTED=1' ${INPUT} -Ou | \
    bcftools query -f '%CHROM\t%POS\t%CHROM:%POS:%REF:%ALT\t%REF\t%ALT[\t%DS]\n' >> ${OUTPUT}
# compress the file and delete the original file
gzip ${OUTPUT}
echo "${OUTPUT}.gz"
