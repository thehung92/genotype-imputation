#!/bin/bash
# setup docker container
docker run --name imputation \
    -v /dragennfs/area7/:/home/public/area7/ \
    -v /dragennfs/area16/:/home/public/area16/ \
    -it registry.vinbdi.org/imputation:v5.2 bash
# cd /home/hungntt/area7/hung_dir/reference_data/DGV4VN
# mkdir -p /home/hungntt/area7/hung_dir/reference_data/DGV4VN/dragen-output
OUTDIR="/home/hungntt/area7/hung_dir/reference_data/DGV4VN/dragen-output"
cd $OUTDIR
INPUT="/home/hungntt/area16/Population/testPop_1030.jointGenotype.hard-filtered.vcf.gz"
OUTPUT="test1030_filter.vcf.gz"
bcftools view -f "PASS" $INPUT |\
  bcftools annotate -x ^FORMAT/GQ,FORMAT/GT --threads 10 -Oz -o $OUTPUT
bcftools index $OUTPUT
echo "done"
# subset files based on list of samples from '/dragennfs/area16/Population/testPop_1014.pass.vcf.gz' and remove GQ for phasing
bcftools query -l /home/hungntt/area16/Population/testPop_1014.pass.vcf.gz > samples_list.txt
INPUT="test1030_filter.vcf.gz"
OUTPUT="test1014_filter.vcf.gz"
bcftools view -S samples_list.txt $INPUT |\
  bcftools annotate -x ^FORMAT/GT -Oz -o $OUTPUT
bcftools index $OUTPUT
#  split chromosome
INPUT="test1014_filter.vcf.gz"
for i in {{1..22},X}; do
  OUTPUT=`basename $INPUT | sed "s/_/.chr${i}_/"`
  bcftools view -r chr${i} $INPUT -Oz -o $OUTPUT &
done

#### upload data to hpc ####
# via sftp
put *QC1* /dragennfs/area7/hung_dir/reference_data/1KGP_LC/
