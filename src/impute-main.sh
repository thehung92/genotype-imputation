#!/bin/bash
# working directory inside container
cd /home/hungntt/area7/hung_dir/run_4
# create dir structure
mkdir data docs output src input


# create symlink to input inside ./input
ln -s /home/hungntt/area7/hung_dir/reference_data/APMRA/APMRA96_annotate.vcf.* ./input/

# >>>>>>>>>>>>>>>>>>>>>>>> process input array >>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<< process input array <<<<<<<<<<<<<<<<<<<<<<<<
# process input array
# exclude multi-allelic variants
INPUT="./input/APMRA96_annotate.vcf.gz"
# mkdir -p ./temp/temp-input/
TEMP="./temp/temp-input/dup-pos.txt"
OUTPUT="./temp/temp-input/APMRA96_QC1.vcf.gz"
OUTPUT2="./input/APMRA96.chr20_QC1.vcf.gz"
Rscript ./src/dup-pos-scr.R $INPUT $TEMP
bcftools view -T ^$TEMP $INPUT -Oz -o $OUTPUT
bcftools index $OUTPUT
# extract chr20
bcftools view -r chr20 $OUTPUT -Oz -o $OUTPUT2

# >>>>>>>>>>>>>>>>>>>>>>>> extract data, symlink to wd >>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<< extract data, symlink to wd <<<<<<<<<<<<<<<<<<<<<<<<
# create symlink to data inside ./data
ln -s /home/hungntt/area7/hung_dir/reference_data/DGV4VN/dragen-output/test1014.chr20_filter.vcf.gz* ./data/
ln -s /home/hungntt/area7/hung_dir/reference_data/1KGP_LC/1KGP3.chr20_QC1.vcf.gz* ./data/
ln -s /home/hungntt/area7/hung_dir/reference_data/SG10K/hg38/SG10K.allchr.hg38_sort.vcf.gz* ./data
# subset sample to create test and ref data using Rscript
ARG1="./data/test1014.chr20_filter.vcf.gz"
ARG2="./input/APMRA96_annotate.vcf.gz"
ARG3="./temp/"
Rscript ./src/filter-sample-scr.R $ARG1 $ARG2 $ARG3
# output: list-ref-samples.txt, list-input-samples.txt, samples-rename.txt
# extract chr20 from SG10K
INPUT="./data/SG10K.allchr.hg38_sort.vcf.gz"
ARG="chr20"
OUTPUT="./data/SG10K.chr20_hg38.vcf.gz"
bcftools view -r $ARG $INPUT -Oz -o $OUTPUT

# >>>>>>>>>>>>>>>>>>>>>>>> input array >>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<< input array <<<<<<<<<<<<<<<<<<<<<<<<
# process input data, subset based on list-input-samples.txt, and rename based on samples-rename.txt
INPUT="./input/APMRA96.chr20_QC1.vcf.gz"
ARG1="./temp/list-input-samples.txt"
ARG2="./temp/samples-rename.txt"
OUTPUT="./temp/temp-input/APMRA94.chr20_QC2.vcf.gz"
bcftools view -S $ARG1 $INPUT -Ov |\
  bcftools reheader -s $ARG2 | bgzip -c > $OUTPUT
# check:  bcftools query -l $OUTPUT | wc -l

# >>>>>>>>>>>>>>>>>>>>>>>> phase VN with shapeit4 >>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<< phase VN with shapeit4 <<<<<<<<<<<<<<<<<<<<<<<<
# process VN reference panel, in the same principle as 1KGP
# split multi-allelic variants, remove AC=0, and phase with shapeit4
INPUT="./data/test1014.chr20_filter.vcf.gz"
OUTPUT="./data/test1014.chr20_QC.vcf.gz"
bcftools norm -m- $INPUT -Ou |\
    bcftools view -i 'INFO/AC>0' -Oz -o $OUTPUT
bcftools index $OUTPUT
OUTPUT2="./data/test1014.chr20_phase.vcf.gz"
ARG1=`basename $OUTPUT | grep -o 'chr[0-9X]*'`
ARG2="/home/hungntt/area7/hung_dir/reference_data/GeneticMap/GRCh38/shapeit_geneticmap_${ARG1}_hg38.txt"
shapeit4 --input $OUTPUT --sequencing --region $ARG1 --map $ARG2 \
    --output $OUTPUT2 --log ${OUTPUT2%.vcf*}.log --thread 8 &&
bcftools index $OUTPUT2 

# >>>>>>>>>>>>>>>>>>>>>>>> Process ref panel >>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<< Process ref panel <<<<<<<<<<<<<<<<<<<<<<<<
# process SG10K panel, phased sequencing data
# merge multi-alleleic, then rm them to get only bialllelic variants
# rm AC<=1 because this do not require merging, output _QC2
INPUT="./data/SG10K.chr20_hg38.vcf.gz"
OUTPUT="./temp/temp-data/SG10K.chr20_QC2.vcf.gz"
bcftools norm -m+any $INPUT -Ou |\
    bcftools view -v 'snps,indels' -M 2 -Ou |\
    bcftools view -i 'INFO/AC>1' -Ou |\
    bcftools +fill-tags -Oz -o ${OUTPUT} -- -t MAF
# due to liftover, check file and exclude GT haploid
INPUT="./temp/temp-data/SG10K.chr20_QC2.vcf.gz"
OUTPUT="./temp/temp-data/SG10K.chr20_QC3.vcf.gz"
bcftools view -e 'GT="hap"' $INPUT -Oz -o $OUTPUT


# subset based on list-ref-samples.txt
INPUT="./data/test1014.chr20_phase.vcf.gz"
ARG="./temp/list-ref-samples.txt"
mkdir -p ./temp/temp-data/
OUTPUT="./temp/temp-data/test920.chr20_phase.vcf.gz"
bcftools view -S $ARG $INPUT -Oz -o $OUTPUT
# QC step 1: remove multi-allelic; rm AC=0; fill MAF tags
OUTPUT2="./temp/temp-data/test920.chr20_QC1.vcf.gz"
bcftools norm -m+any $OUTPUT -Ou |\
  bcftools view -v 'snps,indels' -M 2 -Ou |\
  bcftools +fill-tags -Oz -o $OUTPUT2 -- -t MAF
bcftools index $OUTPUT2
# check:  bcftools query -l $OUTPUT | wc -l

# >>>>>>>>>>>>>>>>>>>>>>>> rm uniq AC=1 before merging >>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<< rm uniq AC=1 before merging <<<<<<<<<<<<<<<<<<<<<<<<
# use Rscript to output list of uniq variants with AC=1
ARG1="./temp/temp-data/test920.chr20_QC1.vcf.gz"
ARG2="./data/1KGP3.chr20_QC1.vcf.gz"
ARG3="./temp/"
Rscript ./src/remove-uniq-variants.R $ARG1 $ARG2 $ARG3
# output: 1KGP3.chr20_uniqAC1.txt       test920.chr20_uniqAC1.txt

# subset to exclude uniqAC1 from VN panel and 1KGP
# test920
INPUT="./temp/temp-data/test920.chr20_QC1.vcf.gz"
ARG="./temp/test920.chr20_uniqAC1.txt"
OUTPUT="./temp/temp-data/test920.chr20_rmuniq.vcf.gz"
bcftools view -T ^$ARG $INPUT -Oz -o $OUTPUT &
# 1KGP
INPUT="./data/1KGP3.chr20_QC1.vcf.gz"
ARG="./temp/1KGP3.chr20_uniqAC1.txt"
OUTPUT="./temp/temp-data/1KGP3.chr20_rmuniq.vcf.gz"
bcftools view -T ^$ARG $INPUT -Oz -o $OUTPUT &
# bcftools view -GH $INPUT | wc -l &
# bcftools view -GH $OUTPUT | wc -l &
# SG10K do not require this step

# >>>>>>>>>>>>>>>>>>>>>>>> merge panel assuming missing=0/0 >>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<< merge panel assuming missing=0/0 <<<<<<<<<<<<<<<<<<<<<<<<
# merge vcf with bcftools merge and assume genotype at missing sites are 0/0
INPUT1="./temp/temp-data/test920.chr20_rmuniq.vcf.gz"
INPUT2="./temp/temp-data/1KGP3.chr20_rmuniq.vcf.gz"
OUTPUT="./temp/temp-data/merge-1KGP3-test920.chr20.vcf.gz"
# bcftools index $INPUT1 & bcftools index $INPUT2 &
bcftools merge --missing-to-ref -m none $INPUT1 $INPUT2 -Oz -o $OUTPUT
# fix missing format and update AN,AC,AF
OUTPUT2="./temp/temp-data/merge-1KGP3-test920.chr20_QC.vcf.gz"
zcat $OUTPUT | sed -e '/^#/!s/\//\|/g' |\
    bcftools +fill-tags -Oz -o $OUTPUT2 -- -t AN,AC,AF

# >>>>>>>>>>>>>>>>>>>>>>>> rm centromere before convert m3vcf >>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<< rm centromere before convert m3vcf <<<<<<<<<<<<<<<<<<<<<<<<
# remove centromere region due to difference in SNP density
# 1KGP
INPUT="./temp/temp-data/1KGP3.chr20_rmuniq.vcf.gz"
ARG="./data/centromere.chr20-hg38.bed"
OUTPUT=`echo $INPUT | sed 's/rmuniq/rmcentro/'` ; echo $OUTPUT
bcftools view -T ^$ARG $INPUT -Oz -o $OUTPUT &
# test920
INPUT="./temp/temp-data/test920.chr20_rmuniq.vcf.gz"
ARG="./data/centromere.chr20-hg38.bed"
OUTPUT=`echo $INPUT | sed 's/rmuniq/rmcentro/'` ; echo $OUTPUT
bcftools view -T ^$ARG $INPUT -Oz -o $OUTPUT &
# merge
INPUT="./temp/temp-data/merge-1KGP3-test920.chr20_QC.vcf.gz"
ARG="./data/centromere.chr20-hg38.bed"
OUTPUT=`echo $INPUT | sed 's/QC/rmcentro/'` ; echo $OUTPUT
bcftools view -T ^$ARG $INPUT -Oz -o $OUTPUT &
# SG10K
INPUT= [not yet run]
ARG="./data/centromere.chr20-hg38.bed"
OUTPUT=`echo $INPUT | sed 's/QC/rmcentro/'` ; echo $OUTPUT
bcftools view -T ^$ARG $INPUT -Oz -o $OUTPUT &

# >>>>>>>>>>>>>>>>>>>>>>>> m3vcf >>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<< m3vcf <<<<<<<<<<<<<<<<<<<<<<<<
# create m3vcf format for genotype imputation based on _rmcentro file
# 1KGP
INPUT="./temp/temp-data/1KGP3.chr20_rmcentro.vcf.gz"
bcftools index $INPUT
ARG=`basename $INPUT | grep -o 'chr[0-9X]*'`
OUTPUT="./data/`basename $INPUT | sed 's/_.*/_panel/'`" ; echo $OUTPUT
minimac3 --refHaps $INPUT --chr $ARG --processReference \
    --prefix $OUTPUT --log --cpus 8 &
# test920
INPUT="./temp/temp-data/test920.chr20_rmcentro.vcf.gz"
bcftools index $INPUT
ARG=`basename $INPUT | grep -o 'chr[0-9X]*'`
OUTPUT="./data/`basename $INPUT | sed 's/_.*/_panel/'`" ; echo $OUTPUT
minimac3 --refHaps $INPUT --chr $ARG --processReference \
    --prefix $OUTPUT --log --cpus 8 &
# merge
INPUT="./temp/temp-data/merge-1KGP3-test920.chr20_rmcentro.vcf.gz"
#bcftools index $INPUT
ARG=`basename $INPUT | grep -o 'chr[0-9X]*'`
OUTPUT="./data/`basename $INPUT | sed 's/_.*/_panel/'`" ; echo $OUTPUT
minimac3 --refHaps $INPUT --chr $ARG --processReference \
    --prefix $OUTPUT --log --cpus 8 &
# SG10K
INPUT="./temp/temp-data/SG10K.chr20_QC3.vcf.gz"
bcftools index $INPUT
ARG=`basename $INPUT | grep -o 'chr[0-9X]*'`
OUTPUT="./data/`basename $INPUT | sed 's/_.*/_panel/'`" ; echo $OUTPUT
minimac3 --refHaps $INPUT --chr $ARG --processReference \
    --prefix $OUTPUT --log --cpus 8 &

# >>>>>>>>>>>>>>>>>>>>>>>> merge panel, reciprocal impute >>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<< merge panel, reciprocal impute <<<<<<<<<<<<<<<<<<<<<<<<
# impute from vn920 to 1kgp, do not require pre-phase
INPUT="./temp/temp-data/test920.chr20_rmcentro.vcf.gz"
ARG="./data/1KGP3.chr20_panel.m3vcf.gz"
OUTPUT=`echo $INPUT | sed 's/_.*/_impute.1KGP3/'` ; echo $OUTPUT
minimac4 --haps $INPUT --refHaps $ARG \
    --ChunkLengthMb 20 --ChunkOverlapMb 3 --allTypedSites \
    --prefix $OUTPUT --log --cpus 8 &
# impute from 1kgp to vn920, do not require pre-phase
INPUT="./temp/temp-data/1KGP3.chr20_rmcentro.vcf.gz"
ARG="./data/test920.chr20_panel.m3vcf.gz"
OUTPUT=`echo $INPUT | sed 's/_.*/_impute.test920/'` ; echo $OUTPUT
minimac4 --haps $INPUT --refHaps $ARG \
    --ChunkLengthMb 20 --ChunkOverlapMb 3 --allTypedSites \
    --prefix $OUTPUT --log --cpus 8 &
# merge *dose.vcf.gz [merge multi-allelic]; then remove them from OUTPUT
INPUT="./temp/temp-data/test920.chr20_impute.1KGP3.dose.vcf.gz"
# bcftools index $INPUT
INPUT2="./temp/temp-data/1KGP3.chr20_impute.test920.dose.vcf.gz"
#bcftools index $INPUT2
OUTPUT="./temp/temp-data/merge-reciprocal.chr20_QC.vcf.gz"
bcftools merge --merge all $INPUT $INPUT2 -Ou |\
    bcftools view -M 2 -Oz -o $OUTPUT
bcftools index $OUTPUT
# >>>>>>>>>>>>>>>>>>>>>>>> m3vcf for reciprocal merge >>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<< m3vcf for reciprocal merge <<<<<<<<<<<<<<<<<<<<<<<<
INPUT="./temp/temp-data/merge-reciprocal.chr20_QC.vcf.gz"
# bcftools index $INPUT
ARG=`basename $INPUT | grep -o 'chr[0-9X]*'`
OUTPUT="./data/`basename $INPUT | sed 's/_.*/_panel/'`" ; echo $OUTPUT
minimac3 --refHaps $INPUT --chr $ARG --processReference \
    --prefix $OUTPUT --log --cpus 8 &


# >>>>>>>>>>>>>>>>>>>>>>>> prephase with corresponding ref panel >>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<< prephase with corresponding ref panel <<<<<<<<<<<<<<<<<<<<<<<<
# prephase array data with _rmcentro file
INPUT="./temp/temp-input/APMRA94.chr20_QC2.vcf.gz"
bcftools index $INPUT
# 1KGP
INPUT="./temp/temp-input/APMRA94.chr20_QC2.vcf.gz"
ARG1=`basename $OUTPUT | grep -o 'chr[0-9X]*'`
ARG2="/home/hungntt/area7/hung_dir/reference_data/GeneticMap/GRCh38/shapeit_geneticmap_${ARG1}_hg38.txt"
ARG3="./temp/temp-data/1KGP3.chr20_rmcentro.vcf.gz"
OUTPUT="./temp/temp-input/`basename $INPUT | sed 's/_.*/_phase./'``basename $ARG3 | sed 's/.chr.*//'`.vcf.gz"; echo $OUTPUT
shapeit4 --input $INPUT --region $ARG1 --map $ARG2 \
    --reference $ARG3 --output $OUTPUT --log ${OUTPUT%.vcf*}.log --thread 4 &
# test920
INPUT="./temp/temp-input/APMRA94.chr20_QC2.vcf.gz"
ARG1=`basename $OUTPUT | grep -o 'chr[0-9X]*'`
ARG2="/home/hungntt/area7/hung_dir/reference_data/GeneticMap/GRCh38/shapeit_geneticmap_${ARG1}_hg38.txt"
ARG3="./temp/temp-data/test920.chr20_rmcentro.vcf.gz"
OUTPUT="./temp/temp-input/`basename $INPUT | sed 's/_.*/_phase./'``basename $ARG3 | sed 's/.chr.*//'`.vcf.gz" ; echo $OUTPUT
shapeit4 --input $INPUT --region $ARG1 --map $ARG2 \
    --reference $ARG3 --output $OUTPUT --log ${OUTPUT%.vcf*}.log --thread 4 &
# merge
INPUT="./temp/temp-input/APMRA94.chr20_QC2.vcf.gz"
ARG1=`basename $OUTPUT | grep -o 'chr[0-9X]*'`
ARG2="/home/hungntt/area7/hung_dir/reference_data/GeneticMap/GRCh38/shapeit_geneticmap_${ARG1}_hg38.txt"
ARG3="./temp/temp-data/merge-1KGP3-test920.chr20_rmcentro.vcf.gz"
OUTPUT="./temp/temp-input/`basename $INPUT | sed 's/_.*/_phase./'``basename $ARG3 | sed 's/.chr.*//'`.vcf.gz" ; echo $OUTPUT
shapeit4 --input $INPUT --region $ARG1 --map $ARG2 \
    --reference $ARG3 --output $OUTPUT --log ${OUTPUT%.vcf*}.log --thread 4 &
# SG10K
INPUT="./temp/temp-input/APMRA94.chr20_QC2.vcf.gz"
ARG1=`basename $OUTPUT | grep -o 'chr[0-9X]*'`
ARG2="/home/hungntt/area7/hung_dir/reference_data/GeneticMap/GRCh38/shapeit_geneticmap_${ARG1}_hg38.txt"
ARG3="./temp/temp-data/SG10K.chr20_QC3.vcf.gz"
OUTPUT="./temp/temp-input/`basename $INPUT | sed 's/_.*/_phase./'``basename $ARG3 | sed 's/.chr.*//'`.vcf.gz" ; echo $OUTPUT
shapeit4 --input $INPUT --region $ARG1 --map $ARG2 \
    --reference $ARG3 --output $OUTPUT --log ${OUTPUT%.vcf*}.log --thread 8 &
# merge-reciprocal
INPUT="./temp/temp-input/APMRA94.chr20_QC2.vcf.gz"
ARG1=`basename $OUTPUT | grep -o 'chr[0-9X]*'`
ARG2="/home/hungntt/area7/hung_dir/reference_data/GeneticMap/GRCh38/shapeit_geneticmap_${ARG1}_hg38.txt"
ARG3="./temp/temp-data/merge-reciprocal.chr20_QC.vcf.gz"
OUTPUT="./temp/temp-input/`basename $INPUT | sed 's/_.*/_phase./'``basename $ARG3 | sed 's/.chr.*//'`.vcf.gz" ; echo $OUTPUT
shapeit4 --input $INPUT --region $ARG1 --map $ARG2 \
    --reference $ARG3 --output $OUTPUT --log ${OUTPUT%.vcf*}.log --thread 8 &

# >>>>>>>>>>>>>>>>>>>>>>>> impute with minimac4 >>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<< impute with minimac4 <<<<<<<<<<<<<<<<<<<<<<<<
#
OUTDIR="./output/"
# 1KGP
INPUT="./temp/temp-input/APMRA94.chr20_phase.1KGP3.vcf.gz"
ARG1=`basename $INPUT | sed 's/.*phase.\(.*\).vcf.*/\1/'` ; echo $ARG1
ARG2=`ls ./data/${ARG1}*"m3vcf.gz"` ; echo $ARG2
# ARG="./data/1KGP3.chr20_panel.m3vcf.gz" ; echo $ARG
OUTPUT=${OUTDIR}`basename $INPUT | sed 's/phase/impute/; s/.vcf.gz//'` ; echo $OUTPUT
minimac4 --haps $INPUT --refHaps $ARG2 \
    --ChunkLengthMb 20 --ChunkOverlapMb 3 --allTypedSites \
    --prefix $OUTPUT --log --cpus 8 &
# test920
INPUT="./temp/temp-input/APMRA94.chr20_phase.test920.vcf.gz"
ARG1=`basename $INPUT | sed 's/.*phase.\(.*\).vcf.*/\1/'`
ARG2=`ls ./data/${ARG1}*"m3vcf.gz"` ; echo $ARG2
OUTPUT=${OUTDIR}`basename $INPUT | sed 's/phase/impute/; s/.vcf.gz//'` ; echo $OUTPUT
minimac4 --haps $INPUT --refHaps $ARG2 \
    --ChunkLengthMb 20 --ChunkOverlapMb 3 --allTypedSites \
    --prefix $OUTPUT --log --cpus 8 &
# merge
INPUT="./temp/temp-input/APMRA94.chr20_phase.merge-1KGP3-test920.vcf.gz"
ARG1=`basename $INPUT | sed 's/.*phase.\(.*\).vcf.*/\1/'`
ARG2=`ls ./data/${ARG1}*"m3vcf.gz"` ; echo $ARG2
OUTPUT=${OUTDIR}`basename $INPUT | sed 's/phase/impute/; s/.vcf.gz//'` ; echo $OUTPUT
minimac4 --haps $INPUT --refHaps $ARG2 \
    --ChunkLengthMb 20 --ChunkOverlapMb 3 --allTypedSites \
    --prefix $OUTPUT --log --cpus 8 &
# SG10K
INPUT="./temp/temp-input/APMRA94.chr20_phase.SG10K.vcf.gz"
ARG1=`basename $INPUT | sed 's/.*phase.\(.*\).vcf.*/\1/'`
ARG2=`ls ./data/${ARG1}*"m3vcf.gz"` ; echo $ARG2
OUTPUT=${OUTDIR}`basename $INPUT | sed 's/phase/impute/; s/.vcf.gz//'` ; echo $OUTPUT
minimac4 --haps $INPUT --refHaps $ARG2 \
    --ChunkLengthMb 20 --ChunkOverlapMb 3 --allTypedSites \
    --prefix $OUTPUT --log --cpus 8 &
# merge-reciprocal
# wait `jobs -p`
INPUT="./temp/temp-input/APMRA94.chr20_phase.merge-reciprocal.vcf.gz"
ARG1=`basename $INPUT | sed 's/.*phase.\(.*\).vcf.*/\1/'`
ARG2=`ls ./data/${ARG1}*"m3vcf.gz"` ; echo $ARG2
OUTPUT=${OUTDIR}`basename $INPUT | sed 's/phase/impute/; s/.vcf.gz//'` ; echo $OUTPUT
minimac4 --haps $INPUT --refHaps $ARG2 \
    --ChunkLengthMb 20 --ChunkOverlapMb 3 --allTypedSites \
    --prefix $OUTPUT --log --cpus 8 &