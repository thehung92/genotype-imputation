#!/bin/bash
# create symlink of HGDP data to input dir
ln -s /home/hungntt/area7/hung_dir/reference_data/HGDP/HGDP.b38_aligned.vcf.gz* ./input/

# >>>>>>>>>>>>>>>>>>>>>>>> process input array >>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<< process input array <<<<<<<<<<<<<<<<<<<<<<<<
# prepare the array data
# remove multi-allelic variants
INPUT="./input/HGDP.b38_aligned.vcf.gz"
OUTPUT="./temp/temp-input/HGDP.b38_QC1.vcf.gz"
bcftools view -M 2 $INPUT -Oz -o $OUTPUT
# extract chr20
INPUT="./temp/temp-input/HGDP.b38_QC1.vcf.gz"
bcftools index $INPUT
OUTPUT="./temp/temp-input/HGDP.chr20_QC1.vcf.gz"
bcftools view -r chr20 $INPUT -Oz -o $OUTPUT
# generate masked.snp.list using R
INPUT="./temp/temp-input/HGDP.chr20_QC1.vcf.gz"
Rscript ./src/mask-snp-scr.R $INPUT ./temp/
# masked snp in vcf file using snp.list
INPUT="./temp/temp-input/HGDP.chr20_QC1.vcf.gz"
ARG="./temp/masked.snps.list"
OUTPUT="./temp/temp-input/HGDP.chr20_QC2.vcf.gz"
bcftools view -T ^$ARG $INPUT -Oz -o $OUTPUT
bcftools index $OUTPUT

# >>>>>>>>>>>>>>>>>>>>>>>> simulate kinh94 to hgdp_array from sequencing data >>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<< simulate kinh94 to hgdp_array from sequencing data <<<<<<<<<<<<<<<<<<<<<<<<
# create hgdp type array [Illumina650k] from vn94 sequencing data [with split multi-allelic]
INPUT="./data/test1014.chr20_QC.vcf.gz"
ARG1="./temp/samples-rename.txt"
ARG2="./temp/list-kinh94.txt"
awk '{print $2}' $ARG1 > $ARG2
OUTPUT="./input/Kinh94.chr20_QC.vcf.gz"
bcftools view -S $ARG2 $INPUT -Oz -o $OUTPUT
# subset with bcftools isec based on chrom:pos:ref:alt
INPUT="./input/Kinh94.chr20_QC.vcf.gz"
# bcftools index $INPUT
ARG="./temp/temp-input/HGDP.chr20_QC1.vcf.gz"
bcftools index $ARG
OUTDIR="./temp/temp-isec"
bcftools isec -c none -n~11 $INPUT $ARG -p $OUTDIR -Oz -w 1
OUTPUT="./temp/temp-isec/0000.vcf.gz"
# merge with HGDP file assuming missing value is 0/0
INPUT1="./temp/temp-isec/0000.vcf.gz"
INPUT2="./temp/temp-input/HGDP.chr20_QC1.vcf.gz"
OUTPUT="./input/HGDP-vn94.chr20_QC1.vcf.gz"
bcftools merge --missing-to-ref -m none $INPUT1 $INPUT2 -Oz -o $OUTPUT
[[ -f  ${OUTPUT}.csi ]] && bcftools index $OUTPUT
# mask snp
INPUT="./input/HGDP-vn94.chr20_QC1.vcf.gz"
ARG="./temp/masked.snps.list"
OUTPUT="./temp/temp-input/HGDP-vn94.chr20_QC2.vcf.gz"
bcftools view -T ^$ARG $INPUT -Oz -o $OUTPUT
bcftools index $OUTPUT

# >>>>>>>>>>>>>>>>>>>>>>>> prephase with corresponding ref panel >>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<< prephase with corresponding ref panel <<<<<<<<<<<<<<<<<<<<<<<<
# function for prephase
_prephase () {
    INPUT=$1
    ARG1=`basename $INPUT | grep -o 'chr[0-9X]*'`
    ARG2="/home/hungntt/area7/hung_dir/reference_data/GeneticMap/GRCh38/shapeit_geneticmap_${ARG1}_hg38.txt"
    ARG3=$2
    OUTPUT="./temp/temp-input/`basename $INPUT | sed 's/_.*/_phase./'``basename $ARG3 | sed 's/.chr.*//'`.vcf.gz"; echo $OUTPUT
    shapeit4 --input $INPUT --region $ARG1 --map $ARG2 \
        --reference $ARG3 --output $OUTPUT --log ${OUTPUT%.vcf*}.log --thread 8
}
# run prephase function with multiple ref panel
_prephase ./temp/temp-input/HGDP-vn94.chr20_QC2.vcf.gz ./temp/temp-data/1KGP3.chr20_rmcentro.vcf.gz &
_prephase ./temp/temp-input/HGDP-vn94.chr20_QC2.vcf.gz "./temp/temp-data/test920.chr20_rmcentro.vcf.gz" &
_prephase ./temp/temp-input/HGDP-vn94.chr20_QC2.vcf.gz "./temp/temp-data/merge-1KGP3-test920.chr20_rmcentro.vcf.gz" &
_prephase ./temp/temp-input/HGDP-vn94.chr20_QC2.vcf.gz "./temp/temp-data/merge-reciprocal.chr20_QC.vcf.gz" &
_prephase ./temp/temp-input/HGDP-vn94.chr20_QC2.vcf.gz "./temp/temp-data/SG10K.chr20_QC3.vcf.gz" &
# run prephase function with HGDP-vn94
# 1KGP
INPUT="./temp/temp-input/HGDP.chr20_QC2.vcf.gz"
ARG1=`basename $OUTPUT | grep -o 'chr[0-9X]*'`
ARG2="/home/hungntt/area7/hung_dir/reference_data/GeneticMap/GRCh38/shapeit_geneticmap_${ARG1}_hg38.txt"
ARG3="./temp/temp-data/1KGP3.chr20_rmcentro.vcf.gz"
OUTPUT="./temp/temp-input/`basename $INPUT | sed 's/_.*/_phase./'``basename $ARG3 | sed 's/.chr.*//'`.vcf.gz"; echo $OUTPUT
shapeit4 --input $INPUT --region $ARG1 --map $ARG2 \
    --reference $ARG3 --output $OUTPUT --log ${OUTPUT%.vcf*}.log --thread 4 &
# test920
INPUT="./temp/temp-input/HGDP.chr20_QC2.vcf.gz"
ARG1=`basename $OUTPUT | grep -o 'chr[0-9X]*'`
ARG2="/home/hungntt/area7/hung_dir/reference_data/GeneticMap/GRCh38/shapeit_geneticmap_${ARG1}_hg38.txt"
ARG3="./temp/temp-data/test920.chr20_rmcentro.vcf.gz"
OUTPUT="./temp/temp-input/`basename $INPUT | sed 's/_.*/_phase./'``basename $ARG3 | sed 's/.chr.*//'`.vcf.gz" ; echo $OUTPUT
shapeit4 --input $INPUT --region $ARG1 --map $ARG2 \
    --reference $ARG3 --output $OUTPUT --log ${OUTPUT%.vcf*}.log --thread 4 &
# merge
INPUT="./temp/temp-input/HGDP.chr20_QC2.vcf.gz"
ARG1=`basename $OUTPUT | grep -o 'chr[0-9X]*'`
ARG2="/home/hungntt/area7/hung_dir/reference_data/GeneticMap/GRCh38/shapeit_geneticmap_${ARG1}_hg38.txt"
ARG3="./temp/temp-data/merge-1KGP3-test920.chr20_rmcentro.vcf.gz"
OUTPUT="./temp/temp-input/`basename $INPUT | sed 's/_.*/_phase./'``basename $ARG3 | sed 's/.chr.*//'`.vcf.gz" ; echo $OUTPUT
shapeit4 --input $INPUT --region $ARG1 --map $ARG2 \
    --reference $ARG3 --output $OUTPUT --log ${OUTPUT%.vcf*}.log --thread 4 &
# SG10K
INPUT="./temp/temp-input/HGDP.chr20_QC2.vcf.gz"
ARG1=`basename $OUTPUT | grep -o 'chr[0-9X]*'`
ARG2="/home/hungntt/area7/hung_dir/reference_data/GeneticMap/GRCh38/shapeit_geneticmap_${ARG1}_hg38.txt"
ARG3="./temp/temp-data/SG10K.chr20_QC3.vcf.gz"
OUTPUT="./temp/temp-input/`basename $INPUT | sed 's/_.*/_phase./'``basename $ARG3 | sed 's/.chr.*//'`.vcf.gz" ; echo $OUTPUT
shapeit4 --input $INPUT --region $ARG1 --map $ARG2 \
    --reference $ARG3 --output $OUTPUT --log ${OUTPUT%.vcf*}.log --thread 8 &
# merge-reciprocal
INPUT="./temp/temp-input/HGDP.chr20_QC2.vcf.gz"
ARG1=`basename $OUTPUT | grep -o 'chr[0-9X]*'`
ARG2="/home/hungntt/area7/hung_dir/reference_data/GeneticMap/GRCh38/shapeit_geneticmap_${ARG1}_hg38.txt"
ARG3="./temp/temp-data/merge-reciprocal.chr20_QC.vcf.gz"
OUTPUT="./temp/temp-input/`basename $INPUT | sed 's/_.*/_phase./'``basename $ARG3 | sed 's/.chr.*//'`.vcf.gz" ; echo $OUTPUT
shapeit4 --input $INPUT --region $ARG1 --map $ARG2 \
    --reference $ARG3 --output $OUTPUT --log ${OUTPUT%.vcf*}.log --thread 8 &

# >>>>>>>>>>>>>>>>>>>>>>>> impute with minimac4 >>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<< impute with minimac4 <<<<<<<<<<<<<<<<<<<<<<<<
# function
_impute () {
    INPUT=$1
    ARG1=`basename $INPUT | sed 's/.*phase.\(.*\).vcf.*/\1/'` ; echo $ARG1
    ARG2=`ls ./data/${ARG1}*"m3vcf.gz"` ; echo $ARG2
    OUTDIR=$2
    OUTPUT=${OUTDIR}`basename $INPUT | sed 's/phase/impute/; s/.vcf.gz//'` ; echo $OUTPUT
    minimac4 --haps $INPUT --refHaps $ARG2 \
        --ChunkLengthMb 20 --ChunkOverlapMb 3 --allTypedSites \
        --prefix $OUTPUT --log --cpus 8
}
export -f _impute
#
ARR=(`ls ./temp/temp-input/HGDP*phase*gz`)
parallel _impute {1} ./output/ ::: ${ARR[@]}
_impute "./temp/temp-input/HGDP-vn94.chr20_phase.1KGP3.vcf.gz" ./output/
#
OUTDIR="./output/"
# 1KGP
INPUT="./temp/temp-input/HGDP.chr20_phase.1KGP3.vcf.gz"
ARG1=`basename $INPUT | sed 's/.*phase.\(.*\).vcf.*/\1/'` ; echo $ARG1
ARG2=`ls ./data/${ARG1}*"m3vcf.gz"` ; echo $ARG2
OUTPUT=${OUTDIR}`basename $INPUT | sed 's/phase/impute/; s/.vcf.gz//'` ; echo $OUTPUT
minimac4 --haps $INPUT --refHaps $ARG2 \
    --ChunkLengthMb 20 --ChunkOverlapMb 3 --allTypedSites \
    --prefix $OUTPUT --log --cpus 8 &
# test920
INPUT="./temp/temp-input/HGDP.chr20_phase.test920.vcf.gz"
ARG1=`basename $INPUT | sed 's/.*phase.\(.*\).vcf.*/\1/'`
ARG2=`ls ./data/${ARG1}*"m3vcf.gz"` ; echo $ARG2
OUTPUT=${OUTDIR}`basename $INPUT | sed 's/phase/impute/; s/.vcf.gz//'` ; echo $OUTPUT
minimac4 --haps $INPUT --refHaps $ARG2 \
    --ChunkLengthMb 20 --ChunkOverlapMb 3 --allTypedSites \
    --prefix $OUTPUT --log --cpus 8 &
# merge
INPUT="./temp/temp-input/HGDP.chr20_phase.merge-1KGP3-test920.vcf.gz"
ARG1=`basename $INPUT | sed 's/.*phase.\(.*\).vcf.*/\1/'`
ARG2=`ls ./data/${ARG1}*"m3vcf.gz"` ; echo $ARG2
OUTPUT=${OUTDIR}`basename $INPUT | sed 's/phase/impute/; s/.vcf.gz//'` ; echo $OUTPUT
minimac4 --haps $INPUT --refHaps $ARG2 \
    --ChunkLengthMb 20 --ChunkOverlapMb 3 --allTypedSites \
    --prefix $OUTPUT --log --cpus 8 &
# SG10K
INPUT="./temp/temp-input/HGDP.chr20_phase.SG10K.vcf.gz"
ARG1=`basename $INPUT | sed 's/.*phase.\(.*\).vcf.*/\1/'`
ARG2=`ls ./data/${ARG1}*"m3vcf.gz"` ; echo $ARG2
OUTPUT=${OUTDIR}`basename $INPUT | sed 's/phase/impute/; s/.vcf.gz//'` ; echo $OUTPUT
minimac4 --haps $INPUT --refHaps $ARG2 \
    --ChunkLengthMb 20 --ChunkOverlapMb 3 --allTypedSites \
    --prefix $OUTPUT --log --cpus 8 &
# merge-reciprocal
INPUT="./temp/temp-input/HGDP.chr20_phase.merge-reciprocal.vcf.gz"
ARG1=`basename $INPUT | sed 's/.*phase.\(.*\).vcf.*/\1/'`
ARG2=`ls ./data/${ARG1}*"m3vcf.gz"` ; echo $ARG2
OUTPUT=${OUTDIR}`basename $INPUT | sed 's/phase/impute/; s/.vcf.gz//'` ; echo $OUTPUT
minimac4 --haps $INPUT --refHaps $ARG2 \
    --ChunkLengthMb 20 --ChunkOverlapMb 3 --allTypedSites \
    --prefix $OUTPUT --log --cpus 8 &