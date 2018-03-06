#!/bin/bash

VCFPATH=$1
TMPPATH=$2
FAMFILE=$3
CHROM=$4

rm -f ${TMPPATH}/mend_stats_${CHROM}.tab
touch ${TMPPATH}/mend_stats_${CHROM}.tab

samples=$(cat ${FAMFILE}  | cut -f 2)
families=$(cat ${FAMFILE} | cut -d'.' -f 1 | uniq)

# TODO do on comet with filtered VCFs
vcf=${VCFPATH}/chr${CHROM}/hipstr.chr${CHROM}.with.1kg.filtered.vcf.gz

# Get Mendelian inheritance
echo "Getting Mend stats ${CHROM}"
for fam in $families
do
    child=$(cat ${FAMFILE} | grep "$fam\.s1" | cut -f 2)
    proband=$(cat ${FAMFILE} | grep "$fam\.p1" | cut -f 2)
    mother=$(cat ${FAMFILE} | grep "$fam\.mo" | cut -f 2)
    father=$(cat ${FAMFILE} | grep "$fam\.fa" | cut -f 2)
    echo ${child},${mother},${father}
    bcftools query -s ${child},${mother},${father} \
	-f "%CHROM\t%INFO/START[\t%SAMPLE\t%GB\t%Q\t%DP]\n" ${vcf} | \
	./get_mend.py >> ${TMPPATH}/mend_stats_${CHROM}.tab
    echo ${proband},${mother},${father}
    bcftools query -s ${proband},${mother},${father} \
	-f "%CHROM\t%INFO/START[\t%SAMPLE\t%GB\t%Q\t%DP]\n" ${vcf} | \
	./get_mend.py >> ${TMPPATH}/mend_stats_${CHROM}.tab
done
