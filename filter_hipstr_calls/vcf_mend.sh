#!/bin/bash

#SBATCH -A csd568
#SBATCH --get-user-env
#SBATCH -p shared
#SBATCH --mem=3G
#SBATCH -t 2000
#SBATCH --job-name mend

source params.sh

CHROM=$1

samples=$(cat ${FAMFILE}  | cut -f 2)
families=$(cat ${FAMFILE} | cut -d'.' -f 1 | uniq)

vcf=${FILTEREDVCF}/hipstr.chr${CHROM}.with.1kg.calls_filtered.vcf.gz
outfile=${SAMPDIR}/mend_stats_${CHROM}.tab

echo "chrom,start,sample,Q,DP,mend,homref" | sed 's/,/\t/g' > ${outfile}

for fam in $families
do
    child=$(cat ${FAMFILE} | grep "$fam\.s1" | cut -f 2)
    proband=$(cat ${FAMFILE} | grep "$fam\.p1" | cut -f 2)
    mother=$(cat ${FAMFILE} | grep "$fam\.mo" | cut -f 2)
    father=$(cat ${FAMFILE} | grep "$fam\.fa" | cut -f 2)
    echo ${child},${mother},${father}
    bcftools query -s ${child},${mother},${father} \
	-f "%CHROM\t%INFO/START[\t%SAMPLE\t%GB\t%Q\t%DP]\n" ${vcf} | \
	./get_mend.py >> ${outfile}
    echo ${proband},${mother},${father}
    bcftools query -s ${proband},${mother},${father} \
	-f "%CHROM\t%INFO/START[\t%SAMPLE\t%GB\t%Q\t%DP]\n" ${vcf} | \
	./get_mend.py >> ${outfile}

done