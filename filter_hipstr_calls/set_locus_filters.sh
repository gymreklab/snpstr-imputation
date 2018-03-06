#!/bin/bash

source params.sh

if [ "x${SLURM_ARRAY_TASK_ID}" == "x" ]
then
    CHROM=$1
else
    CHROM=${SLURM_ARRAY_TASK_ID}
fi

# Combine all filter files to a single one
FILTERFILE=${FILTERDIR}/hipstr_combined_info_chr${CHROM}.tab
echo "chrom,start,end,period,hrun,segdup,hwe.p,numcalls,het,mean_allele" | sed 's/,/\t/g' > ${FILTERFILE}
cat ${HRUN} | sed 's/chr//' | intersectBed -a stdin -b ${SEGDUP} -c | grep -w "^${CHROM}" | \
    intersectBed -a stdin -b ${FILTERDIR}/hipstr_filtered_locstats_chr${CHROM}.tab -wa -wb -f 1 | \
    cut -f 7-9 --complement >> ${FILTERFILE}

# Run filter script
./add_filters.py \
    --vcf ${FILTEREDVCF}/hipstr.chr${CHROM}.with.1kg.calls_filtered.vcf.gz \
    --statsfile ${FILTERFILE} \
    --out ${FINALVCFS}/hipstr.chr${CHROM}.allfilters \
    --min-hwep 0.01 \
    --min-callrate 0.8 \
    --min-het 0.095 \
    --max-hrun-offset -1 \
    --filter-segdup | bgzip -c > ${FINALVCFS}/hipstr.chr${CHROM}.allfilters.vcf.gz
tabix -p vcf -f ${FINALVCFS}/hipstr.chr${CHROM}.allfilters.vcf.gz