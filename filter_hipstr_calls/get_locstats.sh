#!/bin/bash

source params.sh

if [ "x${SLURM_ARRAY_TASK_ID}" == "x" ]
then
    CHROM=$1
else
    CHROM=${SLURM_ARRAY_TASK_ID}
fi

invcf=${FILTEREDVCF}/hipstr.chr${CHROM}.with.1kg.calls_filtered.vcf.gz
outfile=${FILTERDIR}/hipstr_filtered_locstats_chr${CHROM}.tab
#outfile=/oasis/scratch/comet/mgymrek/temp_project/hipstr_filtered_locstats_chr${CHROM}.tab

python get_locstats.py \
    --vcf ${invcf} \
    --unrelated-samples ${PARENTS} \
    > ${outfile}
