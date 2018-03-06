#!/bin/bash

source params.sh

if [ "x${SLURM_ARRAY_TASK_ID}" == "x" ]
then
    CHROM=$1
else
    CHROM=${SLURM_ARRAY_TASK_ID}
fi

invcf=${HIPSTRVCF}/hipstr.chr${CHROM}.with.1kg.filtered.vcf.gz
outvcf=${FILTEREDVCF}/hipstr.chr${CHROM}.with.1kg.calls_filtered.vcf.gz

python /home/mgymrek/workspace/HipSTR-main/scripts/filter_vcf.py \
    --vcf ${invcf} \
    --min-call-qual 0.9 \
    --max-call-flank-indel 0.15 \
    --max-call-stutter 0.15 \
    | bgzip -c> ${outvcf}

tabix -f -p vcf ${outvcf}