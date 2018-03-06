#!/bin/bash

sample=$1

echo ${sample}

python write_bases.py imputed.str.${sample}.vcf.gz $sample $sample.imputeResult.txt
python write_bases.py hipstr.chr.phased.parents.vcf.gz $sample $sample.groundTruth.txt

sort -f ${sample}.groundTruth.txt > ${sample}.ground.sorted.txt
sort -f ${sample}.imputeResult.txt > ${sample}.imputed.sorted.txt
join -1 1 -2 1 ${sample}.ground.sorted.txt ${sample}.imputed.sorted.txt | awk 'NF==5{print}' > ${sample}.diff.txt
rm ${sample}.groundTruth.txt ${sample}.imputeResult.txt ${sample}.ground.sorted.txt ${sample}.imputed.sorted.txt

