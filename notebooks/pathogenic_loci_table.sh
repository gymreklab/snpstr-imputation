#!/bin/bash

source params.sh

for locus in $LOCI
do
    start=$(cat snp_loci_${locus}_all.bed | grep -v locus2 | datamash -f max 7 | cut -f 2 | cut -d':' -f 2)
    chrom=$(cat snp_loci_${locus}_all.bed | grep -v locus2 | datamash -f max 7 | cut -f 2 | cut -d':' -f 1)
    r2=$(cat snp_loci_${locus}_all.bed | grep -v locus2 | datamash -f max 7 | cut -f 7)
    rsid=$(tabix /storage/s1saini/hipstr_allfilters/str_snp/chr${chrom}.str.snp.feb18.vcf.gz ${chrom}:${start}-${start} | cut -f 3)
    echo $locus $rsid $r2

    # TODO LOO results
done
