#!/bin/bash

#SBATCH -A ddp268
#SBATCH -p shared
#SBATCH --get-user-env
#SBATCH -t 0-12:00 # Walltime Days-Hours-Minutes
#SBATCH --mem=32G
#SBATCH -c 1

CHROM=22

rm files*

ls -v beagle*only*vcf.gz.csi | sed 's/.csi//' > files.list
split -l 100 files.list files
rm files.list

for i in `ls -v files*`;
do
bcftools concat -f $i -a -d all -O z -o $i.vcf.gz
bcftools index $i.vcf.gz
done

ls -v files*vcf.gz > files.list
bcftools concat -f files.list -a -d all -O z -o hipstr.chr${CHROM}.phased.vcf.gz
bcftools index hipstr.chr${CHROM}.phased.vcf.gz
