#!/bin/bash

### usage: ./script.sh ssc.vcf.gz ref.samples

### input - 1kg ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
### input - SSC - ssc.chr21.vcf.gz
### change population labels in plot_pca.py

ssc=$1
ref=$2

wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/individual_chromosomes/chr21.1kg.phase3.v5a.vcf.gz
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/individual_chromosomes/chr21.1kg.phase3.v5a.vcf.gz.tbi

bcftools view chr21.1kg.phase3.v5a.vcf.gz --samples-file ${ref} --output-type z --output-file 1000g.ref.vcf.gz
bcftools index 1000g.ref.vcf.gz

bcftools merge -m all 1000g.ref.vcf.gz ${ssc} --output-type z --output chr21.combined.vcf.gz

plink --vcf chr21.combined.vcf.gz --geno 0.3 --make-bed --out chr21.combined
plink --bfile chr21.combined --pca 10 --out pca_10

python plot_pca.py
