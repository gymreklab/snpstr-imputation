#!/bin/bash

### usage: ./script.sh -s ssc.vcf.gz -r ref.samples -c 21


while getopts s:r:c: option
do
        case "${option}"
        in
                s) ssc=${OPTARG};;
                r) ref=${OPTARG};;
                c) chrom=${OPTARG};;
        esac
done

wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/individual_chromosomes/chr${chrom}.1kg.phase3.v5a.vcf.gz
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/individual_chromosomes/chr${chrom}.1kg.phase3.v5a.vcf.gz.tbi
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/integrated_call_samples_v3.20130502.ALL.panel

bcftools view chr${chrom}.1kg.phase3.v5a.vcf.gz --samples-file ${ref} --output-type z --output-file 1000g.ref.vcf.gz
bcftools index 1000g.ref.vcf.gz

bcftools merge -m all 1000g.ref.vcf.gz ${ssc} --output-type z --output chr${chrom}.combined.vcf.gz

plink --vcf chr${chrom}.combined.vcf.gz --geno 0.3 --make-bed --out chr${chrom}.combined
plink --bfile chr${chrom}.combined --pca 10 --out pca_10

num_ssc_samples=`bcftools query -l ${ssc} | wc -l`

python plot_pca.py