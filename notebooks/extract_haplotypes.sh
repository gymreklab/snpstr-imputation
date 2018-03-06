#!/bin/bash

CHROM=$1
START=$2
WINDOW=$3
LOCUS=$4
VCF=/storage/s1saini/hipstr_allfilters/str_snp/chr${CHROM}.str.snp.feb18.vcf.gz

# Get SNPs to use (only ones in LD with the STR)
# Min maf 0.01, min r2=0.01
echo "Getting SNPs to use"
/home/mgymrek/workspace/ssc-imputation/snpstr-ld/snp_str_ld_calculator.py \
  --str-vcf /storage/s1saini/hipstr_rerun/chr${CHROM}/hipstr.chr${CHROM}.with.1kg.filtered.vcf.gz \
  --snp-vcf /storage/resources/datasets/SSC_SNP_v2/shapeit.chr${CHROM}.with.ref.vcf.gz \
  --pairwise-snpstr \
  --region ${CHROM}:$START-$START \
  --max-dist $WINDOW > snp_loci_${LOCUS}_all.bed
cat snp_loci_${LOCUS}_all.bed | grep -v "nan" | grep -v locus2 | \
    awk '($5>=0.01 && $7>0.01) {print $2}'| sed 's/:/\t/' | \
    awk '{print $1 "\t" $2-1 "\t" $2}' > snp_loci_${LOCUS}.bed

# Get allele-r2 to sort on later
echo "Getting allele-r2"
#/home/mgymrek/workspace/ssc-imputation/snpstr-ld/snp_str_ld_calculator.py \
#  --str-vcf /storage/s1saini/hipstr_rerun/chr${CHROM}/hipstr.chr${CHROM}.with.1kg.filtered.vcf.gz \
#  --snp-vcf /storage/resources/datasets/SSC_SNP_v2/shapeit.chr${CHROM}.with.ref.vcf.gz \
#  --pairwise-snpstr --allele-r2 \
#  --region ${CHROM}:$START-$START \
#  --max-dist $WINDOW | grep -v "nan" > snp_loci_alleler2_${LOCUS}.tab
  
# Extract haplotypes
echo "Extracting haplotypes"
bcftools query -R snp_loci_${LOCUS}.bed \
     -f "%ID\t%POS\t%REF\t%ALT\t[%GT\t]\n" \
     $VCF | awk '(length($3)==1)' | sed 's/|/\t/g' > haplotypes_${LOCUS}.tab
bcftools query -r ${CHROM}:${START}-${START} \
    -f "%ID\t%POS\t%REF\t%ALT\t[%GT\t]\n" ${VCF} | \
    awk -v"start=$START" '($2==start)' | \
    sed 's/|/\t/g' >> haplotypes_${LOCUS}.tab
