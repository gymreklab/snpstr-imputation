#!/bin/bash

### usage:
### ./combine.str.snp.sh -s str.vcf.gz -p snp.vcf.gz -o out.vcf.gz
###

while getopts s:p:o: option
do
        case "${option}"
        in
                s) str=${OPTARG};;
                p) snp=${OPTARG};;
                o) out=${OPTARG};;
        esac
done

bcftools index $str
bcftools view $str -m2 --threads 4 -O z -o str.bial.vcf.gz
bcftools index str.bial.vcf.gz

bcftools concat -a str.bial.vcf.gz $snp -O z -o $out
bcftools index $out

rm str.bial.vcf.gz str.bial.vcf.gz.csi
