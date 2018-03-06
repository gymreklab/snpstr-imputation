#!/bin/bash

### usage: 
### ./phase.sh -v vcffile.vcf.gz -f famfile.txt -s /path/to/shapeit -b /path/to/bcftools -p /path/to/plink -m directory_to_genetic_maps -q pythonFamScript
###

while getopts v:f:s:b:p:m:q: option
do
        case "${option}"
        in
                v) VCF=${OPTARG};;
                f) famFile=${OPTARG};;
                s) shapeit=${OPTARG};;
                b) bcftools=$OPTARG;;
                p) plink=$OPTARG;;
                m) maps=$OPTARG;;
                q) pythonScript=$OPTARG;;
        esac
done

for c in `seq 1 22`;
do
$bcftools query -f '%POS\n' -r $c $VCF > chr$c.POS.txt
sed -n '0~80000p' chr$c.POS.txt > chr$c.POS.sed.txt
(echo "1"; cat chr$c.POS.sed.txt; tail -1 chr$c.POS.txt) > chr$c.POS.v2.txt

readarray -t POS < chr$c.POS.v2.txt
echo ${POS[*]}

nlines=`wc -l < chr$c.POS.v2.txt`
let "nlines=$nlines-2"

for i in `seq 0 $nlines`;
do
echo $bcftools view $VCF -r $c:${POS[$i]}-${POS[$i+1]} --output-type z --output-file chr$c.reg$i.vcf.gz -v snps
$bcftools view $VCF -r $c:${POS[$i]}-${POS[$i+1]} --output-type z --output-file chr$c.reg$i.vcf.gz -v snps

echo $plink --noweb --vcf chr$c.reg$i.vcf.gz --make-bed --out chr$c.reg$i
$plink --noweb --vcf chr$c.reg$i.vcf.gz --keep-allele-order --make-bed --out chr$c.reg$i

python $pythonScript $famFile chr$c.reg$i.fam chr$c.reg$i.fam

echo $plink --bfile chr$c.reg$i --me 1 1 --set-me-missing --make-bed --out chr$c.reg$i\_me
$plink --bfile chr$c.reg$i --keep-allele-order --me 1 1 --set-me-missing --make-bed --out chr$c.reg$i\_me

echo $shapeit -B chr$c.reg$i\_me -M $maps/genetic_map_chr$c\_combined_b37.txt --duohmm -W 5 -O shapeit.chr$c.reg$i --thread 16
$shapeit -B chr$c.reg$i\_me -M $maps/genetic_map_chr$c\_combined_b37.txt --duohmm -W 5 -O shapeit.chr$c.reg$i --thread 16
done

mkdir chr$c
mv shapeit* chr$c/
#mv duohmm.chr$c.reg* chr$c/
mv chr$c.* chr$c/
done
