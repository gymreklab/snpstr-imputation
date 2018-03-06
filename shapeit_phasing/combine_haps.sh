#!/bin/sh

### script combines all ShapeIT hap files generated
### usage: 
### ./combine_haps.sh -s /path/to/shapeit -b /path/to/bcftools
###

while getopts s:b: option
do
        case "${option}"
        in
                s) shapeit=${OPTARG};;
                b) bcftools=$OPTARG;;
        esac
done

for c in `seq 1 22`;
do
(ls -v chr$c/shapeit.chr$c.reg*haps| xargs cat) > chr$c/shapeit.chr$c.haps
sort -k 3,3 -V chr$c/shapeit.chr$c.haps > shapeit.chr$c.haps
cp chr22/shapeit.chr22.reg0.sample shapeit.chr$c.sample

$shapeit -convert --input-haps shapeit.chr$c --output-vcf shapeit.chr$c.vcf

bgzip shapeit.chr$c.vcf
$bcftools index shapeit.chr$c.vcf.gz
done

ls -v shapeit*.vcf.gz > filelist.txt
$bcftools concat --file-list filelist.txt --output-type z --output shapeit.all.chr.vcf.gz