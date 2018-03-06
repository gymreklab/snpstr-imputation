#!/bin/bash

### usage:
### ./L1O.sh -s sampleID -p pedigree.fam -l str.snp.data.vcf.gz -t hipstrcalls.groundtruth.vcf.gz -b path/to/beagle -q snps.vcf.gz
### Eg: ./L1O.sh -s SSC00092 -p pedigree.fam -l final.str.snp.vcf.gz -t str-vcf-paired-v2.vcf.gz -b beagle.27Jul16.86a.jar -q shapeit.snps.vcf.gz -d workingDir
###

while getopts s:p:f:t:b:l:q:d: option
do
        case "${option}"
        in
                s) sss=${OPTARG};;
                p) ped=${OPTARG};;
                t) str=$OPTARG;;
                b) beagle=$OPTARG;;
                l) procStr=$OPTARG;;
                q) snp=$OPTARG;;
		d) dir=$OPTARG;;
        esac
done

cd ${dir}

SAMPLEID=$sss
echo $SAMPLEID
famID=`cat $ped | grep $SAMPLEID | head -1 | awk '{print $1}'`
CHILD1=`cat $ped | awk -v var="$famID" '$1==var {print $0}' | awk '$3!=0 {print $0}' | head -1 | awk '{print $2}'`
CHILD2=`cat $ped | awk -v var="$famID" '$1==var {print $0}' | awk '$3!=0 {print $0}' | tail -n 1 | awk '{print $2}'`
FATHERID=`cat $ped | awk -v var="$famID" '$1==var {print $0}' | awk '$3!=0 {print $0}' | head -1 | awk '{print $3}'`
MOTHERID=`cat $ped | awk -v var="$famID" '$1==var {print $0}' | awk '$3!=0 {print $0}' | head -1 | awk '{print $4}'`

cat $ped | awk '$3==0 {print $2}' | grep -v $FATHERID | grep -v $MOTHERID | grep -v $CHILD1 | grep -v $CHILD2 > sampleRef.${SAMPLEID}.txt

bcftools query -f '%ID\n' $str > ID.${SAMPLEID}.txt
bcftools query -f '%CHROM\t%POS\n' $str > str.pos.${SAMPLEID}.txt
bcftools query -f '%CHROM\t%POS\n' $procStr > full.pos.${SAMPLEID}.txt
grep -Fxvf str.pos.${SAMPLEID}.txt full.pos.${SAMPLEID}.txt > snp.pos.${SAMPLEID}.txt

#this is incase related samples are present
#bcftools view $procStr --sample-file sampleRef.${SAMPLEID}.txt --output-type z --output-file ref.${SAMPLEID}.vcf.gz --force-samples
#we are dealing only with parents (all unrelated samples)
bcftools view $procStr --samples ^$SAMPLEID --no-update --output-type z --output-file ref.${SAMPLEID}.vcf.gz --force-samples
bcftools view $snp --samples $SAMPLEID --no-update --output-type z --output-file exclude.${SAMPLEID}.vcf.gz --force-samples
bcftools index -f ref.${SAMPLEID}.vcf.gz

java -Xmx4g -jar $beagle gt=exclude.${SAMPLEID}.vcf.gz ref=ref.${SAMPLEID}.vcf.gz out=imputed.${SAMPLEID}
bcftools index imputed.${SAMPLEID}.vcf.gz
bcftools view imputed.${SAMPLEID}.vcf.gz --include ID=@ID.${SAMPLEID}.txt -O z -o imputed.str.${SAMPLEID}.vcf.gz

python /oasis/scratch/comet/s1saini/temp_project/str-imputation/L1O_Analysis/write_bases.py imputed.str.${SAMPLEID}.vcf.gz $SAMPLEID $SAMPLEID.imputeResult.txt
python /oasis/scratch/comet/s1saini/temp_project/str-imputation/L1O_Analysis/write_bases.py $str $SAMPLEID $SAMPLEID.groundTruth.txt

sort -f ${SAMPLEID}.groundTruth.txt > ${SAMPLEID}.ground.sorted.txt
sort -f ${SAMPLEID}.imputeResult.txt > ${SAMPLEID}.imputed.sorted.txt
join -1 1 -2 1 ${SAMPLEID}.ground.sorted.txt ${SAMPLEID}.imputed.sorted.txt | awk 'NF==5{print}' > ${SAMPLEID}.diff.txt
rm ${SAMPLEID}.groundTruth.txt ${SAMPLEID}.imputeResult.txt ${SAMPLEID}.ground.sorted.txt ${SAMPLEID}.imputed.sorted.txt

echo "done"

rm exclude.${SAMPLEID}.vcf.gz ref.${SAMPLEID}.vcf.gz sampleExclude.${SAMPLEID}.txt sampleRef.${SAMPLEID}.txt ref.${SAMPLEID}.vcf.gz.csi imputed.${SAMPLEID}.vcf.gz imputed.${SAMPLEID}.vcf.gz.csi
rm ID.${SAMPLEID}.txt str.pos.${SAMPLEID}.txt full.pos.${SAMPLEID}.txt snp.pos.${SAMPLEID}.txt imputed.${SAMPLEID}.log