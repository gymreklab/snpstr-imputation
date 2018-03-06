#!/bin/bash

### usage: ./runme.sh phasedstr.vcf.gz snp.vcf.gz snp.str.vcf.gz chrom
### usage: sbatch -p shared -A ddp268 -t 0-12:00 --job-name=chr21loo --mem=5G -c 1 --get-user-env ./runme.sh /oasis/scratch/comet/s1saini/temp_project/phased_feb18/hipstr.chr21.phased.vcf.gz /oasis/scratch/comet/s1saini/temp_project/shapeit/shapeit.chr21.with.ref.reorder.vcf.gz chr21.snp.str.vcf.gz 21


### then ls -v sample.iter*.sh | xargs -i sbatch {}
###

phasedstr=$1
snp=$2
snpstr=$3
chrom=$4

rm ID.SS*.txt
rm imputed.str.SS*.vcf.gz*
rm SS*.diff.txt

bcftools view ${phasedstr} --samples-file ssc_parents_id.txt -O z -o hipstr.chr${chrom}.phased.parents.vcf.gz --force-samples
bcftools index hipstr.chr${chrom}.phased.parents.vcf.gz
bcftools view ${snp} --samples-file ssc_parents_id.txt -O z -o shapeit.chr${chrom}.parents.vcf.gz --force-samples
bcftools index shapeit.chr${chrom}.parents.vcf.gz

phasedstr=hipstr.chr${chrom}.phased.parents.vcf.gz
snp=shapeit.chr${chrom}.parents.vcf.gz

./concat.sh -s ${phasedstr} -p ${snp} -o ${snpstr}

bcftools query -f '%ID\t%CHROM:%POS\n' ${phasedstr} | sort -f > POS_ID.txt

bcftools query -l ${snpstr} > sample.txt

pwd=`pwd`

rm sample.iter*
split -l 20 sample.txt sample.iter
for i in `ls -v sample.iter*`;
do
fname=`echo $i.sh`
echo "#!/bin/bash" > $fname
echo "#SBATCH -A ddp268" >> $fname
echo "#SBATCH -p shared" >> $fname
echo "#SBATCH --get-user-env" >> $fname
echo "#SBATCH -t 0-23:59 # Walltime Days-Hours-Minutes" >> $fname
echo "#SBATCH --mem=5G" >> $fname
echo "#SBATCH -c 1" >> $fname

echo "cat $i | xargs -i ./L1O.sh -s {} -p $pwd/pedigree.fam -l $pwd/${snpstr} -t $pwd/${phasedstr} -b $pwd/beagle.08Jun17.d8b.jar -q $pwd/${snp} -d /scratch/\$USER/\$SLURM_JOBID" >> $fname

echo "wait" >> $fname
echo "cp /scratch/\$USER/\$SLURM_JOBID/*.diff.txt . " >> $fname
echo "cp /scratch/\$USER/\$SLURM_JOBID/imputed.str.*.vcf.gz . " >> $fname
done