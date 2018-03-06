#!/bin/bash

### usage: ./runme.sh phasedstr.vcf.gz snp.vcf.gz snp.str.vcf.gz chrom

###
### ls imputed.str.SS* | sed 's/imputed.str.//' | sed 's/.vcf.gz//' | sdiff sample.txt - | grep "<" | awk '{print $1}' > sample_rerun.txt
### c=2
### ./rerunme.sh /oasis/scratch/comet/s1saini/temp_project/phased_feb18/hipstr.chr$c.phased.vcf.gz /oasis/scratch/comet/s1saini/temp_project/shapeit/shapeit.chr$c.with.ref.reorder.vcf.gz chr$c.snp.str.vcf.gz $c
### ls -v sample.iter*.sh | xargs -i sbatch {}
###

phasedstr=$1
snp=$2
snpstr=$3
chrom=$4

phasedstr=hipstr.chr${chrom}.phased.parents.vcf.gz
snp=shapeit.chr${chrom}.parents.vcf.gz

pwd=`pwd`

rm sample.iter*
split -l 20 sample_rerun.txt sample.iter
for i in `ls -v sample.iter*`;
do
fname=`echo $i.sh`
echo "#!/bin/bash" > $fname
echo "#SBATCH -A ddp268" >> $fname
echo "#SBATCH -p shared" >> $fname
echo "#SBATCH --get-user-env" >> $fname
echo "#SBATCH -t 1-23:59 # Walltime Days-Hours-Minutes" >> $fname
echo "#SBATCH --mem=5G" >> $fname
echo "#SBATCH -c 1" >> $fname

echo "cat $i | xargs -i ./L1O.sh -s {} -p $pwd/pedigree.fam -l $pwd/${snpstr} -t $pwd/${phasedstr} -b $pwd/beagle.08Jun17.d8b.jar -q $pwd/${snp} -d /scratch/\$USER/\$SLURM_JOBID" >> $fname

echo "wait" >> $fname
echo "cp /scratch/\$USER/\$SLURM_JOBID/*.diff.txt . " >> $fname
echo "cp /scratch/\$USER/\$SLURM_JOBID/imputed.str.*.vcf.gz . " >> $fname
done