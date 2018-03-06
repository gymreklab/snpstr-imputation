#!/bin/bash

#SBATCH -A ddp268
#SBATCH -p shared
#SBATCH --get-user-env
#SBATCH -t 0-12:00 # Walltime Days-Hours-Minutes
#SBATCH --mem=5G
#SBATCH -c 1

##### usage: ./thisfile.sh hipstr_calls_21_minimal_unphased.vcf.gz shapeit.chr21.reorder.vcf.gz

chr=11
shapeitInput=shapeit.chr$chr.with.ref.vcf.gz
hipstrInput=/oasis/projects/nsf/csd568/mgymrek/ssc-quads/hipstr_vcfs/final/hipstr.chr${chr}.allfilters.vcf.gz

hipstrfile=hipstr.chr$chr.allfilters.pass.vcf.gz
bcftools view ${hipstrInput} -f PASS -O z -o ${hipstrfile}
bcftools index ${hipstrfile}

shapeitfile=/oasis/scratch/comet/s1saini/temp_project/shapeit/shapeit.chr$chr.with.ref.reorder.vcf.gz

bcftools query -f '%ID\n' $hipstrfile > ID.txt
bcftools query -l $hipstrfile > indv.txt
cp indv.txt sample.txt


pilot=`echo pilot.sh`
rm $pilot
echo "#!/bin/bash" >> $pilot

bcftools query -f '%POS\n' $hipstrfile > POS.txt
nlines=`wc -l < POS.txt`

for j in `seq 1 100 $nlines`;
do

fname=`echo chr$chr.reg$j.sh`
echo "#!/bin/bash" > $fname
echo "#SBATCH -A ddp268" >> $fname
echo "#SBATCH -p shared" >> $fname
echo "#SBATCH --get-user-env" >> $fname
echo "#SBATCH -t 1-23:59 # Walltime Days-Hours-Minutes" >> $fname
echo "#SBATCH --mem=5G" >> $fname
echo "#SBATCH -c 1" >> $fname

echo "DIR=/scratch/\$USER/\$SLURM_JOBID" >> $fname

let innerstart=$j
let innerend=$j+100-1

for i in `seq $innerstart $innerend`;
do
POS=`sed "${i}q;d" POS.txt`
let "start=$POS-50000"
let "end=$POS+50000"

echo "date" >> $fname
echo "bcftools view $hipstrfile -r $chr:$POS -O z -o \${DIR}/str$i.vcf.gz" >> $fname
echo "bcftools index -f \${DIR}/str$i.vcf.gz" >> $fname

echo "bcftools view $shapeitfile -r $chr:$start-$end -O z -o \${DIR}/snp$i.vcf.gz" >> $fname
echo "tabix -p vcf \${DIR}/snp$i.vcf.gz" >> $fname

echo "bcftools concat -a \${DIR}/snp$i.vcf.gz \${DIR}/str$i.vcf.gz -O z -o \${DIR}/str$i.snp.vcf.gz" >> $fname
echo "/oasis/scratch/comet/s1saini/temp_project/bcftools-1.6/bcftools sort \${DIR}/str$i.snp.vcf.gz -O z -o \${DIR}/str$i.snp.sort.vcf.gz" >> $fname
echo "bcftools index -f \${DIR}/str$i.snp.sort.vcf.gz" >> $fname

echo "date" >> $fname

echo "vcftools --gzvcf \${DIR}/str$i.vcf.gz --missing-indv --stdout | sort -k5 -r | awk '(NR>=2) && (\$5>0)' | awk '{print \$1}' > \${DIR}/str$i.missing.samples" >> $fname
echo "vcftools --gzvcf \${DIR}/str$i.vcf.gz --missing-indv --stdout | sort -k5 -r | awk '(NR>=2) && (\$5==0)' | awk '{print \$1}' > \${DIR}/str$i.valid.samples" >> $fname

echo "bcftools view \${DIR}/str$i.snp.sort.vcf.gz --samples-file \${DIR}/str$i.valid.samples | bcftools annotate -x ^FORMAT/GT,^FORMAT/GL -O z -o \${DIR}/str$i.snp.valid.sort.vcf.gz" >> $fname

echo "bcftools view \${DIR}/str$i.snp.sort.vcf.gz --samples-file \${DIR}/str$i.missing.samples | bcftools annotate -x ^FORMAT/GT,^FORMAT/GB | sed 's/\.\:\./\.\/\.\:\./g' | sed 's/\.\/\.\/\./\.\/\./g' > \${DIR}/str$i.snp.missing.sort.delim.vcf" >> $fname
echo "bgzip \${DIR}/str$i.snp.missing.sort.delim.vcf" >> $fname
echo "bcftools index \${DIR}/str$i.snp.missing.sort.delim.vcf.gz" >> $fname

echo "date" >> $fname
echo "timeout 5m java -Xmx5g -jar beagle.r1399.jar gtgl=\${DIR}/str$i.snp.valid.sort.vcf.gz out=\${DIR}/beagle.str$i.valid nthreads=4 usephase=true ped=pedigree.fam" >> $fname
echo "bcftools index -f \${DIR}/beagle.str$i.valid.vcf.gz" >> $fname

echo "missing=\`cat \${DIR}/str$i.missing.samples | wc -l\`" >> $fname
echo "if [ \$missing -eq 0 ]" >> $fname
echo "then" >> $fname
echo "cp \${DIR}/beagle.str$i.valid.vcf.gz \${DIR}/beagle.str$i.phased.vcf.gz" >> $fname
echo "bcftools index \${DIR}/beagle.str$i.phased.vcf.gz" >> $fname
echo "else" >> $fname
echo "/oasis/scratch/comet/s1saini/temp_project/bcftools-1.6/bcftools merge -m id \${DIR}/beagle.str$i.valid.vcf.gz \${DIR}/str$i.snp.missing.sort.delim.vcf.gz -O z -o \${DIR}/str$i.missing.valid.combined.vcf.gz" >> $fname
echo "timeout 5m java -Xmx5g -jar beagle.r1399.jar gt=\${DIR}/str$i.missing.valid.combined.vcf.gz out=\${DIR}/beagle.str$i.phased nthreads=4 usephase=true ped=pedigree.fam" >> $fname
echo "tabix -p vcf \${DIR}/beagle.str$i.phased.vcf.gz" >> $fname
echo "fi" >> $fname
echo "date" >> $fname

echo "./fix_switch_error_gymrek_old.py --phased-vcf \${DIR}/beagle.str$i.phased.vcf.gz --ref-vcf \${DIR}/snp$i.vcf.gz --switch-threshold 0.5 --min-maf 0.1 --check-snps 100 --new-vcf \${DIR}/beagle.str$i.phased.fixed.vcf" >> $fname
echo "bgzip \${DIR}/beagle.str$i.phased.fixed.vcf" >> $fname
echo "bcftools index \${DIR}/beagle.str$i.phased.fixed.vcf.gz" >> $fname

echo "date" >> $fname
echo "bcftools view \${DIR}/beagle.str$i.phased.fixed.vcf.gz --include ID=@ID.txt --samples-file sample.txt -O z -o \${DIR}/beagle.str$i.only.vcf.gz" >> $fname
echo "bcftools index \${DIR}/beagle.str$i.only.vcf.gz" >> $fname

echo "rm \${DIR}/beagle.str$i.phased.vcf.gz.tbi \${DIR}/str$i.snp.vcf.gz \${DIR}/str$i.snp.sort.vcf.gz \${DIR}/str$i.snp.sort.vcf.gz.csi \${DIR}/beagle.str$i.phased.fixed.vcf.gz \${DIR}/beagle.str$i.phased.fixed.vcf.gz.csi \${DIR}/ \${DIR}/beagle.str$i.phased.reorder.vcf.gz \${DIR}/beagle.str$i.phased.reorder.vcf.gz.tbi \${DIR}/snp$i.vcf.gz \${DIR}/snp$i.vcf.gz.tbi \${DIR}/str$i.vcf.gz \${DIR}/str$i.vcf.gz.csi \${DIR}/str$i.missing.samples \${DIR}/out.log \${DIR}/str$i.valid.samples \${DIR}/str$i.missing.vcf.gz \${DIR}/str$i.valid.vcf.gz \${DIR}/str$i.valid.vcf.gz.csi \${DIR}/str$i.missing.vcf.gz.csi \${DIR}/snp$i.missing.vcf.gz \${DIR}/snp$i.valid.vcf.gz \${DIR}/snp$i.missing.vcf.gz.csi \${DIR}/snp$i.valid.vcf.gz.csi \${DIR}/str$i.snp.valid.vcf.gz \${DIR}/str$i.snp.valid.sort.vcf.gz \${DIR}/str$i.snp.missing.vcf.gz \${DIR}/str$i.snp.missing.sort.vcf \${DIR}/str$i.snp.missing.sort.delim.vcf.gz \${DIR}/str$i.snp.missing.sort.delim.vcf.gz.csi \${DIR}/beagle.str$i.valid.vcf.gz \${DIR}/beagle.str$i.valid.warnings \${DIR}/beagle.str$i.valid.log \${DIR}/beagle.str$i.valid.vcf.gz.csi \${DIR}/str$i.missing.valid.combined.vcf.gz \${DIR}/beagle.str$i.phased.vcf.gz \${DIR}/beagle.str$i.phased.warnings \${DIR}/beagle.str$i.phased.log \${DIR}/beagle.str$i.phased.vcf.gz.csi \${DIR}/snp$i.valid.vcf.gz.csi \${DIR}/str$i.snp.valid.vcf.gz \${DIR}/str$i.snp.valid.sort.vcf.gz \${DIR}/str$i.snp.missing.vcf.gz \${DIR}/str$i.snp.missing.sort.vcf \${DIR}/str$i.snp.missing.sort.delim.vcf.gz \${DIR}/str$i.snp.missing.sort.delim.vcf.gz.csi \${DIR}/beagle.str$i.valid.warnings \${DIR}/beagle.str$i.valid.vcf.gz.csi \${DIR}/beagle.str$i.valid.vcf.gz \${DIR}/beagle.str$i.valid.log \${DIR}/str$i.missing.valid.combined.vcf.gz \${DIR}/beagle.str$i.phased.warnings \${DIR}/beagle.str$i.phased.vcf.gz \${DIR}/beagle.str$i.phased.log" >> $fname
echo "date" >> $fname
done

echo "sbatch $fname" >> $pilot
if ! ((j % 500)); then
    echo "sleep 5m" >> $pilot
fi

echo "cd \${DIR}" >> $fname
echo "ls -v beagle*only.vcf.gz.csi | sed 's/.csi//' > mylist.txt" >> $fname
echo "bcftools concat -f mylist.txt -a -d all -O z -o \${SLURM_JOBID}.vcf.gz" >> $fname
echo "bcftools index \${SLURM_JOBID}.vcf.gz" >> $fname
echo "cp \${SLURM_JOBID}.vcf.gz /oasis/scratch/comet/s1saini/temp_project/chr${chr}" >> $fname
echo "cp \${SLURM_JOBID}.vcf.gz.csi /oasis/scratch/comet/s1saini/temp_project/chr${chr}" >> $fname

done


