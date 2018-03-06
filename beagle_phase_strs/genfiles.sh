#!/bin/bash

#SBATCH -A ddp268
#SBATCH -p shared
#SBATCH --get-user-env
#SBATCH -t 0-12:00 # Walltime Days-Hours-Minutes
#SBATCH --mem=32G
#SBATCH -c 4

##### usage: ./thisfile.sh hipstr_calls_21_minimal_unphased.vcf.gz shapeit.chr21.reorder.vcf.gz

hipstrfile=hipstr.chr22.with.1kg.vcf.gz
shapeitInput=shapeit.chr22.with.ref.vcf.gz
chr=22

shapeitfile=shapeit.chr$chr.with.ref.reorder.vcf.gz

bcftools index $hipstrfile
bcftools index $shapeitInput

bcftools query -f '%ID\n' $hipstrfile > ID.txt
bcftools query -l $hipstrfile > indv.txt
cp indv.txt sample.txt
bcftools view $shapeitInput --samples-file indv.txt -O z -o $shapeitfile
bcftools index $shapeitfile



pilot=`echo pilot.sh`
rm $pilot
echo "#!/bin/bash" >> $pilot

bcftools query -f '%POS\n' $hipstrfile > POS.txt
nlines=`wc -l < POS.txt`

for i in `seq 1 $nlines`;
do

POS=`sed "${i}q;d" POS.txt`
let "start=$POS-50000"
let "end=$POS+50000"

fname=`echo chr$chr.reg$i.sh`

echo "#!/bin/bash" > $fname
echo "#SBATCH -A ddp268" >> $fname
echo "#SBATCH -p shared" >> $fname
echo "#SBATCH --get-user-env" >> $fname
echo "#SBATCH -t 0-00:20 # Walltime Days-Hours-Minutes" >> $fname
echo "#SBATCH --mem=32G" >> $fname
echo "#SBATCH -c 4" >> $fname

echo "bcftools view $hipstrfile -r $chr:$POS -O z -o str$i.vcf.gz" >> $fname
echo "bcftools index -f str$i.vcf.gz" >> $fname
echo "vcftools --gzvcf str$i.vcf.gz --missing-indv --stdout | sort -k5 -r | awk '(NR>=2) && (\$5>0)' | awk '{print \$1}' > str$i.missing.samples" >> $fname
echo "vcftools --gzvcf str$i.vcf.gz --missing-indv --stdout | sort -k5 -r | awk '(NR>=2) && (\$5==0)' | awk '{print \$1}' > str$i.valid.samples" >> $fname

echo "bcftools view str$i.vcf.gz --samples-file str$i.missing.samples | bcftools annotate -x ^FORMAT/GT,^FORMAT/GB -O z -o str$i.missing.vcf.gz" >> $fname
echo "bcftools view str$i.vcf.gz --samples-file str$i.valid.samples | bcftools annotate -x ^FORMAT/GL -O z -o str$i.valid.vcf.gz" >> $fname
echo "bcftools index -f str$i.missing.vcf.gz" >> $fname
echo "bcftools index -f str$i.valid.vcf.gz" >> $fname

echo "bcftools view $shapeitfile -r $chr:$start-$end --samples-file str$i.missing.samples -O z -o snp$i.missing.vcf.gz" >> $fname
echo "bcftools view $shapeitfile -r $chr:$start-$end --samples-file str$i.valid.samples -O z -o snp$i.valid.vcf.gz" >> $fname
echo "bcftools index -f snp$i.missing.vcf.gz" >> $fname
echo "bcftools index -f snp$i.valid.vcf.gz" >> $fname

echo "bcftools concat -a snp$i.valid.vcf.gz str$i.valid.vcf.gz -O z -o str$i.snp.valid.vcf.gz" >> $fname
echo "vcf-sort -c str$i.snp.valid.vcf.gz > str$i.snp.valid.sort.vcf" >> $fname
echo "bgzip str$i.snp.valid.sort.vcf" >> $fname

echo "bcftools concat -a snp$i.missing.vcf.gz str$i.missing.vcf.gz -O z -o str$i.snp.missing.vcf.gz" >> $fname
echo "vcf-sort -c str$i.snp.missing.vcf.gz > str$i.snp.missing.sort.vcf" >> $fname
echo "bcftools view str$i.snp.missing.sort.vcf | sed 's/\.\:\./\.\/\.\:\./g' > str$i.snp.missing.sort.delim.vcf" >> $fname
echo "bgzip str$i.snp.missing.sort.delim.vcf" >> $fname
echo "bcftools index -f str$i.snp.missing.sort.delim.vcf.gz" >> $fname

echo "java -Xmx32g -jar beagle.r1399.jar gtgl=str$i.snp.valid.sort.vcf.gz out=beagle.str$i.valid nthreads=4 usephase=true ped=pedigree.fam" >> $fname
echo "bcftools index -f beagle.str$i.valid.vcf.gz" >> $fname

echo "missing=`cat str$i.missing.samples | wc -l`" >> $fname
echo "if [ $missing -eq 0 ]" >> $fname
echo "then" >> $fname
echo "cp beagle.str$i.valid.vcf.gz beagle.str$i.phased.vcf.gz" >> $fname
echo "bcftools index beagle.str$i.phased.vcf.gz" >> $fname
echo "else" >> $fname
echo "bcftools merge -m id beagle.str$i.valid.vcf.gz str$i.snp.missing.sort.delim.vcf.gz -O z -o str$i.missing.valid.combined.vcf.gz" >> $fname
echo "java -Xmx32g -jar beagle.r1399.jar gt=str$i.missing.valid.combined.vcf.gz out=beagle.str$i.phased nthreads=4 usephase=true ped=pedigree.fam" >> $fname
echo "bcftools index beagle.str$i.phased.vcf.gz" >> $fname
echo "fi" >> $fname

echo "bcftools view beagle.str$i.phased.vcf.gz --include ID=@ID.txt --samples-file sample.txt -O z -o beagle.str$i.only.vcf.gz" >> $fname
echo "bcftools index beagle.str$i.only.vcf.gz" >> $fname

echo "rm str$i.vcf.gz str$i.vcf.gz.csi str$i.missing.samples out.log str$i.valid.samples str$i.missing.vcf.gz str$i.valid.vcf.gz str$i.valid.vcf.gz.csi str$i.missing.vcf.gz.csi snp$i.missing.vcf.gz snp$i.valid.vcf.gz snp$i.missing.vcf.gz.csi snp$i.valid.vcf.gz.csi str$i.snp.valid.vcf.gz str$i.snp.valid.sort.vcf.gz str$i.snp.missing.vcf.gz str$i.snp.missing.sort.vcf str$i.snp.missing.sort.delim.vcf.gz str$i.snp.missing.sort.delim.vcf.gz.csi beagle.str$i.valid.vcf.gz beagle.str$i.valid.warnings beagle.str$i.valid.log beagle.str$i.valid.vcf.gz.csi str$i.missing.valid.combined.vcf.gz beagle.str$i.phased.vcf.gz beagle.str$i.phased.warnings beagle.str$i.phased.log beagle.str$i.phased.vcf.gz.csi" >> $fname
chmod 777 $fname

echo "sbatch $fname" >> $pilot
if ! ((i % 500)); then
    echo "sleep 5m" >> $pilot
fi
done
