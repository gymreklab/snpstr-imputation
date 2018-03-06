#!/bin/bash

#SBATCH -A ddp268
#SBATCH -p shared
#SBATCH --get-user-env
#SBATCH -t 0-12:00 # Walltime Days-Hours-Minutes
#SBATCH --mem=5G
#SBATCH -c 1

##### usage: ./thisfile.sh hipstr_calls_21_minimal_unphased.vcf.gz shapeit.chr21.reorder.vcf.gz

chr=19
shapeitInput=shapeit.chr$chr.with.ref.vcf.gz
hipstrInput=/oasis/projects/nsf/csd568/mgymrek/ssc-quads/hipstr_vcfs/final/hipstr.chr${chr}.allfilters.vcf.gz

hipstrfile=hipstr.chr$chr.allfilters.pass.vcf.gz
bcftools view ${hipstrInput} -f PASS -O z -o ${hipstrfile}
bcftools index ${hipstrfile}

shapeitfile=/oasis/scratch/comet/s1saini/temp_project/shapeit/shapeit.chr$chr.with.ref.reorder.vcf.gz

#tabix -p vcf $hipstrfile
#tabix -p vcf $shapeitInput

bcftools query -f '%ID\n' $hipstrfile > ID.txt
bcftools query -l $hipstrfile > indv.txt
cp indv.txt sample.txt
#bcftools view $shapeitInput --samples-file indv.txt -O z -o $shapeitfile
#tabix -p vcf $shapeitfile



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

let innerstart=$j
let innerend=$j+100-1

for i in `seq $innerstart $innerend`;
do
POS=`sed "${i}q;d" POS.txt`
let "start=$POS-50000"
let "end=$POS+50000"

echo "date" >> $fname
echo "bcftools view $hipstrfile -r $chr:$POS -O z -o str$i.vcf.gz" >> $fname
echo "bcftools index -f str$i.vcf.gz" >> $fname

echo "bcftools view $shapeitfile -r $chr:$start-$end -O z -o snp$i.vcf.gz" >> $fname
echo "tabix -p vcf snp$i.vcf.gz" >> $fname

echo "bcftools concat -a snp$i.vcf.gz str$i.vcf.gz -O z -o str$i.snp.vcf.gz" >> $fname
echo "/oasis/scratch/comet/s1saini/temp_project/bcftools-1.6/bcftools sort str$i.snp.vcf.gz -O z -o str$i.snp.sort.vcf.gz" >> $fname
echo "bcftools index -f str$i.snp.sort.vcf.gz" >> $fname

echo "date" >> $fname

echo "vcftools --gzvcf str$i.vcf.gz --missing-indv --stdout | sort -k5 -r | awk '(NR>=2) && (\$5>0)' | awk '{print \$1}' > str$i.missing.samples" >> $fname
echo "vcftools --gzvcf str$i.vcf.gz --missing-indv --stdout | sort -k5 -r | awk '(NR>=2) && (\$5==0)' | awk '{print \$1}' > str$i.valid.samples" >> $fname

echo "bcftools view str$i.snp.sort.vcf.gz --samples-file str$i.valid.samples | bcftools annotate -x ^FORMAT/GT,^FORMAT/GL -O z -o str$i.snp.valid.sort.vcf.gz" >> $fname

echo "bcftools view str$i.snp.sort.vcf.gz --samples-file str$i.missing.samples | bcftools annotate -x ^FORMAT/GT,^FORMAT/GB | sed 's/\.\:\./\.\/\.\:\./g' | sed 's/\.\/\.\/\./\.\/\./g' > str$i.snp.missing.sort.delim.vcf" >> $fname
echo "bgzip str$i.snp.missing.sort.delim.vcf" >> $fname
echo "bcftools index str$i.snp.missing.sort.delim.vcf.gz" >> $fname

echo "date" >> $fname
echo "timeout 5m java -Xmx5g -jar beagle.r1399.jar gtgl=str$i.snp.valid.sort.vcf.gz out=beagle.str$i.valid nthreads=4 usephase=true ped=pedigree.fam" >> $fname
echo "bcftools index -f beagle.str$i.valid.vcf.gz" >> $fname

echo "missing=\`cat str$i.missing.samples | wc -l\`" >> $fname
echo "if [ \$missing -eq 0 ]" >> $fname
echo "then" >> $fname
echo "cp beagle.str$i.valid.vcf.gz beagle.str$i.phased.vcf.gz" >> $fname
echo "bcftools index beagle.str$i.phased.vcf.gz" >> $fname
echo "else" >> $fname
echo "/oasis/scratch/comet/s1saini/temp_project/bcftools-1.6/bcftools merge -m id beagle.str$i.valid.vcf.gz str$i.snp.missing.sort.delim.vcf.gz -O z -o str$i.missing.valid.combined.vcf.gz" >> $fname
echo "timeout 5m java -Xmx5g -jar beagle.r1399.jar gt=str$i.missing.valid.combined.vcf.gz out=beagle.str$i.phased nthreads=4 usephase=true ped=pedigree.fam" >> $fname
echo "tabix -p vcf beagle.str$i.phased.vcf.gz" >> $fname
echo "fi" >> $fname
echo "date" >> $fname

echo "./fix_switch_error_gymrek_old.py --phased-vcf beagle.str$i.phased.vcf.gz --ref-vcf snp$i.vcf.gz --switch-threshold 0.5 --min-maf 0.1 --check-snps 100 --new-vcf beagle.str$i.phased.fixed.vcf" >> $fname
echo "bgzip beagle.str$i.phased.fixed.vcf" >> $fname
echo "bcftools index beagle.str$i.phased.fixed.vcf.gz" >> $fname

echo "date" >> $fname
echo "bcftools view beagle.str$i.phased.fixed.vcf.gz --include ID=@ID.txt --samples-file sample.txt -O z -o beagle.str$i.only.vcf.gz" >> $fname
echo "bcftools index beagle.str$i.only.vcf.gz" >> $fname

echo "rm beagle.str$i.phased.vcf.gz.tbi str$i.snp.vcf.gz str$i.snp.sort.vcf.gz str$i.snp.sort.vcf.gz.csi beagle.str$i.phased.fixed.vcf.gz beagle.str$i.phased.fixed.vcf.gz.csi  beagle.str$i.phased.reorder.vcf.gz beagle.str$i.phased.reorder.vcf.gz.tbi snp$i.vcf.gz snp$i.vcf.gz.tbi str$i.vcf.gz str$i.vcf.gz.csi str$i.missing.samples out.log str$i.valid.samples str$i.missing.vcf.gz str$i.valid.vcf.gz str$i.valid.vcf.gz.csi str$i.missing.vcf.gz.csi snp$i.missing.vcf.gz snp$i.valid.vcf.gz snp$i.missing.vcf.gz.csi snp$i.valid.vcf.gz.csi str$i.snp.valid.vcf.gz str$i.snp.valid.sort.vcf.gz str$i.snp.missing.vcf.gz str$i.snp.missing.sort.vcf str$i.snp.missing.sort.delim.vcf.gz str$i.snp.missing.sort.delim.vcf.gz.csi beagle.str$i.valid.vcf.gz beagle.str$i.valid.warnings beagle.str$i.valid.log beagle.str$i.valid.vcf.gz.csi str$i.missing.valid.combined.vcf.gz beagle.str$i.phased.vcf.gz beagle.str$i.phased.warnings beagle.str$i.phased.log beagle.str$i.phased.vcf.gz.csi snp$i.valid.vcf.gz.csi str$i.snp.valid.vcf.gz str$i.snp.valid.sort.vcf.gz str$i.snp.missing.vcf.gz str$i.snp.missing.sort.vcf str$i.snp.missing.sort.delim.vcf.gz str$i.snp.missing.sort.delim.vcf.gz.csi beagle.str$i.valid.warnings beagle.str$i.valid.vcf.gz.csi beagle.str$i.valid.vcf.gz beagle.str$i.valid.log str$i.missing.valid.combined.vcf.gz beagle.str$i.phased.warnings beagle.str$i.phased.vcf.gz beagle.str$i.phased.log" >> $fname
echo "date" >> $fname
done

echo "sbatch $fname" >> $pilot
if ! ((j % 500)); then
    echo "sleep 5m" >> $pilot
fi
done

