
Scripts and analyses for SNP-STR imputation manuscript

## beagle_phase_strs
Scripts to phase STR loci on XSEDE cluster computing environment.
In the the `genfiles_100.sh` file, edit the chromosome and input VCF files. Execute the file.
Then run the generated `pilot.sh` file to phase 100 STRs per compute job.

## filter_hipstr_calls
Scripts to filter HipSTR calls on the XSEDE cluster computing environment.

## hipstr_str_calling
Scripts to call STRs using HipSTR on Amazon AWS.

    ./launch_aws.sh "aws_access_key" "aws_secret_key" "chrosome" "[parts]" "volume size" "batch_size" "keyname"

## L1O_Analysis
Scripts to perform leave-one-out trials on XSEDE cluster computing environment.
Usage:

    sbatch -p shared \
    -A ddp268 \
    -t 0-12:00 \
    --job-name=chr21loo \
    --mem=5G \
    -c 1 \
    --get-user-env \
    ./runme.sh hipstr.calls.chr21.vcf.gz snps.chr21.vcf.gz chr21.snp.str.vcf.gz 21


## shapeit_phasing
Scripts to phase SNP haplotypes using SHAPEIT

    ./phase.sh -v unphased.snps.vcf.gz \
    -f pedigree.txt \
    -s /path/to/shapeit \
    -b /path/to/bcftools \
    -p /path/to/plink \
    -m directory_to_genetic_maps \
    -q createFam.py

## ssc-pca
Scripts to perform PCA using 1000 genomes reference panel
Input: snps.vcf.gz, 1000genomes.samples.txt
Usage: 

    ./script.sh -s ssc.vcf.gz -r ref.samples -c 21

