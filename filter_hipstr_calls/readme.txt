# Run HipSTR call level filters
./run_filter_calls.sh

# Get locus stats
./run_locstats.sh

# Combine locus stats and set VCF filters
./run_set_locus_filters.sh

# Get sample stats
./run_sampstats.sh

# Get Mendelian inheritance for chr21
sbatch ./vcf_mend.sh 21

# Visualization
./vplot.py \
   --vcf /storage/mgymrek/ssc-imputation/filtered_vcfs/hipstr.chr21.allfilters.vcf.gz \
   --out bad_locus.pdf \
   --locus 21:9573172

