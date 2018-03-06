#!/bin/bash

# Plot haplotype structure
./run_hap_figs.sh

# SNP r2 heatmaps
./run_r2_heatmaps.sh

# Bubble plots
./rub_bubbles.sh

# Pull out info for table
./pathogenic_loci_table.sh

