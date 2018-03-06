#!/bin/bash

source params.sh

for locus in $LOCI
do
    nohup ./plot_r2_heatmap.py $locus &
done
