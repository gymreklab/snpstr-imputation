#!/bin/bash

source params.sh

sbatch \
    -A csd568 \
    --array=1,2,17,22 \
    -p shared \
    --mem=3G \
    -t 2000 \
    --get-user-env \
    --job-name=locus_filters \
    -o ${LOGDIR}/locus_filters_%a.out -e ${LOGDIR}/locus_filters_%a.err \
    ./set_locus_filters.sh