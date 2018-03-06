#!/bin/bash

source params.sh

sbatch \
    -A csd568 \
    --array=1-22 \
    -p shared \
    --mem=3G \
    -t 2000 \
    --get-user-env \
    --job-name=sampstats \
    -o ${LOGDIR}/sampstats_%a.out -e ${LOGDIR}/sampstats_%a.err \
    ./get_sampstats.sh