#!/bin/bash

source params.sh

sbatch \
    -A csd568 \
    --array=1,2,17,22 \
    -p shared \
    --mem=8G \
    -t 2000 \
    --get-user-env \
    --job-name=locstats \
    -o ${LOGDIR}/locstats_%a.out -e ${LOGDIR}/locstats_%a.err \
    ./get_locstats.sh