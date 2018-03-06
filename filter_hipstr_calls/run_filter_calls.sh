#!/bin/bash

source params.sh

sbatch \
    -A csd568 \
    --array=1-2 \
    -p shared \
    --mem=20G \
    -t 2800 \
    --get-user-env \
    --job-name=filter_hipstr_vcfs_take \
    -o ${LOGDIR}/filter_hipstr_vcfs_%a.out -e ${LOGDIR}/filter_hipstr_vcfs_%a.err \
    ./filter_calls.sh