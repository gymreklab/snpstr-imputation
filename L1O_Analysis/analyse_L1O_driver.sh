#!/bin/bash

### Script to analyse the L1O analysis results
### Outputs a file *.csv

### Edit the python file before using this

ls *ground*txt | sed 's/.groundTruth.txt//' > samplesID.txt

python analyse_L1O.py .