#!/bin/sh

ls *.bam > files.list
#cat files.list | xargs -I% -n1 samtools sort % -o %
cat files.list | xargs -I% -n1 samtools index %
