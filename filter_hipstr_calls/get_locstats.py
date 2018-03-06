#!/usr/bin/env python

import argparse
import numpy as np
import scipy.stats
import sys
import vcf

def GetLocusStats(record, samples=[]):
    hwe_p = 0
    het = 0
    # Get genotypes, allele frequencies
    allele_counts = {}
    obs_het = 0
    obs_hom = 0
    total = 0
    for sample in record:
        if len(samples)>0 and sample.sample not in samples: continue
        if sample["GT"] == "." or sample["GT"] == "./.": continue
        gt = map(int, sample["GB"].split("|"))
        if gt[0] == gt[1]: obs_hom += 1
        else:
            obs_het += 1
        total += 1
        for al in gt:
            allele_counts[al] = allele_counts.get(al, 0) + 1
    # Get Allele frequencies
    allele_freqs = {}
    for key in allele_counts.keys():
        allele_freqs[key] = allele_counts[key]*1.0/sum(allele_counts.values())
    # Get expected num homs/hets
    exp_hom_frac = 0
    for al in allele_freqs.keys():
        exp_hom_frac += allele_freqs[al]**2
    # Binomial test for HWE
    hwe_p = scipy.stats.binom_test(obs_het, n=obs_het+obs_hom, p=1-exp_hom_frac)
    # Compute heterozygosity
    het = 1-sum([allele_freqs[al]**2 for al in allele_freqs.keys()])
    # Get mean allele length
    mean_allele = sum([al*allele_freqs[al] for al in allele_freqs])
    return hwe_p, het, mean_allele

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcf", help="VCF file", type=str, required=True)
    parser.add_argument("--unrelated-samples", help="Restrict to these samples to calculate popgen stats", type=str, required=False)
    args = parser.parse_args()

    if args.unrelated_samples:
        samples = [item.strip() for item in open(args.unrelated_samples, "r").readlines()]
    else: samples = []

    reader = vcf.Reader(open(args.vcf, "rb"))
    for record in reader:
        hwe_p, het, mean_allele = GetLocusStats(record, samples=samples)
        num_calls = record.INFO["AN"]/2
        sys.stdout.write("\t".join(map(str, [record.CHROM, record.INFO["START"], record.INFO["END"], \
                                                 hwe_p, num_calls, het, mean_allele]))+"\n")
        sys.stdout.flush()

if __name__ == "__main__":
    main()


