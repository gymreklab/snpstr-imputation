#!/usr/bin/env python
"""
Make Vasarely plot for a given locus from an STR VCF

./vplot.py \
  --vcf /storage/mgymrek/ssc-imputation/filtered_vcfs/hipstr.chr21.allfilters.vcf.gz \
  --out test \
  --locus 21:15871808
"""

import argparse
import sys
import vcf

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import scipy.stats

def Vasarely(obs, exp, fname):
    # First make sure normalized
    obs_total = sum(obs.values())*1.0
    exp_total = sum(exp.values())*1.0
    for gt in obs.keys():
        obs[gt] = obs[gt]/obs_total
    for gt in exp.keys():
        exp[gt] = exp[gt]/exp_total
    
    # Get alleles
    box_w = 1; box_h = box_w
    alleles = set()
    for gt in list(obs.keys()) + list(exp.keys()):
        alleles.add(gt[0])
        alleles.add(gt[1])
    alleles = sorted(list(alleles))
    cm = plt.cm.Greys.from_list("freq",["white","black"])
    n_alleles = len(alleles)
    fig, ax = plt.subplots()

    patches = []
    colors = []
    
    # Get expected (square)
    for i in range(len(alleles)):
        for j in range(len(alleles)):
            x = i*box_w
            y = j*box_h
            rect = mpatches.Rectangle([x,y], box_w, box_h)
            patches.append(rect)
            colors.append(cm(exp.get((alleles[i], alleles[j]), 0)))
    # Get observed (circle)
    for i in range(len(alleles)):
        for j in range(len(alleles)):
            x = i*box_w
            y = j*box_h
            circ = mpatches.Circle([x+box_w/2.0,y+box_h/2.0], box_w*0.4)
            patches.append(circ)
            colors.append(cm(obs.get((alleles[i], alleles[j]), 0)))

    # Plot shapes
    collection = PatchCollection(patches, color=colors)
    ax.add_collection(collection)
    plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
    
    # Plot labels
    for i in range(len(alleles)):
        xaxis_coord = [i*box_w+box_w/2.0, -1*box_h/4.0]
        yaxis_coord = [-1*box_w/5.0, i*box_h+box_h/2.0]
        for coord in [xaxis_coord, yaxis_coord]:
            plt.text(coord[0], coord[1], alleles[i], ha="center", family='sans-serif', size=14)
    
    plt.axis('equal')
    plt.axis('off')
    plt.savefig(fname)

def Vplot(record, outfile, samples=[]):
    """
    Make Vasarely plot of obs vs. expected genotype frequencies
    Also print obs het, exp het, and binomial p-value
    """
    # Get genotypes, allele frequencies
    allele_counts = {}
    obs_gts = {}
    obs_het = 0
    obs_hom = 0
    total = 0
    for sample in record:
        if len(samples)>0 and sample.sample not in samples: continue
        if sample["GT"] == "." or sample["GT"] == "./.": continue
        gt = map(int, sample["GB"].split("|"))
        if gt[0] > gt[1]: gt = gt[::-1] # make smaller allele first
        obs_gts[tuple(gt)] = obs_gts.get(tuple(gt), 0) + 1
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
    # Get Expected genotype frequencies
    exp_gts = {}
    exp_hom_frac = 0
    for a1 in allele_freqs.keys():
        exp_hom_frac += allele_freqs[a1]**2
        for a2 in allele_freqs.keys():
            if a1 > a2: continue
            exp_gts[(a1, a2)] = 2*allele_freqs[a1]*allele_freqs[a2]
    # Binomial test for HWE
    hwe_p = scipy.stats.binom_test(obs_het, n=obs_het+obs_hom, p=1-exp_hom_frac)
    # Make Vasarely plot
    Vasarely(obs_gts, exp_gts, outfile)
    # Print stats to screen
    obs_hom_frac = obs_hom*1.0/(obs_het+obs_hom)
    sys.stdout.write("%s:%s\texp_hom_frac=%s\tobs_hom_frac=%s\tp=%s\n"%(record.CHROM, record.POS, exp_hom_frac, obs_hom_frac, hwe_p))

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcf", help="VCF with STR genotypes.", type=str, required=True)
    parser.add_argument("--out", help="Name of output figure", type=str, required=True)
    parser.add_argument("--locus", help="Locus to plot. chr:start", type=str, required=True)
    parser.add_argument("--unrelated-samples", help="Restrict to these samples to calculate popgen stats", type=str, required=False)
    args = parser.parse_args()

    if args.unrelated_samples:
        samples = [item.strip() for item in open(args.unrelated_samples, "r").readlines()]
    else: samples = []

    chrom, start = args.locus.split(":")
    start = int(start)
    str_reader = vcf.Reader(open(args.vcf, "rb"))
    str_records = str_reader.fetch(chrom, start, start+1)
    for r in str_records:
        if r.POS == start:
            Vplot(r, args.out, samples=samples)
            sys.exit(0)
    sys.stderr.write("ERROR: Could not find locus %s\n"%args.locus)
    sys.exit(1)

if __name__ == "__main__":
    main()
