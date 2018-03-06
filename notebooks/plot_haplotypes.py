#!/usr/bin/env python
"""
Plot halotype structure for a target locus

Usage: ./plot_haplotypes.py <locus> <window>
"""

import matplotlib
matplotlib.use('Agg')
# Allow us to edit fonts in Illustrator
import matplotlib
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib import gridspec

import pandas as pd
import os
import numpy as np
import scipy.stats
import sys
import matplotlib.pyplot as plt

try:
    locus = sys.argv[1]
    WINDOW = sys.argv[2]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)
OUTPATH="pdfs/"

if locus == "SCA12":
    CHROM=5
    START=146258291
if locus == "VLDLR":
    CHROM=9
    START=2622145
if locus == "RUNX2a":
    CHROM=6
    START=45390419
if locus == "RUNX2b":
    CHROM=6
    START=45390419
if locus == "ATN1":
    CHROM=12
    START=7045880
if locus == "ATXN3":
    CHROM=14
    START=92537355
if locus == "JPH3":
    CHROM=16
    START=87637889
if locus == "DM1":
    CHROM=19
    START=46273457
if locus == "SCA6":
    CHROM=19
    START=13318669
if locus == "ATXN7":
    CHROM=3
    START=63898361
if locus == "PAPBN1":
    CHROM=17
    START=68735170
if locus == "ATXN2":
    CHROM=12
    START=112036754
if locus == "Gilbert":
    CHROM=2
    START=234668880
if locus == "HOXD13":
    CHROM=2
    START=176957786

####### Extract haplotype info #######
cmd="./extract_haplotypes.sh %s %s %s %s"%(CHROM, START, WINDOW, locus)
os.system(cmd)

######## Plot haplotypes #############
sys.stderr.write("Plot for haplotypes...\n")
# Read in haplotypes
numhaps = 1916*2
colnames = ["id","pos","ref","alt"] + ["hap_%s"%i for i in range(numhaps)]
haplotypes = pd.read_csv("haplotypes_%s.tab"%(locus), sep="\t",
                        names=colnames, usecols=range(len(colnames)))
haplotypes["vartype"] = haplotypes.apply(lambda x: ["SNP","STR"][int(len(x["ref"])>1)], 1)
haplotypes.index = ["pos"+str(haplotypes["pos"].values[i]) for i in range(haplotypes.shape[0]-1)] + ["STR"]

# Annotate STR lengths
ref = haplotypes[haplotypes["vartype"]=="STR"]["ref"].values[0]
alt = haplotypes[haplotypes["vartype"]=="STR"]["alt"].values[0].split(",")
str_allele_lengths = [len(ref)] + [len(item) for item in alt]
str_allele_lengths = [item-len(ref) for item in str_allele_lengths]
for i in range(numhaps):
    col = "hap_%s"%i
    gtlen = str_allele_lengths[haplotypes[haplotypes["vartype"]=="STR"][col].values[0]]
    haplotypes.loc["STR", col] = gtlen

# Reaad in allele-r2
hapcols = colnames[4:]
haplotype_filt = haplotypes[hapcols].transpose()

allsnps = [item for item in haplotype_filt.columns if "pos" in item]
ar2 = pd.read_csv("snp_loci_alleler2_%s.tab"%(locus), sep="\t")
ar2["pos"] = ar2["locus2"].apply(lambda x: "pos"+x.split(":")[1])
ar2 = ar2[ar2["pos"].apply(lambda x: x in allsnps)]
best_ar2 = ar2.groupby("pos", as_index=False).agg({"r2": max}).sort_values("r2", ascending=True)

# Get haplotype matrix, sort by allele-r2
haplotype_filt = haplotype_filt.sort_values(by="STR")
haplotype_filt = haplotype_filt.sort_values(["STR"]+list(best_ar2["pos"].values), ascending=False)

def PlotHapmaptrix(hapmatrix, ar2, allele, allsnps, fname):
    box_w =  1.0/len(allsnps)
    box_h = box_w
    hap_height = hapmatrix.shape[0]*0.0025
    legend_height = 0.5
    fig = plt.figure()
    fig.set_size_inches(3, hap_height + legend_height)
    gs = gridspec.GridSpec(2, 1, height_ratios=[hap_height, legend_height]) 
    ax = fig.add_subplot(gs[0])
    # Plot SNPs
    imx = ax.imshow(hapmatrix, cmap=plt.cm.Greys.from_list("snp", ["lightgray","black"]), 
              aspect="auto", extent=(0, hapmatrix.shape[1], box_h, hapmatrix.shape[0]-box_h))
    ax2 = fig.add_subplot(gs[1])
    # Plot snp allele r2
    cm = plt.cm.Blues.from_list("freq",["white","blue"])
    patches = []
    colors = []
    for i in range(len(allsnps)):
        r2 = ar2[(ar2["pos"] == allsnps[i]) & (ar2["allele"]==allele)]["r2"].values[0]
        x = i*box_w
        y = 0
        rect = mpatches.Rectangle([x, y], box_w, box_h)
        patches.append(rect)
        colors.append(cm(r2))
    collection = PatchCollection(patches, color=colors, edgecolor="black")
    ax2.add_collection(collection)

    ax.set_yticks([]);
    ax.set_yticklabels([]);
    ax.set_xticks([]);
    ax.set_xticklabels([]);
    ax2.set_ylim(bottom=0, top=box_h)
    ax2.set_yticks([]);
    ax2.set_yticklabels([]);
    ax2.set_xticks([]);
    ax2.set_xticklabels([]);
    ax.set_title("STR allele %s"%allele)
    fig.subplots_adjust(hspace=0)
    fig.savefig(fname)
    
for allele in sorted(list(set(str_allele_lengths))):
    hapmatrix = np.matrix(haplotype_filt[haplotype_filt["STR"]==allele][allsnps])
    if hapmatrix.shape[0]>= 10:
        sys.stderr.write("%s:%s\n"%(allele, hapmatrix.shape))
        fname = os.path.join(OUTPATH, "PathogenicHaplotypes_%s_%s.pdf"%(locus,allele))
        PlotHapmaptrix(hapmatrix, ar2, allele, allsnps, fname)


