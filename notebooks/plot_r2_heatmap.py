#!/usr/bin/env python
"""
usage: ./plot_r2_heatmap.py <locus>
"""

import sys
try:
    locus = sys.argv[1]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

OUTPATH="pdfs/"

import pandas as pd
import os

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib import gridspec

import vcf
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load list of loci to restrict to
usesnps = pd.read_csv("snp_loci_%s.bed"%locus, sep="\t", names=["chrom","start","snppos"])
# Load allele r2
ar2 = pd.read_csv("snp_loci_alleler2_%s.tab"%locus, sep="\t")
ar2["snppos"] = ar2["locus2"].apply(lambda x: int(x.split(":")[1]))
ar2 = pd.merge(ar2, usesnps, on=["snppos"])

chrom=ar2["locus1"].values[0].split(":")[0]
# Convert snp pos to rsid
reader = vcf.Reader(open("/storage/s1saini/hipstr_allfilters/str_snp/chr%s.str.snp.feb18.vcf.gz"%chrom, "rb"))
rsids = []
for i in range(ar2.shape[0]):
    start = int(ar2["locus2"].values[i].split(":")[1])
    records = reader.fetch(chrom, start-1, start)
    rid = None
    for r in records:
        rid = r.ID.replace("_","")
        break
    rsids.append(rid)
ar2["rsid"] = rsids

# Make matrix of allele num vs. SNP
alleles = sorted(list(set(ar2["allele"])))
snps = sorted(list(set(ar2["snppos"])))
rsids = [ar2[ar2["snppos"]==item]["rsid"].values[0] for item in snps]

data = np.zeros((len(alleles), len(snps)))
for i in range(len(alleles)):
    for j in range(len(snps)):
        r2 = ar2[(ar2["allele"]==alleles[i]) & (ar2["snppos"]==snps[j])]["r2"].values[0]
        data[i,j] = r2

fig = plt.figure()
fig.set_size_inches((25, 5))
ax = fig.add_subplot(111)
sns.heatmap(data, cmap=plt.cm.Blues, vmin=0, vmax=1, yticklabels=alleles, xticklabels=rsids, ax=ax, cbar=False);
ax.set_xticklabels(ax.get_xticklabels(), rotation=80, size=10);
ax.set_yticklabels(ax.get_yticklabels(), size=12, rotation=0);
ax.set_ylabel("bp relative to hg19", size=20);
fig.savefig(os.path.join(OUTPATH, "Pathogenic_%s_alleler2.pdf"%locus))
