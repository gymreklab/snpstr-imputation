#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import scipy.stats
import sys
from cyvcf2 import VCF, Writer

MIN_HRUN_PERIOD = 5

def GetFilters(x, args, nsamp):
    filters = []
    if x["hwe.p"] < args.min_hwep: filters.append("HWE")
    if x["numcalls"] < args.min_callrate*nsamp: filters.append("Callrate")
    if x["het"] < args.min_het: filters.append("Het")
    if x["period"] >= MIN_HRUN_PERIOD and x["hrun"] > x["period"] - args.max_hrun_offset:
        filters.append("Hrun")
    if args.filter_segdup and x["segdup"]>0: filters.append("Segdup")
    if len(filters) == 0:
        return "."
    else: return ";".join(filters)

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcf", help="VCF file", type=str, required=True)
    parser.add_argument("--statsfile", help="File with chrom, start, locus stats", type=str, required=True)
    parser.add_argument("--out", help="Prefix for output files", type=str, required=True)
    parser.add_argument("--min-hwep", help="Minimum HWE p-value", type=float, default=0)
    parser.add_argument("--min-callrate", help="Minimum call rate", type=float, default=0)
    parser.add_argument("--min-het", help="Minimum heterozygosity", type=float, default=0)
    parser.add_argument("--max-hrun-offset", help="For periods 5+, discard if the ref has " \
                            "homopolymer run > period+offset", type=int, default=100000)
    parser.add_argument("--filter-segdup", help="Filter loci overlapping a segdup", action="store_true")
    args = parser.parse_args()

    # Get VCF reader
    reader = VCF(args.vcf)

    # Load locus filters
    sys.stderr.write("Getting filters...\n")
    locstats = pd.read_csv(args.statsfile, sep="\t")
    locstats["FILTER"] = locstats.apply(lambda x: GetFilters(x, args, len(reader.samples)), 1)
    locstats.to_csv(args.out + ".tab", sep="\t", index=False)

    # Get filter dictionary
    sys.stderr.write("Getting filter dictionary...\n")
    filter_dict = dict(zip(list(locstats["start"]), list(locstats["FILTER"])))

    # Set filter field
    sys.stderr.write("Setting filter field in VCFs...\n")
    adict = {
        "HWE": "HWE less than %s"%args.min_hwep,
        "Callrate": "Callrate less than %s"%args.min_callrate,
        "Het": "Het less than %s"%args.min_het,
        "Hrun": "Hrun greater than %s"%args.max_hrun_offset,
        "Segdup": "Locus in a segmental duplication",
        "MissingInfo": "No stats provided for the locus",
        }
    for f in adict:
        reader.add_filter_to_header({"ID": f, "Description": adict[f]})
    writer = Writer("/dev/stdout", reader)
    for record in reader:
        filters = filter_dict.get(record.INFO["START"], "MissingInfo")
        if filters != ".":
            record.FILTER = filters.split(";")
        else: record.FILTER = "PASS"
        writer.write_record(record)
    writer.close()
    reader.close()

if __name__ == "__main__":
    main()


