#!/usr/bin/env python

"""
Calculate LD between an STR and SNP variant

Note on outputs:
- freq_het: gives allele freq if using allele-r2, else heterozygosity
- maf: gives MAF of the SNP, or if 2nd VCF is STRs same as freq_het

# Test
./snp_str_ld_calculator.py \
  --str-vcf /storage/s1saini/hipstr_rerun/chr5/hipstr.chr5.with.1kg.filtered.vcf.gz \
  --snp-vcf /storage/resources/datasets/SSC_SNP_v2/shapeit.chr5.with.ref.vcf.gz \
  --str-locus 5:153681115 \
  --snp-locus-rsid rs11740474 \
  --use-info-start

./snp_str_ld_calculator.py \
  --str-vcf /storage/s1saini/hipstr_rerun/chr5/hipstr.chr5.with.1kg.filtered.vcf.gz \
  --snp-vcf /storage/resources/datasets/SSC_SNP_v2/shapeit.chr5.with.ref.vcf.gz \
  --pairwise-snpstr

./snp_str_ld_calculator.py \
  --str-vcf hipstr.1kg.EUR.filtered.vcf.gz \
  --str-vcf2 1kg.EUR.wgs.imputed.vcf.gz \
  --allele-r2 --mincount 3

./snp_str_ld_calculator.py \
  --str-vcf /storage/s1saini/hipstr_rerun/chr5/hipstr.chr5.with.1kg.filtered.vcf.gz \
  --snp-vcf /storage/resources/datasets/SSC_SNP_v2/shapeit.chr5.with.ref.vcf.gz \
  --loci-file-rsid /storage/mgymrek/ssc-imputation/pgc/bychrom/ldfile_5.tab \
  --samples ~/workspace/ssc-imputation/metadata/ssc_parent_ids.txt \
  --use-info-start

"""

import argparse
import numpy as np
import sys
import vcf
import scipy.stats

START_BUFFER = 10

def PrintLine(str_locus, snp_locus, ld):
    for ldval in ld:
        if ldval is None: return
        ld_r, ld_pval = ldval[0]
        allele = ldval[1]
        freq = ldval[2]
        maf = ldval[3]
        ld_r2 = ld_r**2
        kldiv = ldval[4]
        sys.stdout.write("\t".join(map(str, [str_locus, snp_locus, allele, freq, maf, kldiv, ld_r2, ld_pval]))+"\n")
        sys.stdout.flush()

def CalcLD_r(str_record, snp_record, samples=[], str2=False, allele_r2=False, mincount=0, minmaf=0, usefilter=False):
    try:
        if min([snp_record.aaf[0], 1-snp_record.aaf[0]]) < minmaf: return [None] # Assume SNP biallelic
    except:
        return [None]
    if usefilter:
        if not(str_record.FILTER is None or str_record.FILTER == "PASS" or len(str_record.FILTER) == 0): return [None]
        if not(snp_record.FILTER is None or snp_record.FILTER == "PASS" or len(snp_record.FILTER) == 0): return [None]
    sample_to_gts = {}
    all_str_alleles = set()
    allele_counts = {}
    allele_counts2 = {} # if using str2
    if None in str_record.ALT: allelelens = [0]
    else: allelelens = [0] + [len(item)-len(str_record.REF) for item in str_record.ALT]
    for sample in str_record:
        if len(samples)>0 and sample.sample not in samples: continue
        sample_to_gts[sample.sample] = {"STR": None, "SNP": None}
        if None not in sample.gt_alleles:
            alleles = map(lambda x: allelelens[int(x)], sample.gt_alleles)
        else: continue
        for a in alleles:
            all_str_alleles.add(a)
            allele_counts[a] = allele_counts.get(a, 0) + 1
        sample_to_gts[sample.sample]["STR"] = alleles
    if None in snp_record.ALT: allelelens2 = [0]
    else: allelelens2 = [0] + [len(snp_record.ALT[i])-len(snp_record.REF) for i in range(len(snp_record.ALT))]
    for sample in snp_record:
        if sample.sample not in sample_to_gts.keys():
            continue
        if None not in sample.gt_alleles:
            if str2: # 2nd STR is imputed. Get diff from reference
                gt = map(int, sample.gt_alleles)
                for i in range(len(gt)):
                    a = allelelens2[gt[i]]
                    allele_counts2[a] = allele_counts2.get(a, 0) + 1
                sample_to_gts[sample.sample]["SNP"] = [allelelens2[gt[0]], allelelens2[gt[1]]]
            else:
                sample_to_gts[sample.sample]["SNP"] = map(int, sample.gt_alleles)
    kldiv = None
    maf = None
    if str2:
        kldiv = GetKLDivergence(allele_counts, allele_counts2, mincount)
    if allele_r2:
        ldresults = []
        for a in all_str_alleles:
            str_data = []
            snp_data = []
            if allele_counts[a] < mincount: continue
            for sample in sample_to_gts:
                if len(samples)>0 and sample not in samples: continue
                if sample_to_gts[sample]["STR"] is not None and sample_to_gts[sample]["SNP"] is not None:
                    # Enforce allele count
                    if allele_counts[sample_to_gts[sample]["STR"][0]] < mincount: continue
                    if allele_counts[sample_to_gts[sample]["STR"][1]] < mincount: continue
                    if str2 and allele_counts2[sample_to_gts[sample]["SNP"][0]] < mincount: continue
                    if str2 and allele_counts2[sample_to_gts[sample]["SNP"][1]] < mincount: continue
                    str_data.append(sum(map(lambda x: int(x==a), sample_to_gts[sample]["STR"])))
                    if str2:
                        snp_data.append(sum(map(lambda x: int(x==a), sample_to_gts[sample]["SNP"])))
                    else:
                        snp_data.append(sum(sample_to_gts[sample]["SNP"]))
            if len(snp_data) == 0: continue
            if str2:
                maf = allele_counts2.get(a, 0)*1.0/sum(allele_counts2.values())
            else:
                maf = min([snp_record.aaf[0], 1-snp_record.aaf[0]])
            ld = (scipy.stats.pearsonr(str_data, snp_data), a, allele_counts[a]*1.0/sum(allele_counts.values()), maf, kldiv)
            ldresults.append(ld)
        return ldresults
    else:
        str_data = []
        snp_data = []
        if str2:
            maf = GetHeterozygosity(allele_counts2, mincount=mincount)
        else: maf = min([snp_record.aaf[0], 1-snp_record.aaf[0]])
        for sample in sample_to_gts:
            if len(samples)>0 and sample not in samples: continue
            if sample_to_gts[sample]["STR"] is None or sample_to_gts[sample]["SNP"] is None: continue
            # Enforce allele counts
            if allele_counts[sample_to_gts[sample]["STR"][0]] < mincount: continue
            if allele_counts[sample_to_gts[sample]["STR"][1]] < mincount: continue
            # Get data
            str_data.append(sum(sample_to_gts[sample]["STR"]))
            snp_data.append(sum(sample_to_gts[sample]["SNP"]))
        het = GetHeterozygosity(allele_counts, mincount=mincount)
        if len(str_data) == 0:
            return []
        return [(scipy.stats.pearsonr(str_data, snp_data), "locus", het, maf, kldiv)]

def GetHeterozygosity(allele_counts, mincount=0):
    counts = []
    for a in allele_counts:
        if allele_counts[a] < mincount: continue
        counts.append(allele_counts[a])
    freqs = [item*1.0/sum(counts) for item in counts]
    return 1-sum([item**2 for item in freqs])

def CalcLD(str_reader, snp_reader, str_locus, snp_locus, use_info_start=False, samples=[], allele_r2=False, mincount=0, minmaf=0, usefilter=False):
    # Find STR record.
    chrom, start = str_locus.split(":")
    start = int(start)
    str_record = None
    if use_info_start:
        records = str_reader.fetch(chrom, start-START_BUFFER, start+START_BUFFER)
        for r in records:
            if r.INFO["START"] == start:
                str_record = r
                break
        if str_record is None:
            sys.stderr.write("ERROR: couldn't find STR record for %s\n"%start)
            return [None]
    else:
        records = str_reader.fetch(chrom, start, start+1)
        for r in records:
            str_record = r
            if r is None:
                sys.stderr.write("ERROR: couldn't find STR record for %s\n"%start)
                return [None]
            break
        if str_record is None or str_record.POS != start:
            sys.stderr.write("ERROR: couldn't find STR record for %s\n"%start)
            return [None]
    # Find SNP record
    snp_chrom, snp_start = snp_locus.split(":")
    snp_start = int(snp_start)
    records = snp_reader.fetch(snp_chrom, snp_start-1, snp_start)
    snp_record = None
    for r in records:
        snp_record = r
        break
    if snp_record is None:
        sys.stderr.write("Could not find SNP locus %s\n"%snp_locus)
        return [None]
    if snp_record.POS != snp_start:
        sys.stderr.write("ERROR: couldn't find SNP record for %s\n"%snp_start)
        return [None]
    return CalcLD_r(str_record, snp_record, samples=samples, allele_r2=allele_r2, mincount=mincount, minmaf=minmaf, usefilter=usefilter)

def GetKLDivergence(allele_counts, allele_counts2, mincount, pcount=1):
    p = []
    q = []
    alleles = set(allele_counts.keys()).union(set(allele_counts2.keys()))
    for a in alleles:
        acount = allele_counts.get(a, 0)
        if acount < mincount: acount = 0
        p.append(acount + pcount) # add pseudocount after removing outliers
        acount2 = allele_counts2.get(a, 0)
        if acount2 < mincount: acount2 = 0
        q.append(acount2 + pcount)
    psum = sum(p)
    qsum = sum(q)
    p = [item*1.0/psum for item in p]
    q = [item*1.0/qsum for item in q]
    return sum([p[i]*np.log(p[i]*1.0/q[i]) for i in range(len(p))])

def FindSnpLocus(snp_reader, str_locus, snp_locus_rsid, snp_region_start, snp_region_end):
    chrom, start = str_locus.split(":")
    start = int(start)
    try:
        records = snp_reader.fetch(chrom, snp_region_start, snp_region_end)
    except:
        sys.stderr.write("Could not fetch records for %s\n"%snp_locus_rsid)
        return None
    for r in records:
        if r.ID == snp_locus_rsid or str(r.POS) in snp_locus_rsid:
            return "%s:%s"%(r.CHROM, r.POS)
    return None

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--str-vcf", help="VCF with STR genotypes", type=str, required=True)
    parser.add_argument("--snp-vcf", help="VCF with SNP genotype", type=str, required=False)
    parser.add_argument("--str-vcf2", help="Compare two STR VCFs (e.g. true vs. imputed)", type=str, required=False)
    parser.add_argument("--str-locus", help="chr:start of STR locus", type=str, required=False)
    parser.add_argument("--snp-locus", help="chr:start of SNP locus", type=str, required=False)
    parser.add_argument("--snp-locus-rsid", help="rsid of SNP locus", type=str, required=False)
    parser.add_argument("--loci-file", help="File with chr,start(STR),start(SNP)", type=str, required=False)
    parser.add_argument("--loci-file-rsid", help="File with chr,start(STR),rsid. Optional columns snp_start, snp_end", type=str, required=False)
    parser.add_argument("--pairwise-snpstr", help="Calculate pairwise LD for all SNPs within maxdist of the STR", action="store_true")
    parser.add_argument("--max-dist", help="Don't consider snp/str more than this many bp apart", type=int, default=100000)
    parser.add_argument("--use-info-start", help="Match STR start on INFO/START (not POS)", action="store_true")
    parser.add_argument("--samples", help="Only consider samples in this file. (e.g. founders)", type=str, required=False)
    parser.add_argument("--allele-r2", help="Calculate r2 *per allele* rather than per locus", action="store_true")
    parser.add_argument("--mincount", help="Remove STR genotypes with an allele of count < this", type=int, default=0)
    parser.add_argument("--min-maf", help="Don't consider SNPs below this MAF. Only works properly when 2nd VCF is SNPs", type=float, default=0.0)
    parser.add_argument("--usefilter", help="Filter things not passing in VCF", action="store_true")
    parser.add_argument("--region", help="Restrict pairwise analysis to a certain region", type=str, required=False)
    args = parser.parse_args()

    # Output header
    sys.stdout.write("\t".join(["locus1","locus2","allele","freq_het","MAF","KL","r2","pval"])+"\n")
    # Open readers
    snp_reader = None
    str_reader2 = None
    str_reader = vcf.Reader(open(args.str_vcf, "rb"))
    if args.snp_vcf is not None:
        snp_reader = vcf.Reader(open(args.snp_vcf, "rb"))
    if args.str_vcf2 is not None:
        str_reader2 = vcf.Reader(open(args.str_vcf2, "rb"))

    # Get samples
    samples = set()
    if args.samples is not None:
        samples = set([item.strip() for item in open(args.samples, "r").readlines()])

    ###### Case 1: Single SNP/STR ##########
    if args.str_locus:
        str_locus = args.str_locus
        snp_locus = None
        if args.snp_locus:
            snp_locus = args.snp_locus
            snpid = snp_locus
            ld = CalcLD(str_reader, snp_reader, str_locus, snp_locus, \
                        use_info_start=args.use_info_start, samples=samples, \
                        allele_r2=args.allele_r2, mincount=args.mincount,
                        minmaf=args.min_maf, usefilter=args.usefilter)
        elif args.snp_locus_rsid:
            snpid = args.snp_locus_rsid
            str_start = int(str_locus.split(":")[1])
            snp_region_start = str_start - args.max_dist
            snp_region_end = str_start + args.max_dist
            snp_locus = FindSnpLocus(snp_reader, args.str_locus, args.snp_locus_rsid, \
                                     snp_region_start, snp_region_end)
            ld = CalcLD(str_reader, snp_reader, args.str_locus, snp_locus, \
                        use_info_start=args.use_info_start, samples=samples, \
                        allele_r2=args.allele_r2, mincount=args.mincount, minmaf=args.min_maf, usefilter=args.usefilter)
            if snp_locus is None:
                sys.stderr.write("ERROR: Couldn't find SNP locus within %s of %s\n"%(args.max_dist, args.str_locus))
        else:
            sys.stderr.write("ERROR: No SNP locus specified. Use --snp-locus or --snp-locus-rsid\n")
        PrintLine(str_locus, snpid, ld)

    ###### Case 2: List of pairs from a file ##########
    # Keep track of snp id positions in case we see same one many times
    snp_rsid_to_pos = {}
    if args.loci_file_rsid is not None or args.loci_file is not None:
        has_rsid = False
        fname = args.loci_file
        if args.loci_file_rsid is not None:
            fname = args.loci_file_rsid
            has_rsid = True
        with open(fname, "r") as f:
            for line in f:
                items = line.strip().split()
                chrom = items[0]
                str_start = int(items[1])
                str_locus = "%s:%s"%(chrom, str_start)
                snp_loci = []
                if has_rsid:
                    rsids = items[2].split(",")
                    try:
                        snp_region_start = int(items[3])
                    except IndexError:
                        snp_region_start = str_start-max_dist
                    try:
                        snp_region_end = int(items[4])
                    except IndexError:
                        snp_region_end = str_start+max_dist
                    for rsid in rsids:
                        snp_locus = None
                        try:
                            snp_locus = snp_rsid_to_pos[rsid]
                            if snp_locus is None: continue
                        except KeyError:
                            snp_locus = FindSnpLocus(snp_reader, str_locus, rsid, \
                                                     snp_region_start, snp_region_end)
                            snp_rsid_to_pos[rsid] = snp_locus
                        if snp_locus is None:
                            sys.stderr.write("Could not find %s\n"%rsid)
                            snp_rsid_to_pos[rsid] = None
                            continue
                        snp_loci.append(snp_locus)
                else:
                    snp_start = int(items[2])
                    snp_loci = ["%s:%s"%(chrom, item) for item in snp_start.split(",")]
                for i in range(len(snp_loci)):
                    snp_locus = snp_loci[i]
                    snpid = snp_locus
                    if args.loci_file_rsid is not None:
                        snpid = rsids[i]
                    ld = CalcLD(str_reader, snp_reader, str_locus, snp_locus, \
                                use_info_start=args.use_info_start, \
                                samples=samples, allele_r2=args.allele_r2, mincount=args.mincount, minmaf=args.min_maf, usefilter=args.usefilter)
                    PrintLine(str_locus, snpid, ld)

    ###### Case 3: All SNP-STR pairwise ##########
    if args.pairwise_snpstr:
        if args.region != None:
            str_records = str_reader.fetch(args.region)
        else: str_records = str_reader
        for str_record in str_records:
            if args.usefilter:
                if not(str_record.FILTER is None or str_record.FILTER == "PASS" or len(str_record.FILTER) == 0): continue
            str_locus = "%s:%s"%(str_record.CHROM, str_record.POS)
            region_start = max([0, str_record.POS - args.max_dist])
            region_end = str_record.POS + args.max_dist
            try:
                snp_records = snp_reader.fetch(str_locus.split(":")[0], region_start, region_end)
            except ValueError:
                sys.stderr.write("ERROR fetching SNP records for STR locus %s\n"%str_locus)
                continue
            for snp_record in snp_records:
                snp_locus = "%s:%s"%(snp_record.CHROM, snp_record.POS)
                if snp_locus == str_locus: continue
                ld = CalcLD_r(str_record, snp_record, samples=samples, allele_r2=args.allele_r2, mincount=args.mincount, minmaf=args.min_maf, usefilter=args.usefilter)
                PrintLine(str_locus, snp_locus, ld)

    ###### Case 4: Compare two STR VCFs ######
    if str_reader2 is not None:
        str_records = str_reader
        if args.region != None:
            str_records = str_reader.fetch(args.region)
        for str_record in str_records:
            if len(str_record.REF) == 1: continue # SNP
            if args.usefilter and not(str_record.FILTER is None or str_record.FILTER == "PASS" or len(str_record.FILTER) == 0): continue
            str_locus = "%s:%s"%(str_record.CHROM, str_record.POS)
            str_record2 = None
            records = str_reader2.fetch(str_record.CHROM, str_record.POS-START_BUFFER, str_record.POS+1+START_BUFFER)
            for r in records:
                if r.ID == str_record.ID:
                    str_record2 = r
                    break
            if str_record2 is None:
                continue
            str_locus2 = "%s:%s"%(str_record2.CHROM, str_record2.POS)
            ld = CalcLD_r(str_record, str_record2, samples=samples, str2=True, allele_r2=args.allele_r2, mincount=args.mincount, minmaf=args.min_maf, usefilter=args.usefilter)
            PrintLine(str_locus, str_locus2, ld)

if __name__ == "__main__":
    main()