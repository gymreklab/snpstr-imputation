#!/usr/bin/env python

import sys

def GetMend(child, mother, father):
    child_gt = child.split("|")
    mother_gt = mother.split("|")
    father_gt = father.split("|")
    mend = False
    if child_gt[0] in mother_gt and child_gt[1] in father_gt: mend = True
    elif child_gt[1] in mother_gt and child_gt[0] in father_gt: mend = True
    else: mend = False
    allref = (mother == "0|0" and father == "0|0" and child == "0|0")
    return mend, allref

line = sys.stdin.readline()
while line != "":
    items = line.strip().split()
    if "." in items:
        line = sys.stdin.readline()
        continue
    newitems = []
    newitems.append(items[0]) # chrom
    newitems.append(items[1]) # start
    newitems.append(items[2]) # child
    newitems.append(min([float(x) for x in [items[4],items[8],items[12]]])) # min Q
    newitems.append(min([float(x) for x in [items[5],items[9],items[13]]])) # min DP
    ismend, allref = GetMend(items[3], items[7], items[11])
    newitems.append(ismend)
    newitems.append(allref)
    sys.stdout.write("\t".join(map(str, newitems))+"\n")
    line = sys.stdin.readline()
