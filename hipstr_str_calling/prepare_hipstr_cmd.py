####
### Not needed anymore
####


import os
import sys

bamFiles = []

chrm = []
startPos = []
endPod = []
samToolsIn = []

f = open(sys.argv[1])
for l in f.readlines():
    a,b,c = l.strip().split(" ")
    chrm.append(a)
    startPos.append(b)
    endPod.append(c)
    samToolsIn.append(a+":"+b+"-"+c)
f.close()

regionsSamTools = " ".join(samToolsIn)

for file in os.listdir("bam"):
    if file.endswith(".bam"):
        bamFiles.append(file)

bamFilesCat = "bam/"+",bam/".join(bamFiles)

delCmd = "rm bedFile.bed"
os.popen(delCmd).read()
for i in range(len(startPos)):
    commd = 'cat HipSTR_reference.hg19.bed | grep "chr'+chrm[i]+'\t'+startPos[i]+'"'
    bedOutput = " ".join(os.popen(commd).read().strip().split("\t"))
    with open("bedFile.bed", "a") as myfile:
        myfile.write(bedOutput+"\n")

delCmd = "rm hipstrCommand.sh"
os.popen(delCmd).read()
with open("hipstrCommand.sh", "a") as myfile:
        myfile.write("#!/bin/sh"+"\n\n")

numSample = 160
for i in range(len(bamFiles)/numSample + 1):
    bamFilesCat = "bam/"+",bam/".join(bamFiles[(i*numSample):((i+1)*numSample)])
    if(len(bamFilesCat) > 4):
        hipstrCmd = sys.argv[3]+" --bams "+ bamFilesCat + " --fasta ../fasta/ --regions bedFile.bed --str-vcf str_calls"+str(i)+".vcf.gz --snp-vcf "+sys.argv[2]+"\n"
        with open("hipstrCommand.sh", "a") as myfile:
            myfile.write(hipstrCmd)

commd = "chmod 777 hipstrCommand.sh"
os.popen(commd).read()
