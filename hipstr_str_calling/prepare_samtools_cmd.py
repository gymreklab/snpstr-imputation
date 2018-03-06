###
### Usage: python prepare_samtools_cmd.py bedFile.bed famID.txt /samtools/samtools.exe
###

import sys
import os

delCmd = "rm samToolsCommand.sh"
os.popen(delCmd).read()

chrm = []
startPos = []
endPod = []
samToolsIn = []

f = open(sys.argv[1])
for l in f.readlines():
    a,b,c = l.strip().split()[0:3]
    d = str(filter(str.isdigit, a))
    b = str(int(b)-5000)
    c = str(int(c)+5000)
    chrm.append(d)
    startPos.append(b)
    endPod.append(c)
    samToolsIn.append(d+":"+b+"-"+c)
f.close()

regionsSamTools = " ".join(samToolsIn)

f = open(sys.argv[2])
famID = []
samID = []
for l in f.readlines():
    a,b=l.strip().split()
    famID.append(a)
    samID.append(b)
f.close()

samToolsCommands = []
for i in range(len(famID)):
    cmdd = sys.argv[3]+" view -b s3://sscwgs/"+famID[i]+"/BAM/Sample_"+samID[i]+"/analysis/"+samID[i]+".final.bam "+regionsSamTools+" > bam/"+samID[i]+".bam"
    samToolsCommands.append(cmdd)

with open("samToolsCommand.sh", "a") as myfile:
        myfile.write("#!/bin/sh"+"\n\n")

for i in range(len(samToolsCommands)):
    with open("samToolsCommand.sh", "a") as myfile:
        myfile.write(samToolsCommands[i]+"\n")


commd = "chmod 777 samToolsCommand.sh"
os.popen(commd).read()
