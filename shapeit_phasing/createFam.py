### python create_fam.py fam160Detail.txt chr20.fam chr20.fam.txt

import sys
import pandas as pd
import numpy as np



header = ["fam","sam","fid","mid","sex","pheno"]
famDetail = pd.read_csv(sys.argv[1],names=header,delim_whitespace=True)
famFile = pd.read_csv(sys.argv[2],names=header,delim_whitespace=True)

famNP = famFile.values
for i in range(famNP.shape[0]):
    if famNP[i][1] in famDetail['sam'].values:
        famNP[i][0] = famDetail[famDetail['sam']==famNP[i][1]].fam.values[0]
        famNP[i][2] = famDetail[famDetail['sam']==famNP[i][1]].fid.values[0]
        famNP[i][3] = famDetail[famDetail['sam']==famNP[i][1]].mid.values[0]
        famNP[i][4] = famDetail[famDetail['sam']==famNP[i][1]].sex.values[0]


df = pd.DataFrame(famNP)
df.to_csv(sys.argv[3],sep=" ",header=False,index=False)
print ".fam file updated"