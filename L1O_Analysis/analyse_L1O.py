import pandas as pd
import numpy as np
from scipy.stats.stats import pearsonr
import sys

dir = sys.argv[1]
files = sys.argv[2]
pos_file = sys.argv[3]
output = sys.argv[4]

sampleHead = ["sampleID"]
samples = pd.read_csv(dir+"/"+files,names=sampleHead)

header = ['str','gral1','gral2','imal1','imal2']
mergeData = pd.read_csv(dir+"/"+samples.values[0][0],names=header,delim_whitespace=True)

for i in samples.values[1:]:
    data = pd.read_csv(dir+"/"+i[0],names=header,delim_whitespace=True)
    mergeData = mergeData.append(data)

sumMerged = pd.DataFrame({'str':mergeData['str'], 'im':(mergeData['imal1'].values + mergeData['imal2'].values), 'gr':(mergeData['gral1'].values + mergeData['gral2'].values)})

### remove calls with allele counts < 3
counts = sumMerged.groupby(['str','gr']).count().reset_index()
counts.columns = ['str','gr','counts']

sumMerged = sumMerged.merge(counts, how="inner", on=["str","gr"])
sumMerged = sumMerged[sumMerged.counts>=3]

def rVal(group, col1, col2):
    c1 = group[col1]
    c2 = group[col2]
    return pd.Series({'r': pearsonr(c1, c2)[0], 'pVal': pearsonr(c1, c2)[1]})

corr = sumMerged.groupby('str').apply(rVal, 'gr', 'im').reset_index()

### find number of samples
numSamples = sumMerged.groupby(['str']).count().reset_index()[['str','counts']]
numSamples.columns = ['str','numSamples']

droppedNa = mergeData.dropna(axis=0)
concord = list()
for i in droppedNa.values:
    listA = set(i[1:3].astype(int))
    listB = set(i[3:5].astype(int))
    concord.append( (2-(max(len(listA-listB) , len(listB-listA))))/2.0 )

concordance = pd.DataFrame({'str':droppedNa['str'], 'concord':concord})
concordance = concordance.groupby('str').mean().reset_index()


idHeader = ['str', 'pos']
idData = pd.read_csv(dir+"/"+pos_file,names=idHeader,delim_whitespace=True)

finalData = corr.merge(concordance, how="inner", on="str")
finalData = finalData.merge(numSamples, how="inner",on="str")
finalData = finalData.merge(idData, how="inner",on="str")

finalData.to_csv(dir+"/"+output, index=False)