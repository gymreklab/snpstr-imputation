from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pylab
import matplotlib
import sys

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

header = ["sample","family"]+["pca"+str(i) for i in range(10)]
pca_eigvec = pd.read_csv("pca_10.eigenvec", delim_whitespace= True, names=header)
sample_pop = pd.read_csv("integrated_call_samples_v3.20130502.ALL.panel", delim_whitespace=True)

pca_eigvec = pca_eigvec.merge(sample_pop, on = "sample", how="left")
pca_eigvec.fillna('SSC', inplace=True)

label = list(pca_eigvec['super_pop'])
cmap = get_cmap(np.unique(label).shape[0])

colors = dict([(np.unique(label)[i], np.random.rand(3,1)) for i in range(np.unique(label).shape[0])])

eigvec = pca_eigvec.values[:,2:12]

df = pd.DataFrame(dict(x=eigvec[:,0], y=eigvec[:,1], label=label))

groups = df.groupby('label')

# Plot
fig, ax = plt.subplots()
ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
for name, group in groups:
    ax.plot(group.x, group.y, marker='.', linestyle='', ms=10, label=name, color=colors[name])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('PC1')
plt.ylabel('PC2')

handles,labels = ax.get_legend_handles_labels()

leg=ax.legend(handles, labels, numpoints=1, loc='best', fancybox=True, framealpha=0.5)

plt.savefig("pca.pdf")

