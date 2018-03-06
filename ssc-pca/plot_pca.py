from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pylab
import matplotlib

header = ["sam","fam"]+["pca"+str(i) for i in range(10)]
pca_eigvec = pd.read_csv("pca_10.eigenvec", sep=" ",names=header)

label = ['AFR']*70 + ['AMR']*50 + ['EAS']*50 + ['EUR']*50 + ['SAS']*50 + ['SSC']*2076
colors = {'AFR':'brown', 'AMR':'black', 'EAS':'green', 'EUR':'blue', 'SAS':'purple', 'SSC':'red'}

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

#handles = [handles[0], handles[2], handles[3], handles[4], handles[1], handles[5]]
#labels = [labels[0], labels[2], labels[3], labels[4], labels[1], labels[5]]

leg=ax.legend(handles, labels, numpoints=1, loc='best', fancybox=True, framealpha=0.5)
#leg.get_frame().set_alpha(0.5)

plt.savefig("pca.pdf")

