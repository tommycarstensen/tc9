import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
import itertools

path = '/lustre/scratch114/projects/uganda_gwas/users/dg11/f2/split2/f2_freq_substable_norm.txt'

labels = pops = ("GWD", "MSL", "ESN","YRI", "Baganda", "Banyarwanda", "RwandeseUgandan", "Barundi", "Banyankole", "LWK", "Amhara", "Oromo", "Somali", "Wolayta", "Zulu", "ACB", "ASW", "CEU", "GBR", "TSI", "IBS", "FIN", "ITU", "STU", "GIH", "PJL", "MXL", "PEL", "PUR", "CLM",  "BEB", "CDX", "CHB", "CHS" , "JPT", "KHV")

colors = ('#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf')
multiples = (4, 10, 1, 2, 5, 4, 5, 5)
assert sum(multiples) == len(pops)
# generate list of list of colors and flatten it
colors = list(itertools.chain(*[ multiples[i] * [colors[i]] for i in range(len(colors)) ]))

n = len(pops)
a = np.empty((n, n))
with open(path) as f:
    for line in f:
        l = line.rstrip().split()
        cnt = int(l[2])
        i1 = pops.index(l[1].split(':')[0])
        i2 = pops.index(l[1].split(':')[1])
        a[i1][i2] = cnt
        a[i2][i1] = cnt

## normalize
for i in range(n):
#    a[:,i] /= max(a[:,i])
    a[i] /= max(a[i])

figs, axs = plt.subplots(n, 1)

minorLocator = MultipleLocator(2000)
minorLocator = MultipleLocator(0.2)

for i, ax in enumerate(axs):
#    for x, (y, color, labelx) in enumerate(zip(a[i], colors, pops)):
    axs[i].bar(
        np.arange(n), a[i], width=0.8, bottom=0.0, align='center', color=colors,
        alpha=0.9, edgecolor=None, linewidth=0)
#        ax.set_xticks(False)
    ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    ax.set_xticklabels(pops)
    labely = pops[i]
    ax.tick_params(
        axis='y', which='both', labelsize='xx-small',
        left='off', right='off')
#        ax.yaxis.label(labelsy[irow], fontsize='small')
#        ax.yaxis.set_title(labelsy[irow], fontsize='small')
    ax.yaxis.set_label_position("right")
    ax.set_ylabel(labely, fontsize='x-small', rotation='horizontal', horizontalalignment='left')
#    ax.yaxis.set_minor_locator(minorLocator)
#    ax.yaxis.grid(True, which='minor')
#    ax.set_yticks((10000,))
    ax.set_yticks(())
#    ax.set_ylim((0, 15000))
    ax.set_ylim((0, 1))
    ax.set_xlim((-0.5, n-0.5))
#        ax.tick_params(axis='y', which='minor', left='on')
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

axs[-1].tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
axs[-1].set_xticks(range(len(labels)))
axs[-1].set_xticklabels(pops, fontsize='small', rotation=90)
for ax in axs:
    ax.spines['top'].set_visible(False)
#axs[-1].set_xticks(label)
#ax1.set_xticklabels([daysofweek[i][0] for i in xval])
#ax1.legend()

plt.savefig('f2_freq_substable_norm.png', dpi=600, bbox_inches='tight')
