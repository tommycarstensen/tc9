import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

path = 'finestructure_results.txt'

a = np.genfromtxt(path, delimiter='\t', usecols=(1,2,3,4,5))
nrows = a.shape[0]

labelsy = np.genfromtxt(path, delimiter='\t', usecols=(0), dtype=str)

labels = ('TSI', 'CHB', 'Hadza', '''Ju/'hoan North''', 'Mbuti Pygmy')
colors = ('#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00')

#daysofweek = {1:('Sunday','r'), 
#              2:('Monday','g'), 
#              3:('Tuesday','b'), 
#              4:('Wednesday','yellow'), 
#              5:('Thursday','k'), 
#              6:('Friday', 'magenta'), 
#              7:('Saturday', 'orange')}

figs, axs = plt.subplots(nrows, 1)

minorLocator = MultipleLocator(0.01)

for irow, ax in zip(range(nrows), axs):
    for x, (y, color, label) in enumerate(zip(a[irow], colors, labelsy)):
        axs[irow].bar(
            x, y, width=0.8, bottom=0.0, align='center', color=color,
            alpha=0.9, label=label)
#        ax.set_xticks(False)
        ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
        ax.tick_params(axis='y', which='both', labelsize='xx-small')
#        ax.yaxis.label(labelsy[irow], fontsize='small')
#        ax.yaxis.set_title(labelsy[irow], fontsize='small')
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(labelsy[irow], fontsize='x-small', rotation='horizontal', horizontalalignment='left')
        ax.yaxis.set_minor_locator(minorLocator)
        ax.yaxis.grid(True, which='minor')
        ax.set_yticks((0, 0.05))
#        ax.tick_params(axis='y', which='minor', left='on')

axs[-1].tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
axs[-1].set_xticks(range(5))
axs[-1].set_xticklabels(labels, fontsize='small')
for ax in axs:
    ax.spines['top'].set_visible(False)
#axs[-1].set_xticks(label)
#ax1.set_xticklabels([daysofweek[i][0] for i in xval])
#ax1.legend()
plt.savefig('finestructure_barplot.png', dpi=600, bbox_inches='tight')
