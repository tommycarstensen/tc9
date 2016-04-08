## http://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
# hide first label
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as plticker
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt

import re
import gzip
import numpy as np

d_labels = {
    'snps,mnps':'SNPs', 'indels':'Indels',
    'snps':'snps', 'snps,mnps,indels,other':'snps,mnps,indels,other', 'mnps,indels':'mnps,indels',
    'other,indels':'other,indels', 'snps,mnps,indels':'snps,mnps,indels',
}
d_tranches = {
    'snps,mnps':'../pipeline_UG3.3_ug2g_agv_NA12878/out_VariantRecalibrator/SNP.tranches',
    'indels':'../pipeline_UG3.3_ug2g_agv_NA12878/out_VariantRecalibrator/INDEL.tranches',
    'snps':'../pipeline_UG3.3_ug2g_agv_NA12878/out_VariantRecalibrator/SNP.tranches',
    'snps,mnps,indels,other':'../pipeline_UG3.3_ug2g_agv_NA12878/out_VariantRecalibrator/INDEL.tranches',
    'mnps,indels':'../pipeline_UG3.3_ug2g_agv_NA12878/out_VariantRecalibrator/INDEL.tranches',
    'other,indels':'../pipeline_UG3.3_ug2g_agv_NA12878/out_VariantRecalibrator/INDEL.tranches',
    'snps,mnps,indels':'../pipeline_UG3.3_ug2g_agv_NA12878/out_VariantRecalibrator/INDEL.tranches',
    }
subdir = 'pipeline_HC3.2'
d_tranches = {
    'snps,mnps':'../{}/out_VariantRecalibrator/SNP.tranches'.format(subdir),
    'indels':'../{}/out_VariantRecalibrator/INDEL.tranches'.format(subdir),
}

def thousands(x, pos):
    return '{0:,}'.format(int(x/1000))

#for mode in 'snps,mnps indels snps snps,mnps,indels,other mnps,indels other,indels snps,mnps,indels'.split():
#for mode in ('snps,mnps',):
#for mode in 'snps,mnps indels snps'.split():

fig = plt.figure()
ax1 = host_subplot(111, axes_class=AA.Axes)
## copy x-axis
ax2 = ax1.twiny()
ax2b = ax1.twiny()
## copy y-axis
ax3 = ax1.twinx()
ax3b = ax1.twinx()

offset=30
offset=35
new_fixed_axis = ax2b.get_grid_helper().new_fixed_axis
ax2b.axis['top'] = new_fixed_axis(loc="top", axes=ax2b, offset=(0, offset))
offset=55
offset=60
new_fixed_axis = ax3b.get_grid_helper().new_fixed_axis
ax3b.axis['right'] = new_fixed_axis(loc="right", axes=ax3b, offset=(offset, 0))
# toggle on ticks and labels
ax2b.axis['top'].toggle(all=True)
ax3b.axis['right'].toggle(all=True)

plt.subplots_adjust(top=0.80, right=0.80)

ax2.set_xlabel('False SNP calls (thousands)', color='b')
ax3.set_ylabel('True SNP calls (thousands)', color='b')
ax2b.set_xlabel('False indel calls (thousands)', color='g')
ax3b.set_ylabel('True indel calls (thousands)', color='g')
ax1.set_xlabel(r'FDR = FP/(TP+FP)')
ax1.set_ylabel(r'TPR = TP/P = TP/(TP+FN)')

l_baseline = []
l_fp = []
l_tp = []
for mode in 'snps,mnps indels'.split():

    print(mode)

    d_minVQSLod2targetTruthSensitivity = {}
    with open(d_tranches[mode]) as f:
        for line in f:
            if line[0] == '#':
                continue
            l = line.split(',')
            if l[0] == 'targetTruthSensitivity':
                index = l.index('minVQSLod')
                continue
            targetTruthSensitivity = l[0]
            minVQSLod = float(l[index])
            if not targetTruthSensitivity in ('100.00', '99.90', '99.50', '99.00', '90.00'):
                continue
            d_minVQSLod2targetTruthSensitivity[minVQSLod] = targetTruthSensitivity
#            print('targetTruthSensitivity', targetTruthSensitivity)

    with open('out_vcfeval/{}/{}/vcfeval.log'.format(subdir, mode)) as f:
        baseline = 0
        for line in f:
            match = re.match(r'.+baseline contains ([0-9]+) variants', line)
            if match:
                baseline += int(match.group(1))
    l_baseline.append(baseline)

    l_minVQSLod = list(sorted(d_minVQSLod2targetTruthSensitivity.keys()))
    d_minVQSLod2xy = {}
    minVQSLod = l_minVQSLod.pop()
    with gzip.open(
        'out_vcfeval/{}/{}/weighted_roc.tsv.gz'.format(subdir, mode), 'rt') as f:
        x = [0]
        y = [0]
        for line in f:
            if line[0] == '#':
                continue
            score, tp, fp = map(float, line.rstrip().split('\t'))
            tpr = tp/baseline
            x.append(fp)
            y.append(tpr)
            if score <= minVQSLod:
                d_minVQSLod2xy[minVQSLod] = [x[-1], y[-1]]
                minVQSLod = l_minVQSLod.pop()
                print('targetTruthSensitivity', d_minVQSLod2targetTruthSensitivity[minVQSLod], fp, tpr)
    l_fp.append(fp)
    l_tp.append(tp)

    print('tp', tp, 'fp', fp, 'tp+fp', tp+fp, 'tpr', tpr, 'fdr', fp/(tp+fp), 'baseline', baseline)

    ## fp2fdr
    x = [_/(tp+fp) for _ in x]
    for minVQSLod in d_minVQSLod2xy.keys():
        d_minVQSLod2xy[minVQSLod][0] /= tp+fp

    ax1.plot(x, y, label=d_labels[mode], linewidth=3)
    for minVQSLod in d_minVQSLod2xy.keys():
        x, y = d_minVQSLod2xy[minVQSLod]
        targetTruthSensitivity = d_minVQSLod2targetTruthSensitivity[minVQSLod]
        ax1.annotate(
            targetTruthSensitivity, xy=(x, y),
            xytext=(-0.1, 0.1), textcoords='offset points',
            ha='right', va='bottom', fontsize='xx-small')
    x, y = zip(*d_minVQSLod2xy.values())
    ax1.plot(x, y, linestyle='', marker='o', color='r')

## colors and alignment
for label in ax3.yaxis.get_ticklabels():
    label.set_horizontalalignment('right')
    label.set_color('b')
for label in ax3b.yaxis.get_ticklabels():
    label.set_horizontalalignment('right')
    label.set_color('g')

for tick in ax3.yaxis.get_major_ticks():
    tick.label1On = True
    tick.label1.set_color('g')
    tick.label2On = True
    tick.label2.set_color('g')
ax3b.yaxis.label.set_color('g')
ax3b.axis['right'].label.set_color('g')

for tl in ax3b.get_yticklabels():
   tl.set_color('g')

## legend
ax1.legend(loc='lower right', shadow=True)

ax1.set_ylim([0,1])
ax1.set_xlim([0,max([l_fp[i]/(l_tp[i]+l_fp[i]) for i in range(2)])])
ax2.set_xlim(0,l_fp[0])
ax3.set_ylim(0,l_baseline[0])
ax2b.set_xlim(0,l_fp[1])
ax3b.set_ylim(0,l_baseline[1])

ax1.grid(b=True, which='major', color='k', linestyle='--')

## set_major_locator
ax1.xaxis.set_major_locator(MaxNLocator(prune='lower'))
ax1.yaxis.set_major_locator(MaxNLocator(prune='both'))
ax2.xaxis.set_major_locator(MaxNLocator(prune='lower'))
ax3.yaxis.set_major_locator(MaxNLocator(prune='lower'))
ax3b.yaxis.set_major_locator(MaxNLocator(prune='lower'))
ax3.yaxis.set_major_locator(MultipleLocator(500000))
ax3b.yaxis.set_major_locator(MultipleLocator(100000))
ax2.xaxis.set_major_locator(MultipleLocator(100000))
ax2b.xaxis.set_major_locator(MultipleLocator(100000))
ax1.xaxis.set_major_locator(MultipleLocator(0.05))

## set_major_formatter
majorFormatter = FuncFormatter(thousands)
ax2.xaxis.set_major_formatter(majorFormatter)
ax3.yaxis.set_major_formatter(majorFormatter)
ax2b.xaxis.set_major_formatter(majorFormatter)
ax3b.yaxis.set_major_formatter(majorFormatter)

## tick_params
#ax2.tick_params(direction='out', colors='b')
#ax3.tick_params(direction='out', colors='b', pad=40)
#ax2b.tick_params(direction='out', colors='g')
#ax3b.tick_params(direction='out', colors='g', pad=30)

plt.savefig('fdr.png', dpi=600)
plt.close()
plt.clf()
print(max(y))
