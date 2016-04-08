#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#plt.rcParams['backend'] = 'TkAgg'

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
# hide first label
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as plticker

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

def thousands(x, pos):
    return '{0:,}'.format(int(x/1000))

#for mode in 'snps,mnps indels snps snps,mnps,indels,other mnps,indels other,indels snps,mnps,indels'.split():
#for mode in ('snps,mnps',):
#for mode in 'snps,mnps indels snps'.split():

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

    with open('out_vcfeval/{}/vcfeval.log'.format(mode)) as f:
        baseline = 0
        for line in f:
            match = re.match(r'.+baseline contains ([0-9]+) variants', line)
            if match:
                baseline += int(match.group(1))

    l_minVQSLod = list(sorted(d_minVQSLod2targetTruthSensitivity.keys()))
    d_minVQSLod2xy = {}
    minVQSLod = l_minVQSLod.pop()
    with gzip.open(
        'out_vcfeval/{}/weighted_roc.tsv.gz'.format(mode), 'rt') as f:
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
#                print('targetTruthSensitivity', d_minVQSLod2targetTruthSensitivity[minVQSLod], fp, tpr)

#    fn = 0
#    with gzip.open('out_vcfeval/{}/fn.vcf.gz'.format(mode)) as f:
#        for line in f:
#            if line[0] != '#':
#                fn += 1

#    print('fn', fn, 'tp', tp, 'fn+tp', fn+int(float(tp)), 'baseline', baseline)
#    assert fn+int(float(tp)) == baseline

    print('tp', tp, 'fp', fp, 'tp+fp', tp+fp, 'tpr', tpr, 'fdr', fp/(tp+fp), 'baseline', baseline)

    ## fp2fdr
    x = [_/(tp+fp) for _ in x]
    for minVQSLod in d_minVQSLod2xy.keys():
        d_minVQSLod2xy[minVQSLod][0] /= tp+fp

#    def fdr2fp(fdr, pos):
#        return str(int(fdr*(fp+tp)))  #0:fp+tp
#        return str(int(fdr*fp))  #0:fp
#    formatter = FuncFormatter(fdr2fp)
#
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax3 = ax2.twinx()

    plt.subplots_adjust(top=0.9, right=0.75)

    ax2.set_xlabel('False calls (thousands)')
    ax3.set_ylabel('True calls (thousands)')
    ax1.set_xlabel(r'FDR = FP/(TP+FP)')
    ax1.set_ylabel(r'TPR = TP/P = TP/(TP+FN)')

    ax1.legend(loc='lower right', shadow=True)

    ax1.set_ylim([0,1])
    ax1.set_xlim([0,fp/(tp+fp)])
    ax2.set_xlim(0,fp)
    ax3.set_ylim(0,baseline*ax1.get_ylim()[1])
    ## set_major_formatter
    majorFormatter = FuncFormatter(thousands)
    ax2.xaxis.set_major_formatter(majorFormatter)
    ax3.yaxis.set_major_formatter(majorFormatter)
    ## set_major_locator
    ax1.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax1.yaxis.set_major_locator(MaxNLocator(prune='both'))
    ax2.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax3.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    ## Many true positive SNPs...
    if 'snps' in mode:
        ax3.yaxis.set_major_locator(MultipleLocator(500000))
    else:
        ax3.yaxis.set_major_locator(MultipleLocator(100000))
    ## Many false positive indels...
    if 'indels' in mode:
        ax2.xaxis.set_major_locator(MultipleLocator(100000))
        ax1.xaxis.set_major_locator(MultipleLocator(0.1))
    else:
        ax2.xaxis.set_major_locator(MultipleLocator(100000))
        ax1.xaxis.set_major_locator(MultipleLocator(0.02))
#    ax1.xaxis.set_major_locator(MaxNLocator(prune='lower'))
#    ax1.yaxis.set_major_locator(MaxNLocator(prune='lower'))
#    ax2.xaxis.set_major_locator(MaxNLocator(prune='lower'))
#    ax3.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    ##
    ax2.tick_params(direction='out')
    ax3.tick_params(direction='out')
    for label in ax3.yaxis.get_ticklabels():
        label.set_horizontalalignment('right')
    if 'snps' in mode:
        ax3.tick_params(pad=40)
    else:
        ax3.tick_params(pad=30)
    ## grid
    ax1.grid(b=True, which='major', color='k', linestyle='--')
    ax1.plot(x, y, label=d_labels[mode], linewidth=3)
    for minVQSLod in d_minVQSLod2xy.keys():
        x, y = d_minVQSLod2xy[minVQSLod]
        targetTruthSensitivity = d_minVQSLod2targetTruthSensitivity[minVQSLod]
        ax1.annotate(
            targetTruthSensitivity, xy=(x, y),
            xytext=(-0.1, 0.1), textcoords='offset points',
            ha='right', va='bottom', fontsize='xx-small')
    x, y = zip(*d_minVQSLod2xy.values())
    ax1.plot(x, y, label=d_labels[mode], linestyle='', marker='o', color='r')

    plt.savefig('fdr_{}.png'.format(mode), dpi=600)
    plt.close()
    plt.clf()

    print(max(y))
