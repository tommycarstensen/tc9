import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import matplotlib
import re
import gzip
#import matplotlib.font_manager as font_manager
#from matplotlib import pyplot as plt, dates as mdates
import numpy as np

l_colors = [
    '#809CCE',  # dark pastel blue
    '#42C312',  # dark pastel green
    '#9668DC',  # dark pastel purple
    '#CFCFC3',  # pastel gray
    '#000000',  # black
    ]

l_colors = [
'#7fc97f',
'#beaed4',
'#fdc086',
'#ffff99',
'#386cb0',
]

l_colors = [
'#b2182b',
'#ef8a62',
'#fddbc7',
'#999999',
'#4d4d4d',
]

l_colors = [
'#e41a1c',
'#377eb8',
'#4daf4a',
'#984ea3',
'#ff7f00',
'#000000',
]

l_colors = [
'#e41a1c',
'#377eb8',
'#4daf4a',
'#984ea3',
]

l_colors = [
'#a6cee3',
'#1f78b4',
'#b2df8a',
'#33a02c',
'#fb9a99',
'#e31a1c',
'#fdbf6f',
'#ff7f00',
'#cab2d6',
'#6a3d9a',
'#ffff99',
'#b15928',
]

plt.rcParams['backend'] = 'TkAgg' 
#plt.rcParams['backend'] = 'Agg' 

type_plot = 'FDR'

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
#ax1.set_yticklabels([0,.2,.4,.6,.8,1.])
#ax2.set_yticklabels([10000,12000,14000,15000,16000])

for mode in ('snps', 'indels'):
#    plt.rc('font', family='sans-serif')
#    plt.rc('font', serif='LiberationSans-Regular')
#    plt.rcParams['font.family'] = 'LiberationSans-Regular'
#
#    path='/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf'
#    prop = font_manager.FontProperties(fname=path)
#    plt.rcParams['font.family'] = prop.get_name()
    print(mode)
    for i, affix in enumerate((
#        'NIST.samtools.{}.QUAL'.format(mode),
#        'NIST.UG.{}.QUAL'.format(mode),
#        'NIST.FB.{}.QUAL'.format(mode),
#        'NIST.HC.{}.QUAL'.format(mode),
#        'NIST.samtools.{}.VQSLOD'.format(mode),
        'NIST.UG.{}.VQSLOD'.format(mode),
#        'NIST.FB.{}.VQSLOD'.format(mode),
        'NIST.HC.{}.VQSLOD'.format(mode),
        'NIST.HC.0.{}.VQSLOD'.format(mode),
        'NIST.HC.1.{}.VQSLOD'.format(mode),
#        'NIST.HC.2.{}.VQSLOD'.format(mode),
        'NIST.HC.5.{}.VQSLOD'.format(mode),
#        'NIST.HC.10.{}.VQSLOD'.format(mode),
#        'NIST.HC.20.{}.VQSLOD'.format(mode),
        'NIST.HC.30.{}.VQSLOD'.format(mode),
        'NIST.HC.ERC_NONE.singlesample.{}.VQSLOD'.format(mode),
        )):

        ## VariantAnnotator failure
        if mode == 'indels' and 'VQSLOD' in affix and ('samtools' in affix or 'FB' in affix):
            continue

        print(affix)
        with open('out_vcfeval/{}/vcfeval.log'.format(affix)) as f:
            for line in f:
                match = re.match(r'.+baseline contains ([0-9]+) variants', line)
                if match:
                    baseline = true = int(match.group(1))
                    break

        x = [0]
        y = [0]
        with gzip.open('out_vcfeval/{}/weighted_roc.tsv.gz'.format(affix), 'rt') as f:
            for line in f:
                if line[0] == '#':
                    continue
                score, tp, fp = line.rstrip().split('\t')
                tpr = float(tp)/baseline
                fdr = int(fp)/(int(fp)+float(tp))
                x.append(int(fp))
                y.append(tpr)
#                print(affix, 'fp', fp, 'tp', tp, 'y=tpr', tpr, 'x=fdr', fdr)
#                y.append(tp)
        if 'VQSLOD' in affix:
            style = '-'
        else:
            style = '--'
        color = l_colors[i%len(l_colors)]
#zcat union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio_AddRTGPlatGenConf_filtNISTclustergt9_RemNISTfilt_RemPartComp_RemRep_RemPartComp_v0.2.bed.gz | awk '$1==20{sum+=$3-$2+1} END{print sum}'
#51798343
#FDR
        if type_plot == 'FDR':
            x = [_/(float(tp)+int(fp)) for _ in x]
#FPR
        if type_plot == 'FPR':
            x = [_/(51798343-baseline) for _ in x]
#        plt.plot(x, y, style, label=affix, color=color, linewidth=3)
        ax1.plot(x, y, style, label=affix, color=color, linewidth=3)

#    plt.xlabel('FDR = FP/(TP+FP)')
    if type_plot == 'FDR':
#        plt.xlabel('FDR = FP/(TP+FP)')
        ax1.set_xlabel('FDR = FP/(TP+FP)')
    if type_plot == 'FPR':
        ax1.set_xlabel('FPR = FP/N = FP/(FP+TN)')
#    plt.ylabel('TPR = TP/P = TP/(TP+FN)')
    ax1.set_ylabel('TPR = TP/P = TP/(TP+FN)')
    ax2.set_ylabel('TP')
#    plt.legend(loc='lower right', shadow=True)
    ax1.legend(loc='lower right', shadow=True)
    plt.suptitle(mode[:-1].upper()+'s', fontsize=14, fontweight='bold')
#    pyplt.set_ylim(0,1)
#    plt.ylim(0,1)

    axes = plt.gca()
    axes.set_ylim([0,1])
    ax1.set_ylim([0,1])
    axes.set_xlim([0,1])
    ax2.set_ylim([0,baseline])

#    plt.grid(b=True, which='major', color='k', linestyle='--')
    ax1.grid(b=True, which='major', color='k', linestyle='--')
    plt.savefig('{}_{}.png'.format(type_plot, mode), dpi=600)
#    plt.savefig('{}_{}.jpg'.format(type_plot, mode), dpi=600)
#    plt.savefig('{}_{}.eps'.format(type_plot, mode), dpi=600)
    plt.close()
    plt.clf()
