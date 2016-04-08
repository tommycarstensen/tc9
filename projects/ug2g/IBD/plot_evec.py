import numpy as np
import matplotlib.pyplot as plt
import collections

with open('../../../metadata/pop.dic') as f:
    d_sample2pop = {line.split()[0]:line.split()[1] for line in f}

file_evec = 'out_EIGENSOFT/merged.ug2g1986.uganda_gwas4435.LDprunedMAF05.evec'
#file_evec = 'out_EIGENSOFT/merged.ug2g343.uganda_gwas343.LDprunedMAF05.evec'
with open(file_evec, 'r') as f:
    f.readline()
    l = []
    for line in f.readlines():
        l.append(line.rstrip().split())

with open('../../../metadata/APP2EGAN.tsv') as f:
    d_EGAN2APP = {line.split()[1]:line.split()[0] for line in f}

l = list(zip(*l))
samples = [s.split(':')[0] for s in l[0]]
d_sample2color = {}
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
'#000000',
]
pops = []
cnt = collections.Counter()
for sample in samples:
    try:
        sample = d_EGAN2APP[sample]
        pop = 'Illumina HiSeq 2000/2500'
    except KeyError:
        pop = 'Illumina Omni2.5-8'
        pass
#    pop = d_sample2pop[sample]
    pops.append(pop)
    cnt[pop] += 1

d_pop2color = {pop: l_colors[i] for i, (pop, v) in enumerate(cnt.most_common())}
colors = [d_pop2color[pop] for pop in pops]

## Shrink current axis's height by 10% on the bottom
#box = plt.get_position()
#plt.set_position([box.x0, box.y0 + box.height * 0.1,
#                 box.width, box.height * 0.9])

for PC1, PC2 in ((1,2),(3,4),(5,6),(7,8),(9,10)):
    print(PC1, PC2)
    x = l[PC1]
    y = l[PC2]
    plt.xlabel('PC {}'.format(PC1))
    plt.ylabel('PC {}'.format(PC2))
    for i, (popi, _) in enumerate(cnt.most_common()):
        xi = []
        yi = []
        for j, popj in enumerate(pops):
            if popi != popj:
                continue
            xi.append(x[j])
            yi.append(y[j])
        plt.plot(xi, yi, 'o', label=popi, alpha=0.9, color=l_colors[i])
#    plt.scatter(x, y, c=colors, alpha=0.9)
#    plt.legend()

#    ## above
#    plt.legend(bbox_to_anchor=(0, 1.05, 1, .105), loc=3, ncol=2, borderaxespad=0., shadow=True, mode='expand', numpoints=1, prop={'size':7})
    ## right

#    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., shadow=True, numpoints=1, prop={'size':7})

# Put a legend below current axis
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.10),
          fancybox=True, shadow=True, ncol=2, numpoints=1)

    plt.savefig('{}_PCs{}_{}.png'.format(file_evec.replace('.','_'), PC1, PC2), dpi=300)
    plt.close()
