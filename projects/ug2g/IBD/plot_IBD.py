import numpy as np
import matplotlib.pyplot as plt
#import collections
import gzip

affix='merged.ug2g1977.uganda_gwas4435.LDprunedMAF05'
d_maxIBD = {}
d_sumIBD = {}
d_cntIBD = {}
with gzip.open('{}.genome.gz'.format(affix), 'rt') as f:
    seqseq = []
    seqchip = []
    chipchip = []
    ## skip header
    f.readline()
    for line in f:
        l = line.rstrip().split()
        IBD = float(l[9])
        if l[0][0] == 'E' and l[2][0] == 'E':
            seqseq.append(IBD)
        elif l[0][0] == 'E' or l[2][0] == 'E':
            seqchip.append(IBD)
        else:
            chipchip.append(IBD)
        for sampleID in (l[0], l[2],):
            try:
                maxIBD = d_maxIBD[sampleID]
            except KeyError:
                maxIBD = 0
                d_sumIBD[sampleID] = 0
                d_cntIBD[sampleID] = 0
            if IBD > maxIBD:
                d_maxIBD[sampleID] = IBD
            d_sumIBD[sampleID] += IBD
            d_cntIBD[sampleID] += 1

l_maxIBD = list(d_max_IBD.values())
l_avgIBD = [d_sumIBD[sampleID]/d_cntIBD[sampleID] for sampleID in d_sumIBD.keys()]

plt.xlabel('IBD maximum')
plt.ylabel('Count')
plt.hist(l_maxIBD, bins=50, range=(0,1))
#plt.legend(loc='upper right', shadow=True, numpoints=1)
plt.savefig('{}.IBDmax.png'.format(affix), dpi=300)
#plt.hist([seqseq,chipchip,seqchip], label=['both sequence', 'both chip', 'one sequence, one chip'], bins=50, range=(0.1,1), normed=True)
#plt.savefig('{}.upper.IBD.png'.format(affix), dpi=300)
plt.close()


plt.xlabel('IBD average')
plt.ylabel('Count')
plt.hist(l_avgIBD, bins=50, range=(0,1))
#plt.legend(loc='upper right', shadow=True, numpoints=1)
plt.savefig('{}.IBDav.png'.format(affix), dpi=300)
#plt.hist([seqseq,chipchip,seqchip], label=['both sequence', 'both chip', 'one sequence, one chip'], bins=50, range=(0.1,1), normed=True)
#plt.savefig('{}.upper.IBD.png'.format(affix), dpi=300)
plt.close()

stop

plt.xlabel('IBD')
plt.ylabel('Count')
plt.hist([seqseq,chipchip,seqchip], label=['both sequence', 'both chip', 'one sequence, one chip'], bins=50, range=(0,0.05))
plt.legend(loc='upper right', shadow=True, numpoints=1)
plt.savefig('{}.IBD.png'.format(affix), dpi=300)
#plt.hist([seqseq,chipchip,seqchip], label=['both sequence', 'both chip', 'one sequence, one chip'], bins=50, range=(0.1,1), normed=True)
#plt.savefig('{}.upper.IBD.png'.format(affix), dpi=300)
plt.close()
