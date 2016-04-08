import matplotlib.pyplot as plt
from numpy.random import normal
import gzip

with gzip.open('merged.ug2g1977.uganda_gwas4435.LDprunedMAF05.genome.gz', 'rt') as f:
    for line in f:
        break
    x = []
    for line in f:
        x.append(float(line.split()[9]))
    print(len(x))
    plt.hist(x, bins=100, range=[0.0, 1.0], log=True)
#    plt.title("IBD")
    plt.xlabel("IBD")
#    plt.ylabel("Frequency (%)")
    plt.ylabel('Count')
    plt.savefig('merged.ug2g1977.uganda_gwas4435.LDprunedMAF05.genome.png')
