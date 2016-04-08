import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['backend'] = 'TkAgg' 
#plt.rcParams['backend'] = 'Agg' 

with open('het.dat') as fd:
    x = [float(het) for het in list(zip(*[line.split() for line in fd]))[1]]

#plt.plot(x, y, style, label=affix, color=color, linewidth=3)
#        print(len(x), len(y), affix, x[0], y[0], x[-1], y[-1])

#    plt.xlabel('FDR = FP/(TP+FP)')

plt.xlabel('Heterozygosity')
plt.ylabel('Count')
#plt.suptitle(mode[:-1].upper()+'s', fontsize=14, fontweight='bold')
#axes = plt.gca()
#axes.set_ylim([0,1])
#plt.grid(b=True, which='major', color='k', linestyle='--')
plt.hist(x, bins=200)
plt.savefig('het.png', dpi=600)
plt.close()
plt.clf()
