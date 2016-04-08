import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
#import statistics

pos, x, y, exp_freq_a1 = np.loadtxt(
    'comparison.1000Gp1_vs_1000Gp1_agv_ug2g.txt',
    unpack=True,
    dtype={
        'names': ('pos', 'x', 'y', 'exp_freq_a1'),
        'formats': ('<S30', np.float, np.float, np.float)},
    )

for exp_freq_a1_min, exp_freq_a1_max in ((0,0.01),(0.01,0.05),(0.05,0.95),(0.95,0.99),(0.99,1)):

    print(exp_freq_a1_min, exp_freq_a1_max)

    x2,y2 = list(zip(*((xi,yi) for xi,yi,exp_freq_a1i in zip(x,y,exp_freq_a1) if (exp_freq_a1i > exp_freq_a1_min and exp_freq_a1i <= exp_freq_a1_max))))
#    x2 = numpy.log10(x2)
#    y2 = numpy.log10(x2)
    print(len(x2), len(y2), np.mean(x2), np.mean(y2))

    heatmap, xedges, yedges = np.histogram2d(x2, y2, bins=100)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    extent = [0, 1, 0, 1]

    fig = plt.figure()
    plt.clf()
    ax = fig.add_subplot(111)
    ax.set_xlabel(
        '1000Gp1, ({}+/-{})'.format(
#            statistics.mean(x2), statistics.stdev(x2)))
            round(np.mean(x2), 2), round(np.std(x2), 2)))
    ax.set_ylabel(
        '1000Gp1+AGVP+UG2G, ({}+/-{})'.format(
#            statistics.mean(y2), statistics.stdev(y2)))
            round(np.mean(y2), 2), round(np.std(y2), 2)))

    fig.suptitle('r2_type0; {} < exp_freq_a1 <= {}'.format(exp_freq_a1_min, exp_freq_a1_max))

    plt.imshow(heatmap.T, extent=extent, origin='lower', norm=LogNorm())
#    plt.imshow(heatmap.T, extent=extent, origin='lower')
    plt.colorbar()
    plt.savefig('comparison.1000Gp1_vs_1000Gp1_agv_ug2g.{}.{}.png'.format(exp_freq_a1_min, exp_freq_a1_max), dpi=300)
