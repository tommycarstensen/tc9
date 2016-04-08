import argparse
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import datetime
from matplotlib.ticker import ScalarFormatter
from scipy.interpolate import spline
import numpy as np


def main():

    args = parse_args()

    fig1, ax1 = plt.subplots()

    n_samples_max = 0
    linestyles = ['-', '--', '-.', ':']
    for path, color, linestyle, linewidth, label, path_extract in zip(
        args.info, args.colors, args.linestyles, args.linewidths, args.labels,
        args.extract):
        if path_extract:
            with open(path_extract) as f:
                set_extract = set(f.read().rstrip().split())
        d_r2 = {'cnt': {}, 'sum': {}, 'sumsq': {}}
        for dirpath, dirnames, filenames in os.walk(path):
            for filename in filenames:
                if not filename.endswith('info'):
                    continue
                filepath = os.path.join(dirpath, filename)
                print(filepath)
                ## Get sample count.
                with open(filepath[:-5]+'_samples') as f_samples:
                    for n_samples, line in enumerate(f_samples):
                        pass
                n_samples -= 2
                n_samples_max = max(n_samples_max, n_samples)
                with open(filepath) as f:
                    line = f.readline()
                    l = line.rstrip().split()
                    i_freq = l.index('exp_freq_a1')
                    i_r2 = l.index('r2_type0')
                    for line in f:
                        l = line.split()

                        r2 = float(l[i_r2])
                        if r2 == -1:
                            continue

                        if path_extract:
                            if not ':'.join((l[2],l[3],l[4])) in set_extract:
#                                if l[8] != '3':
#                                    print('skipping', l[0], l[1], ':'.join((l[2],l[3],l[4])), list(set_extract)[:2])
                                continue

                        freq = float(l[i_freq])
                        freq = round(round(freq * n_samples * 2, 0) / (n_samples * 2), 3)
#                        print(r2, freq)
                        try:
                            d_r2['sum'][freq] += r2
                            d_r2['cnt'][freq] += 1
                            d_r2['sumsq'][freq] += (r2 * r2)
                        except KeyError:
                            d_r2['sum'][freq] = r2
                            d_r2['cnt'][freq] = 1
                            d_r2['sumsq'][freq] = r2 * r2

        x = []
        y = []
        for freq in sorted(d_r2['cnt'].keys()):
            x.append(freq)
            y.append(d_r2['sum'][freq]/d_r2['cnt'][freq])
#           stddev = math.sqrt( ((sumxx-(sumx**2)/n))/n )
#        x_smooth = np.linspace(min(x), max(x), 1000)
#        y_smooth = spline(x, y, x_smooth)
        ax1.plot(
            x, y,
#            x_smooth, y_smooth,
            label=label, linestyle=linestyles[linestyle], color=color,
            linewidth=linewidth,
            )

    ## Set xticks.
    t = (1, 0.5, 0.2)
    div = 1
    xticks = []
    while True:
        for i in t:
            if i / div < 1 / (n_samples_max * 2):
                break
            else:
                xticks.append(i / div)
                continue
            break
        if i / div < 1 / (n_samples_max * 2):
            break
        div *= 10
    ax1.set_xticks(xticks)

    ax1.set_xscale('log')
#    ax1.set_xticks([0, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.00])
    ax1.set_xlim(0, 1)
#    ax1.legend(loc='lower right')
#    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#    ax1.set_xlabel('exp_freq_a1')
#    ax1.set_ylabel('r2_type0')
    ax1.set_xlabel('Allele frequency')
    ax1.set_ylabel('Squared correlation between input and masked/imputed genotypes')
    ax1.get_xaxis().set_major_formatter(ScalarFormatter())
    ax1.grid('on')

    handles, labels = ax1.get_legend_handles_labels()
    lgd = ax1.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1))
    lgd = ax1.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(
        '{}.png'.format(args.affix),
        bbox_extra_artists=(lgd,), bbox_inches='tight')

    print('xtics', xticks)
    return


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('--info', nargs='+')
    parser.add_argument('--out', help='output affix')
    parser.add_argument('--title')
    parser.add_argument('--sample_count')
    parser.add_argument('--affix')
    parser.add_argument('--colors', nargs='*')
    parser.add_argument('--linestyles', nargs='*', type=int)
    parser.add_argument('--linewidths', nargs='*', type=float)
    parser.add_argument('--labels', nargs='*', type=str)
    parser.add_argument('--extract', nargs='*')

    args = parser.parse_args()

    if not args.colors:
        args.colors = len(args.info)*[None]
    if not args.linestyles:
        args.linestyles = len(args.info)*[0]
    if not args.linewidths:
        args.linewidths = len(args.info)*[None]
    if not args.labels:
        args.labels = args.info
    if not args.extract:
        args.extract = len(args.info)*[None]
    assert len(args.colors) == len(args.info)
    assert len(args.linestyles) == len(args.info)
    assert len(args.linewidths) == len(args.info)
    assert len(args.labels) == len(args.info)

    if not args.affix:
        args.affix = datetime.date.today().isoformat()

    return args


if __name__ == '__main__':
    main()
