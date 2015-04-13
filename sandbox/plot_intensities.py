#!/bin/env python3

#Tommy Carstensen, Wellcome Trust Sanger Institute, January 2015

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import gzip
import argparse
import math
import os


def main():

    args = parse_arguments()
    print('looping intensities')
    with gzip.open(args.intensities, 'rt') as f:
        for line in f:
            l_samples_int = [s[:-1] for s in line.split()[3:-1:2]]
            break
        for i, line in enumerate(f):
            l = line.split()
            rsID = l[0]
#            print(rsID, args.rsID[0], len(args.rsID))
            if rsID in args.rsID:
                print('found', rsID, 'on line', i)
                if os.path.isfile(
                    os.path.join(args.out, '{}.png'.format(rsID))):
                    print('already plotted')
                    continue
                if args.prob:
                    prob = args.prob[args.rsID.index(rsID)]
                    del args.prob[args.rsID.index(rsID)]
                else:
                    prob = None
                plot_rsID(args, l, l_samples_int, prob)
                args.rsID.remove(rsID)
                if len(args.rsID) == 0:
                    break

    return


def plot_rsID(args, l, l_samples_int, prob):

    rsID = l[0]
    x, y = zip(*[[float(l[i]), float(l[i + 1])] for i in range(3, len(l), 2)])
    alleles = l[2]
    l_colors, chrom, pos = parse_bfiles(args, l_samples_int, alleles, rsID)
    if not l_colors:
        return
    fig, ax = plt.subplots()
    d_labels = {
        args.colors[0]: '{}{}'.format(alleles[0], alleles[0]),
        args.colors[1]: '{}{}'.format(alleles[0], alleles[1]),
        args.colors[2]: '{}{}'.format(alleles[1], alleles[1]),
        args.colors[3]: '00',
#        '#ffffff': 'Failed QC',
        }
    n = 0
    for color in args.colors:
        print(color)
        xj = []
        yj = []
        cj = []
        for xi, yi, ci in zip(x, y, l_colors):
            if ci != color:
                continue
            xj += [xi]
            yj += [yi]
            cj += [ci]
            n += 1
        ax.scatter(xj, yj, c=color, alpha=0.5, label=d_labels[color])
    lim_max = math.ceil(max(1, max(x), max(y)))
    ax.set_xlim((0, lim_max))
    ax.set_ylim((0, lim_max))
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax.set_aspect(abs(x1 - x0) / abs(y1 - y0))
    plt.xlabel('Intensity {}'.format(alleles[0]))
    plt.ylabel('Intensity {}'.format(alleles[1]))
    plt.legend(loc='upper right', shadow=True)
    title = '{}\n, chrom={}, pos={}, n={}'.format(rsID, chrom, pos, n)
    if prob:
        title += ', p={}'.format(prob)
    plt.title(title)
    plt.grid(b=True, which='major', linestyle='--')
    plt.savefig(os.path.join(args.out, '{}.png'.format(rsID)), dpi=150)
    plt.close()
    plt.clf()

    print('done', rsID)

    return


def parse_bfiles(args, l_samples_int, alleles, rsID):

    d_flip = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

    i, rsID_bim, chrom, pos, A1, A2 = parse_bim_line(args, rsID)
    if rsID_bim != rsID:
        return None, None, None

    with open('{}.fam'.format(args.bfile)) as fam:
        set_samples_fam = set([line.split()[0] for line in fam.readlines()])
        n_samples = len(set_samples_fam)
    n_bytes_per_SNP = math.ceil(n_samples / 4)

    # https://www.cog-genomics.org/plink2/formats
    magic_number = bytearray([108, 27])
    mode = bytearray([1])
    with open('{}.bed'.format(args.bfile), 'rb') as bed:
        bed.seek(len(magic_number) + len(mode) + i * n_bytes_per_SNP)
        bytesSNP = bed.read(n_bytes_per_SNP)
    int_bytesSNP2 = int.from_bytes(bytesSNP, 'little')
    l_colors = []
    i = -1
    if any([
        A1 == alleles[0] and A2 == alleles[1],
        A1 == d_flip[alleles[0]] and A2 == d_flip[alleles[1]], 
        (A1 == alleles[0] or A1 == d_flip[alleles[0]]) and A2 == '0']):
        bool_reverse = False
    else:
        bool_reverse = True
    for sample in l_samples_int:
        if not sample in set_samples_fam:
            l_colors += ['#ffffff']
            continue
        i += 1
        shift2 = 2 * i
        bits = int_bytesSNP2 >> shift2 & 0b11
        ## hom00
        if bits == 0:
            if not bool_reverse:
                l_colors += [args.colors[0]]
            else:
                l_colors += [args.colors[2]]
        ## het
        elif bits == 2:
            l_colors += [args.colors[1]]
        ## hom11
        elif bits == 3:
            if not bool_reverse:
                l_colors += [args.colors[2]]
            else:
                l_colors += [args.colors[0]]
        ## miss
        elif bits == 1:
            l_colors += [args.colors[3]]
        else:
            stop

    return l_colors, chrom, pos


def parse_bim_line(args, rsID):

    print('looping bim')
    with open('{}.bim'.format(args.bfile)) as bim:
        for i, line in enumerate(bim):
            if line.split()[1] == rsID:
                break

    chrom = line.split()[0]
    pos = line.split()[3]
    A1 = line.split()[4]
    A2 = line.split()[5]

    return i, line.split()[1], chrom, pos, A1, A2


def parse_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument('--intensities')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--position', type=int)
    group.add_argument('--rsID', nargs='+')
    parser.add_argument('--prob', nargs='*')
    parser.add_argument('--bfile')
    parser.add_argument('--out', required=True)
    parser.add_argument(
        '--colors', nargs='+',
        default=['#ff0000', '#00ff00', '#0000ff', '#000000'])
    args = parser.parse_args()

    assert len(args.rsID) == len(args.prob) or len(args.prob) == 0

    return args

if __name__ == '__main__':
    main()
