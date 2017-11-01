#!/bin/env python3

#Tommy Carstensen, Wellcome Trust Sanger Institute, January 2015

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import gzip
import argparse
import math
import os
import operator


def main():

    args = parse_arguments()

    print('looping intensities')
    if args.fileformat == 'flat':
        flat(args)
    else:
        matrix(args)

    return


def flat(args):

    l = []
#    with open(args.bfile+'.bim') as f:
#        for i, line in enumerate(f):
#            l.append((operator.itemgetter(1, 3)(line.rstrip().split()), i))
#    for j
#            stop
#    assert len(args.rsID) == 1, args.rsID
    d_xy = {rsID: [[], []] for rsID in args.rsID}
    l_samples_int = []
    with open_file(args.intensities) as f:
        for line in f:
            l = line.rstrip().split()
            rsID = l[0]
            if rsID in args.rsID:
#                ## Append alleles.
#                d_alleles[rsID] = (l[2], l[3])
                ## Append intensities.
                d_xy[rsID][0].append(float(l[7]))
                d_xy[rsID][1].append(float(l[8]))
                ## Append sample ID.
                if len(d_xy[rsID][0]) > len(l_samples_int):
                    l_samples_int.append(l[1])

    ## Append alleles.
    d_alleles = {}
    with open(args.bfile+'.bim') as f:
        for line in f:
            l = line.rstrip().split()
            rsID = l[1]
            if rsID in args.rsID:
                d_alleles[rsID] = (l[4], l[5])

    prob = None
    for rsID in args.rsID:
        assert d_xy[rsID] != [[],[]]
        alleles = d_alleles[rsID]
        x = d_xy[rsID][0]
        y = d_xy[rsID][1]
        plot_rsID(args, rsID, x, y, alleles, l_samples_int, prob)

    return


def matrix(args):

    with open_file(args.intensities) as f:
        for line in f:
            # uganda_gwas
            l_samples_int = [s[:-1] for s in line.split()[3:-1:2]]
#            # vaccgene
#            l_samples_int = [s[:-2] for s in line.split()[3:-1:2]]
            break
        for i, line in enumerate(f):
            l = line.split()
            rsID = l[0]
#            print(rsID, l[1], l[2], l[3])
#            print(rsID, args.rsID[0], len(args.rsID))
            if rsID in args.rsID:
                print('found', rsID, 'on line', i)
                if os.path.isfile(
                    os.path.join('{}.png'.format(rsID))):
                    print('already plotted')
                    pass
                else:
                    if args.prob:
                        prob = args.prob[args.rsID.index(rsID)]
                        del args.prob[args.rsID.index(rsID)]
                    else:
                        prob = None
                    rsID = l[0]
                    x, y = zip(*[[float(l[i]), float(l[i + 1])] for i in range(3, len(l), 2)])
                    alleles = l[2]
                    plot_rsID(args, rsID, x, y, alleles, l_samples_int, prob)
                args.rsID.remove(rsID)
                if len(args.rsID) == 0:
                    break

    return


def open_file(path):

    if os.path.splitext(path)[1] == '.gz':
        f = gzip.open(path, 'rt')
    else:
        f = open(path)

    return f


def plot_rsID(args, rsID, x, y, alleles, l_samples_int, prob):

#    l_colors = len(x)*['#000000']
#    chrom = 4
#    pos = 145023404
#    MAF = 'xx'
    l_colors, chrom, pos, MAF = parse_bfiles(
        args, l_samples_int, alleles, rsID)
    if not l_colors:
        return
    fig, ax = plt.subplots()
    plt.rcParams.update({'axes.titlesize': 'small'})
    d_labels = {
        args.colors[0]: '{}{}'.format(alleles[0], alleles[0]),
        args.colors[1]: '{}{}'.format(alleles[0], alleles[1]),
        args.colors[2]: '{}{}'.format(alleles[1], alleles[1]),
        args.colors[3]: '00',
#        '#ffffff': 'Failed QC',
        }
    n = 0
    print(x, y, l_colors)
    assert len(x) == len(l_colors)
    for color in args.colors:
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
    print(len(x), len(y), len(l_colors))
    while '#ffffff' in l_colors:
        l_colors.remove('#ffffff')
    print(len(x), len(y), len(l_colors))
    lim_max = math.ceil(max(1, max(x), max(y)))
    ax.set_xlim((0, lim_max))
    ax.set_ylim((0, lim_max))
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax.set_aspect(abs(x1 - x0) / abs(y1 - y0))
    plt.xlabel('Intensity {}'.format(alleles[0]))
    plt.ylabel('Intensity {}'.format(alleles[1]))
    plt.legend(loc='upper right', shadow=True)
    title = '{}\nchrom={}, pos={}, n={}, MAF={:.3f}'.format(
        rsID, chrom, pos, n, MAF)
    if prob:
        title += ', p={}'.format(prob)
    if args.text:
        title += '\n{}'.format(args.text.replace('\t',' '))
    plt.title(title)
    plt.grid(b=True, which='major', linestyle='--')
    plt.savefig(os.path.join('{}.png'.format(rsID)), dpi=150)
    plt.close()
    plt.clf()

    print('done', rsID)

    return


def parse_bfiles(args, l_samples_int, alleles, rsID):

    d_flip = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

    i_SNP, rsID_bim, chrom, pos, A1, A2 = parse_bim_line(args, rsID)
    if rsID_bim != rsID:
        print(rsID_bim, rsID)
        return None, None, None

    with open('{}.fam'.format(args.bfile)) as fam:
        set_samples_fam = set([line.split()[0] for line in fam.readlines()])
    with open('{}.fam'.format(args.bfile)) as fam:
        l_samples_fam = [line.split()[0] for line in fam.readlines()]
    n_samples = len(set_samples_fam)
    if args.keep:
        with open('{}'.format(args.keep)) as fam:
            set_samples_fam -= set([line.split()[0][9:] for line in fam.readlines()])
    n_bytes_per_SNP = math.ceil(n_samples / 4)

    # https://www.cog-genomics.org/plink2/formats
    magic_number = bytearray([108, 27])
    mode = bytearray([1])
    with open('{}.bed'.format(args.bfile), 'rb') as bed:
        bed.seek(len(magic_number) + len(mode) + i_SNP * n_bytes_per_SNP)
        bytesSNP = bed.read(n_bytes_per_SNP)
    int_bytesSNP2 = int.from_bytes(bytesSNP, byteorder='little', signed=False)
    if any([
        A1 == alleles[0] and A2 == alleles[1],
        A1 == d_flip[alleles[0]] and A2 == d_flip[alleles[1]], 
        (A1 == alleles[0] or A1 == d_flip[alleles[0]]) and A2 == '0']):
        bool_reverse = False
    else:
        bool_reverse = True

    l_colors, MAF = parse_genotypes(
        args, bool_reverse, l_samples_int, set_samples_fam, l_samples_fam,
        int_bytesSNP2)

    return l_colors, chrom, pos, MAF


def parse_genotypes(
    args, bool_reverse, l_samples_int, set_samples_fam, l_samples_fam,
    int_bytesSNP2):

    l_colors = []

    ## Number of allele observations
    cnt1 = 0
    cnt2 = 0
    assert len(set_samples_fam & set(l_samples_int)) > 0
    for sample in l_samples_int:
        if not sample in set_samples_fam:
##            l_colors += ['#ffffff']
            l_colors += [None]
            continue
        i = l_samples_fam.index(sample)
        shift2 = 2 * i
        bits = int_bytesSNP2 >> shift2 & 0b11
        ## hom00
        if bits == 0:
            if not bool_reverse:
                l_colors += [args.colors[0]]
            else:
                l_colors += [args.colors[2]]
            cnt1 += 2
            cnt2 += 2
        ## het
        elif bits == 2:
            l_colors += [args.colors[1]]
            cnt1 += 2
            cnt2 += 1
        ## hom11
        elif bits == 3:
            if not bool_reverse:
                l_colors += [args.colors[2]]
            else:
                l_colors += [args.colors[0]]
            cnt1 += 2
        ## miss
        elif bits == 1:
            l_colors += [args.colors[3]]
            print('missing', sample)
        else:
            stop

    print(cnt1, cnt2)
    MAF = min(cnt2/cnt1, 1-cnt2/cnt1)

    return l_colors, MAF


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
    parser.add_argument('--prob', nargs='*', required=False, help='String to be printed on figure')
    parser.add_argument('--bfile')
    parser.add_argument('--keep')
#    parser.add_argument('--out', required=True)
    parser.add_argument(
        '--colors', nargs='+',
        default=['#ff0000', '#00ff00', '#0000ff', '#000000'])
    parser.add_argument('--text', '--title', required=False)
    parser.add_argument('--fileformat', choices=('flat', 'matrix'), default='flat')

    args = parser.parse_args()

    if args.prob:
        assert len(args.rsID) == len(args.prob)

    assert os.path.isfile(args.bfile+'.bed')
    assert os.path.isfile(args.bfile+'.bim')
    assert os.path.isfile(args.bfile+'.fam')

    return args


if __name__ == '__main__':
    main()
