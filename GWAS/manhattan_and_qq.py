#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, July 2013, June 2017

import argparse
import os
import urllib.request
import sys
import fileinput
import gzip
import pyensembl
import sys
import fileinput
import itertools
import operator
import math
import pysam
import matplotlib.pyplot as plt
import statsmodels.api as sm
import numpy as np
##from scipy.stats import chi2
import scipy
from pyliftover import LiftOver


def main():

    args = argparser()

    os.environ['PYENSEMBL_CACHE_DIR'] = args.pyensembl
    ensembl_object = pyensembl.EnsemblRelease(release=75)

#    download_ensembl(args, ensembl_object)

    ## If dbSNP argument is not a file, then it's a column of rsIDs.
    if not os.path.isfile(args.dbSNP):
        args.dbSNP = int(args.dbSNP)

    ## http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=4
    cycle_colors = itertools.cycle(
        (('#a6cee3', '#1f78b4'), ('#b2df8a', '#33a02c')))

    pos_init_chrom = 0
    pos_prev = 0
    fi = fileinput.FileInput(
        files=args.input, openhook=fileinput.hook_compressed)
    l_x = []
    l_y = []
    l_c = []
    l_prob = []
    x_ticks = []
    ## Gene annotations a local maxima.
    annotations = []
    y_max = 0
    d_pos_init_chrom = {}
    for chrom, split_lines in itertools.groupby(
        split_line(fi), operator.itemgetter(args.chrom)):
        print('Looping over chromosome', chrom, file=sys.stderr)
        pos_init_chrom += pos_prev + 1500000  # todo: make argument
        d_pos_init_chrom[chrom] = pos_init_chrom
        colors = next(cycle_colors)
        x_ticks.append((None, None))
        pos = 0
        for l in split_lines:
            af = float(l[args.af])
            if min(af, 1 - af) < args.threshold_maf:
                continue
            pos, pos_prev = int(l[args.pos]), pos
            ## Assert that input is sorted.
            assert pos >= pos_prev, (pos, pos_prev)
            prob = float(l[args.prob])
            ref = l[args.ref]
            alt = l[args.alt]
            x = pos + pos_init_chrom
            try:
                y = -math.log10(prob)
            except:
                ## GEMMA seems to print a p_lrt of 0.000000e+00, when the numbers get too small.
                ## e.g. bilirubin 2:234664586
                assert l[args.prob] == '0.000000e+00'
                y = y+0.01
            y_max = max(y, y_max)
            l_x.append(x)
            l_y.append(y)
            l_prob.append(prob)
            ## Place a chromosome tick halfway through the chromosome basepair range.
            x_ticks[-1] = (pos_init_chrom + pos / 2, chrom)
            if prob > args.threshold_p:
                l_c.append(colors[0])
            else:
                l_c.append(colors[1])
                gene_names = ensembl_object.gene_names_at_locus(chrom, pos)
                gene_ids = ensembl_object.gene_ids_at_locus(chrom, pos)
                protein_ids = ensembl_object.protein_ids_at_locus(chrom, pos)
                transcript_ids = ensembl_object.transcript_ids_at_locus(
                    chrom, pos)
##                print('protein_ids', protein_ids)
##                print('gene_names', gene_names)
##                print('gene_ids', gene_ids)
##                print('transcript_ids', transcript_ids)
####                for gene_id in gene_ids:
####                    locus = ensembl_object.locus_of_gene_id(gene_id)
######                    print(gene_id, locus.start, locus.end, locus.strand)
                if os.path.isfile(args.dbSNP):
                    rsID = parse_dbSNP(args, chrom, pos, ref, alt)
                else:
                    rsID = l[args.dbSNP]
                annotation = {
                    'chrom': chrom,
                    'x': x,
                    'y': y,
                    'prob': prob,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'rsID': rsID,
                    'gene_ids': gene_ids,
                    'gene_names': gene_names,
                    'af': af,
                    }
##                print(
##                    prob, af, chrom, pos, rsID, ref, alt,
##                    ','.join(gene_ids), ','.join(gene_names),
##                    sep='\t', file=sys.stdout)

                ## Don't append to the cluster, if it doesn't have an rsID.
                ## Pick something with a lower probability instead then.
                if not rsID:
                    pass
                ## No clusters yet.
                ## Create first cluster.
                elif not annotations:
                    annotations.append(annotation)
                ## Not in the vicinity of previous cluster.
                ## Append new cluster.
                elif (
                    chrom != annotations[-1]['chrom'] or
                    pos - annotations[-1]['pos'] > 1000000):  # todo: make arg!
                    annotations.append(annotation)
                ## In the vicinity of previous cluster.
                ## Probability lower than current local minimum.
                ## Overwrite previous cluster.
                elif prob < annotations[-1]['prob']:
                    annotations[-1] = annotation

    print('annotations', annotations, file=sys.stderr)

    plot_qq(args, l_y, l_prob)
    plt.clf()
    plot_manhattan(
        args, annotations, l_x, l_y, l_c, x_ticks, y_max, d_pos_init_chrom)

    return


def plot_qq(args, l_y, l_prob):

    e = -np.log10(ppoints(len(l_y)))
#    o = -np.log10(sorted(l_prob))
    o = sorted(l_y, reverse=True)
    print('plt.scatter(qq)', file=sys.stderr)
    plt.scatter(e, o, s=2)
    inflation = get_lambda(l_prob)
    with open('{}.qq.txt'.format(args.out), 'w') as f:
        f.write(str(inflation)+'\n')
    plt.text(
        0.76, 0.02, '$\lambda$ = {:.2f}\nn = {:,}'.format(inflation, len(l_y)), transform=plt.gca().transAxes, fontsize=7,
        verticalalignment='bottom',
        bbox=dict(boxstyle='round', facecolor='blue', alpha=0.5),
        )
    plt.title(args.title, fontsize='small')
    plt.ylabel(r'Observed -log$_{10}$($p$)')
    plt.xlabel(r'Expected -log$_{10}$($p$)')
    plt.xlim(0, math.ceil(max(args.min_y, o[0], e[0])))
    plt.ylim(0, math.ceil(max(args.min_y, o[0], e[0])))
    plt.gca().set_aspect('equal')
    plt.plot([0, 1], [0, 1], transform=plt.gca().transAxes, ls='--', c='k')
    print('plt.savefig( {}.qq.png )'.format(args.out), file=sys.stderr)
    plt.savefig('{}.qq.png'.format(args.out))

##    qqplot(args, [np.array(l_prob)], [args.out])
##    sys.exit()

##    import numpy as np
##    expected = reversed(map(operator.neg, map(math.log10, ppoints(len(l_y), 0.5))))
##    expected = reversed(map(operator.neg, map(math.log10, )))
##    sm.qqplot(np.array(l_prob), line='45')  # wrong!!!

##    sm.qqplot(np.array(l_y), line='45')
##    observed = np.array(list(
##        map(operator.neg, map(math.log10, sorted(l_prob, reverse=True)))))
##    sm.qqplot(observed, line='45')
##    print(list(observed)[:10])
##    print(type(observed))
##    print(list(expected)[:10])
##    print(type(expected))

    return


def get_lambda(p):
    pm = np.median(p)
    chi = scipy.stats.chi2.ppf(1. - pm, 1)
    return chi / scipy.stats.chi2.ppf(0.5, 1)


def ppoints(n, a=0.5):
    return (np.arange(n) + 1 - a)/(n + 1 - 2*a)


def plot_manhattan(
    args, annotations, l_x, l_y, l_c, x_ticks, y_max, d_pos_init_chrom):

    y_max = max(int(y_max + 3), args.min_y)

    if args.EFO:
        ## Just make some assumptions about builds here for now.
        ## https://en.wikipedia.org/wiki/Reference_genome
        lo = LiftOver('hg38', 'hg19')
        with open(args.EFO) as f:
            cnt = collections.Counter()
            for line in f:
                cnt[line.split('\t')[7]] += 1
            trait_most_common = cnt.most_common(1)[0][0]
        with open(args.EFO) as f:
            ## Skip header.
            for line in f:
                break
            for line in f:
                l = line.split('\t')
#                ## Try to weed out all the garbage present in the GWAS catalog.
#                if not l[7] == trait_most_common:
#                    continue
                CHR_ID = l[11]
                ## Skip if missing data.
                if CHR_ID == '':
                    continue
                try:
                    CHR_POS = int(l[12])
                ## Continue if CHR_POS is not an integer.
                except ValueError:
                    continue
                rsID = l[21]
                y = PVALUE_MLOG = min(y_max, float(l[28]))
#                if y < -math.log10(args.threshold_p):
#                    continue
                try:
                    x = d_pos_init_chrom[CHR_ID] + lo.convert_coordinate(
                        'chr{}'.format(CHR_ID), CHR_POS)[0][1]
                except KeyError:
                    assert CHR_ID == 'X'
                    continue
                except IndexError:
                    print('IndexError', CHR_ID, CHR_POS, lo.convert_coordinate('chr{}'.format(CHR_ID), CHR_POS), file=sys.stderr)
                    continue
#                l_x.append(x)
#                l_y.append(y)
#                l_c.append('#FF0000')
                ## Colour most frequently occuring trait red.
                if l[7] == trait_most_common:
                    plt.vlines(x, 0, y, colors='#FF0000', linewidth=0.5, linestyle='--')
                ## Colour less frequently occuring traits orange,
                ## because these might be junk in the GWAS catalog.
                else:
                    plt.vlines(x, 0, y, colors='#FF8000', linewidth=0.5, linestyle='--')

    n = len(l_y)

    plt.ylabel(r'-log$_{10}$($p$)')

#    plt.axhline(-math.log10(0.05 / n), color='0.8', linewidth=0.5)
#    plt.axhline(-math.log10(5 * 10 ** -8), color='0.5', linewidth=0.5)
    plt.axhline(-math.log10(args.threshold_p), color='0.2', linewidth=0.5, linestyle='--')
    try:
        plt.ylim((0, y_max))  # todo: make argument
    except:
        pass

    print('plt.scatter(manhattan)', file=sys.stderr)
    plt.scatter(l_x, l_y, c=l_c, s=3)

    plt.title(args.title, fontsize='small')

    for annotation in annotations:
#        if annotation['prob'] > 0.05 / n:
        if annotation['prob'] > args.threshold_p:
            continue
        print('\t'.join(
            [str(annotation[k]) for k in sorted(annotation.keys())]))
        plt.annotate(
            '\n'.join((
                'p={:.1E}'.format(annotation['prob']),
                'pos={:,}'.format(annotation['pos']),
                'MAF={:.3f}'.format(min(annotation['af'], 1 - annotation['af'])),
                annotation['rsID'],
                ','.join(annotation['gene_names']),
                )),
            xy=(annotation['x'], annotation['y']),
##            xytext=(),
            fontsize='xx-small',
            horizontalalignment='center',
            verticalalignment='bottom',
            rotation=30,
            )

    plt.xticks(
        *zip(*x_ticks),
        rotation=-75, size=6, fontsize=6)

    print('plt.savefig( {}.manhattan.png )'.format(args.out), file=sys.stderr)
    plt.savefig('{}.manhattan.png'.format(args.out), dpi=600)

    return


def download_ensembl(args, ensembl_object):

    if not os.path.isdir(
        os.path.join(args.pyensembl, 'pyensembl/GRCh37/ensembl75')):
        print(
            'Downloading Ensembl release, which might take a few minutes.',
            file=sys.stderr)
        ensembl_object.download()
        ensembl_object.index()

    return


def split_line(fi):
    for l in map(str.split, fi):
        yield l


def parse_dbSNP(args, chrom, pos, ref, alt):

    ## todo: check this function... I might miss variants in dbSNP...

    tbx = pysam.TabixFile(args.dbSNP)
    row = None
    for row in tbx.fetch(chrom, pos - 1 - 1, pos, parser=pysam.asTuple()):
        print(row, file=sys.stderr)
        if any([
            ## SNP.
            row[3] == ref and alt in row[4].split(','),
            ## SNP in dbSNP and MNP in preselected.txt
            all([
                row[3] == ref[0],
                len(set(alt[0].split(',')) & set(row[4].split(','))) > 0,
                len(ref) in map(len, alt.split(',')),
                len(row[3]) == 1,
                ]),
            ## Insertion.
            all([
                int(row[1]) == pos,
                len(row[3]) == 1,
                ref == '-',
                ## One or more ALTs overlap (e.g. rs3835252).
                len(set(x[1:] for x in row[4].split(',')) & set(alt.split(','))) >= 1,
                ]),
            ## Deletion.
            all([
                int(row[1]) == pos,
                len(row[4]) == 1,
                alt == '-',
                row[3][:1] == ref,
                len(row[3]) > 1,
                len(row[4]) == 1,
                ]),
            ## Deletion.
            all([
                int(row[1]) + 1 == pos,
                len(row[3]) == len(ref) + 1,
                set(map(len, row[4].split(','))) == set([1]),
                alt == '-',
                row[3][1:] == ref,
                ]),
            ]):
            rsID = row[2]
            break
    ## Not found in dbSNP.
    else:
        rsID = ''

    return rsID


def argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument('--title', required=True)
    parser.add_argument('--input', required=True, nargs='+')
    parser.add_argument('--threshold_p', default=1 * 10 ** -9, type=float)
    parser.add_argument('--threshold_maf', default=0, type=float)
    parser.add_argument('--min_y', default=50, type=int)
    parser.add_argument(
        '--dbSNP', '--rsID', required=True,
        help='Path to dbSNP VCF or column of rsIDs')
    parser.add_argument(
        '--pyensembl', default=os.getcwd(), help='Cache location of pyensembl')
    parser.add_argument('--out', required=False)
    ## Columns.
    parser.add_argument('--chrom', default=0, type=int)
    parser.add_argument('--pos', default=1, type=int)
    parser.add_argument('--prob', default=2, type=int)
    parser.add_argument('--ref', type=int)
    parser.add_argument('--alt', type=int)
    parser.add_argument('--af', type=int)
    parser.add_argument('--EFO')

    args = parser.parse_args()

    if not args.out:
        args.out = args.title
    if args.threshold_maf:
        args.title += ', MAF > {}'.format(args.threshold_maf)

    return args


if __name__ == '__main__':
    main()
