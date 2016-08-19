#!/bin/env python3

# Tommy Carstensen
# Wellcome Trust Sanger Institute
# 2016Jun, 2016Aug19

import argparse
from Bio.bgzf import BgzfWriter
import gzip
import sys
import os
import math


def main():

    '''This script assumes identical ordering of samples across VCF files.'''

    args = parse_args()

    ## make out dir
    if not os.path.isdir(os.path.dirname(args.out)):
        os.mkdir(os.path.dirname(args.out))

    ## parse samples to be removed
    if args.remove:
        with open(args.remove) as f:
            set_remove = set([line.rstrip() for line in f])
    else:
        set_remove = set()

    if args.conversion == 'PL2GL':
        format1 = 'PL'
        func = PL2GL
    elif args.conversion == 'GP2GL':
        format1 = 'GP'
        func = GP2GL

#    with BgzfWriter(args.out, 'wb') as fd_out:
    with gzip.open(args.out, 'wt') as fd_out:
        for i_vcf, vcf in enumerate(args.vcf):
            with gzip.open(vcf, 'rt') as fd_vcf:

                ## loop metadata lines and header line
                for line in fd_vcf:
                    ## end of metadata lines
                    if line[0:2] != '##':
                        ## break loop at the header line
                        break
                    else:
                        ## Print metadata lines of the first VCF file.
                        if i_vcf == 0:
                            print(line, end='', file=fd_out)

                n = n_samples = len(line.rstrip().split('\t'))-9

                ## Get sample indexes.
                if not args.remove:
                    indexes_GT = list(range(9, n+1))
                else:
                    indexes_GT = [i for i in range(9, n+9) if line.rstrip().split()[i] not in set_remove]

                ## Print header line.
                if i_vcf == 0:
                    l = list(range(0, 9)) + indexes_GT
                    ## print header line
                    print(
                        '\t'.join([line.split()[i] for i in l]),
                        sep='\t', end='\n', file=fd_out)

                ## loop records
                for line in fd_vcf:
                    l = line.rstrip().split('\t')

                    ## Skip non-biallelic SNPs and indels
                    if args.biallelic and len(l[4].split(',')) > 1:
                         continue

#                    s = '\t'.join(l[:9])+':GL'
                    s = '\t'.join(l[:8])+'\t'+'GT:GL'
                    ## index PL or GP or ...
                    i1 = l[8].split(':').index(format1)
                    print(s, sep='\t', file=fd_out, end='\t')
                    for i_GT in indexes_GT:

                        s1 = l[i_GT].split(':')[i1]
                        if s1 == '.':
                            s2 = '.'
                        else:
#                            l_probs = [
#                                pow(10, -int(log10likelihood)/10) for log10likelihood in s_PL.split(',')]
#                            sum_prob = sum(l_probs)
#                            s_GL = ','.join(map(str, (round(prob/sum_prob, 4) for prob in l_probs)))
                            s2 = func(s1)
#                        print(
#                            ':'.join((l[i_GT], s_GL)), sep='\t', end='\t',
#                            file=fd_out)
                        print(l[0]+":"+s2, sep='\t', end='\t', file=fd_out)
                    ## Add newline character to end of line.
                    print(end='\n', file=fd_out)

    return


def PL2GL(s1):

    s2 = ','.join(
        str(-float(PL)/10) for PL in s1.split(','))

    return s2


def GP2GL(s1):

    s2 = ','.join(
        str(-abs(math.log10(float(GP)))) if GP != '0' else '-99' for GP in s1.split(','))

    return s2


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', nargs = '+')
    parser.add_argument('--out', required=True)
    parser.add_argument('--remove')
    parser.add_argument('--biallelic', action='store_true')
    parser.add_argument('--conversion', choices=('PL2GL', 'GP2GL'), default='PL2GL',)
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()
