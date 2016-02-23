#!/bin/env python3

import argparse
from Bio.bgzf import BgzfWriter
import gzip
import sys

def main():

    args = parse_args()

#    with BgzfWriter(args.out, 'wb') as fd_out:
    with gzip.open(args.out, 'wt') as fd_out:
        for i_vcf, vcf in enumerate(args.vcf):
            with gzip.open(vcf, 'rt') as fd_vcf:
                ## loop metadata lines and header line
                for line in fd_vcf:
                    ## end of metadata lines
                    if line[0:2] != '##':
                        ## print header line
                        print(line, end='', file=fd_out)
                        ## break loop at the header line
                        break
                    if i_vcf == 0:
                        ## print metadata line
                        print(line, end='', file=fd_out)
                ## loop records
                for line in fd_vcf:
                    l = line.rstrip().split('\t')
#                    s = '\t'.join(l[:9])+':GL'
                    s = '\t'.join(l[:8])+'\t'+'GL'
                    ## index PL
                    i_PL = l[8].split(':').index('PL')
                    print(s, sep='\t', file=fd_out, end='\t')
                    for i_GT in range(9, len(l)):
                        s_PL = l[i_GT].split(':')[i_PL]
                        if s_PL == '.':
                            s_GL = '.'
                        else:
#                            l_probs = [
#                                pow(10, -int(log10likelihood)/10) for log10likelihood in s_PL.split(',')]
#                            sum_prob = sum(l_probs)
#                            s_GL = ','.join(map(str, (round(prob/sum_prob, 4) for prob in l_probs)))
                            s_GL = ','.join(
                                str(-float(PL)/10) for PL in s_PL.split(','))
#                        print(
#                            ':'.join((l[i_GT], s_GL)), sep='\t', end='\t',
#                            file=fd_out)
                        print(s_GL, sep='\t', end='\t', file=fd_out)
                    print(end='\n', file=fd_out)

    return


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', nargs = '+')
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()
