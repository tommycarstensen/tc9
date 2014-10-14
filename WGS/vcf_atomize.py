#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, October 2014

## Use this script to normalize variants; i.e.
## 1) Split multiallelic variants into biallelics.
## 2) Convert MNPs to SNPs.
## 3) Trim and left align INDELs relative to reference sequence.
## 4) Remove duplicates.

## currently lots of limitations
## currenty doesn't trim and left align
## currently only works for unphased vcf files (just use re.split for haps)
## currently only works for a single sample

import argparse
import os
import sys
import itertools
import re
import collections
import gzip


def main():

    d_args = parse_args()

    loop_vcf(d_args)

    return


def loop_vcf(d_args):

    pattern = re.compile(r'([/|])')

    d_fai = read_fai(d_args['ref']+'.fai')

    with open_file(d_args['vcf']) as fd_vcf, \
         open_file(d_args['out'], 'w') as fd_out, \
         open_file(d_args['ref']) as fd_ref:
        for line in fd_vcf:
            ## skip meta information lines
            if line[:2] == '##':
                print(line, end='', file=fd_out)
                continue
            ## skip header line
            print(line, end='', file=fd_out)
            break
        for line in fd_vcf:
            l = line.split('\t')
            REF = l[3]
            ALT = l[4]
            ## 1) simple biallelic SNP
            if len(REF) == 1 and len(ALT) == 1:
                print(line, end='', file=fd_out, flush=True)
                continue
            l_ALT = ALT.split(',')
            ## parse vcf
            CHROM = l[0]
            POS = int(l[1])
            ID = l[2]
            QUAL = l[5]
            FILTER = l[6]
            INFO = l[7]
            FORMAT = l[8]
            GT = l[9].split(":",1)[0]

##            ## parse ref
##            byte_init = d_fai[CHROM]['start']
##            quotient = (POS-1)//60
##            fd_ref.seek(-1+byte_init+POS+quotient)
##            read = fd_ref.read(len(REF)+(POS-1+len(REF))//60-quotient).replace('\n','')
##            assert REF == read

            ## 2) MNP or multiallelic SNP
            if all([len(ALT_i) == len(REF) for ALT_i in l_ALT]):
                REFALT_tuples = list(zip(*[REF]+l_ALT))
                type_variant = 'MNP'
            ## 3) INDEL
            else:
                REFALT_tuples = [tuple([REF]+l_ALT,)]

            for POS_new, REFALT_tuple in enumerate(REFALT_tuples, POS):
                d_GT_old2new = collections.OrderedDict(
                    (i, REFALT_tuple.index(REFALT))
                    for i, REFALT in enumerate(REFALT_tuple))
                REF_new = REFALT_tuple[0]

                for i, i_new in enumerate(
                    sorted(set(d_GT_old2new.values()))):
                    ALT_new = REFALT_tuple[i_new]
                    if REF_new == ALT_new:
                        continue
                    line_new = '{CHROM}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\t{FORMAT}'.format(
                        CHROM=CHROM, POS=POS_new, ID=ID,
                        REF=REF_new, ALT=ALT_new, QUAL='.', FILTER=FILTER,
                        INFO='.', FORMAT='GT')
                    bool_print = False
                    ## loop over samples
                    for GT in iter(s.split(":",1)[0] for s in l[9:]):
                        if GT == './.':
                            GT_new = GT
                        else:
                            GT1, sep, GT2 = re.split(pattern, GT)
                            GT1 = int(GT1)
                            GT2 = int(GT2)
                            if GT1 == 0 and GT2 == 0:
                                GT_new = '0{}0'.format(sep)
                            elif d_GT_old2new[GT1] == i_new and d_GT_old2new[GT2] == i_new:
                                GT_new = '1{}1'.format(sep)
                            elif d_GT_old2new[GT1] == i_new:
                                GT_new = '1{}0'.format(sep)
                            elif d_GT_old2new[GT2] == i_new:
                                GT_new = '0{}1'.format(sep)
                            else:
                                GT_new = '0/0'
                        if GT_new != './.':
                            bool_print = True
                        line_new += '\t{}'.format(GT_new)
                    if bool_print:
                        print(line_new, end='\n', file=fd_out, flush=True)

    return


def open_file(file, mode='r'):

    if file == '-':
        if mode == 'w':
            fd = sys.stdout
        else:
            fd = sys.stdin
    elif os.path.splitext(file)[-1] == '.gz':
        fd = gzip.open(file, mode+'t')
    else:
        fd = open(file, mode)

    return fd


def read_fai(path_fai):

    d = {}
    with open(path_fai) as file_fai:
        for line_fai in file_fai:
            l = line_fai.rstrip().split()
            chrom = l[0]
            byte_length = int(l[1])
            byte_start = int(l[2])
            bytes_per_line_excl_line_break = int(l[3])
            bytes_per_line_incl_line_break = int(l[4])
            d[chrom] = {'length':byte_length,'start':byte_start}

    return d


def isfile(arg):

##    if not os.path.isfile(arg) and not os.path.islink(arg):
    if arg == '-':
        pass
    elif not os.path.exists(arg):
        raise argparse.ArgumentTypeError('{} is neither a file nor a symbolic link.'.format(arg))
    else:
        pass

    return arg


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', type=isfile, default='-')
    parser.add_argument('--ref', required=True, type=isfile)
    parser.add_argument('--out', type=isfile, default='-')
    parser.add_argument('--remove_monomorphic', action='store_true')

    return dict(vars(parser.parse_args()))


if __name__ == '__main__':
    main()
