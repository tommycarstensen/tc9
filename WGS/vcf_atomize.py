#!/usr/bin/env python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, October-November 2014

## Use this script to normalize variants; i.e.
## 1) Split multiallelic variants into biallelic variants.
## 2) Convert MNPs to SNPs.
## 3) Trim and left align INDELs relative to reference sequence.
## 4) Remove duplicates.

import argparse
import os
import sys
import itertools
import re
import collections
import gzip
import contextlib
import operator


def main():

    args = parse_args()

    with open_file(args.vcf) as fd_vcf, \
         open_file(args.out, 'w') as fd_out, \
         open_file(args.ref) as fd_ref:

        print_header(fd_vcf, fd_out)
        print_body(fd_vcf, fd_out, fd_ref, args)

    return


def print_body(fd_vcf, fd_out, fd_ref, args):

    pattern = re.compile(r'([/|])')

    d_fai = read_fai(args.ref+'.fai')
    
    for line in fd_vcf:
        print_new_line(line, fd_out, pattern, d_fai, fd_ref, args)

    return


def print_new_line(line, fd_out, pattern, d_fai, fd_ref, args):

    REF, ALT, CHROM, POS, ID, QUAL, FILTER, INFO, FORMAT, l = split_line_vcf(
        line)

    l_ALT = ALT.split(',')

    ## 1) simple biallelic SNP
    if len(REF) == 1 and len(ALT) == 1:
        REFALT_tuples = [(POS, REF, ALT, (0,1), (0,1))]
        type_variant = 'SNP'
        pass
    ## 2) MNP or multiallelic SNP
    elif len(REF) > 1 and all([len(ALT_i) == len(REF) for ALT_i in l_ALT]):
        ## multiallelic SNP
        if len(REF) == 1:
            type_variant == 'SNP'
            print(REF, ALT)
            stop
        ## MNP
        else:
            type_variant = 'MNP'
            REFALT_tuples = []
            for POS_new, REFALT in enumerate(zip(*[REF]+l_ALT), POS):
                if not args.keep_multiallelics:
                    REFALT_tuples += [
                        (POS_new, REFALT[0], REFALT[i], (0,1), (0,i))
                        for i in range(1, len(REFALT))]
                else:
                    REFALT_tuples += [(
                        POS_new, REFALT[0], ','.join(REFALT[1:]),
                        tuple(range(len(REFALT))),
                        tuple(range(len(REFALT))))]
    ## 4) INDEL or multiallelic SNP or complex (e.g. MNP and INDEL)
    else:
        if (
            len(REF) == 1 and
            all([len(l_ALT[i]) == len(REF) for i in range(len(l_ALT))])):
            type_variant = 'SNP_multiallelic'
        elif not any([len(ALT_i) > 1 and len(ALT_i) == len(REF) for ALT_i in l_ALT]):
            type_variant = 'INDEL'
        else:
            type_variant = 'complex'
            
        if args.keep_multiallelics:
            REFALT_tuples = [(
                POS, REF, ALT,
                tuple(range(len(l_ALT)+1)),
                tuple(range(len(l_ALT)+1)))]
        else:
            REFALT_tuples = [
                (POS, REF, l_ALT[i], (0, 1), (0, i+1))
                for i in range(len(l_ALT))]
        type_variant = 'INDEL'

    for (
        i, (POS_new, REF_new, ALT_new, GT_tuple_new, GT_tuple_old)
        ) in enumerate(REFALT_tuples):
        if REF_new == ALT_new:
            continue
        POS_new, REF_new, ALT_new = trim_and_left_align(
            d_fai, fd_ref, CHROM, POS_new, REF_new, ALT_new)
        REFALT_tuples[i] = (POS_new, REF_new, ALT_new, GT_tuple_new, GT_tuple_old)

    for (
        POS_new, REF_new, ALT_new, GT_tuple_new, GT_tuple_old
        ## sort by POS_new and GT_tuple_old
        ) in sorted(REFALT_tuples, key = operator.itemgetter(0, 4)):

        if not args.keep_multiallelics:
            assert len(ALT_new.split(',')) == 1
        line_new = '{CHROM}\t{POS}\t{ID}\t{REF}\t{ALT}\t\
{QUAL}\t{FILTER}\t{INFO}\t{FORMAT}'.format(
    CHROM=CHROM, POS=POS_new, ID=ID, REF=REF_new, ALT=ALT_new,
    QUAL=QUAL, FILTER=FILTER, INFO=INFO, FORMAT='GT')
        bool_print = False
        ## loop over samples
        for GT in iter(s.split(":",1)[0] for s in l[9:]):
            if GT == './.':
                GT_new = GT
            else:
                GT1, sep, GT2 = re.split(pattern, GT)
                GT1 = int(GT1)
                GT2 = int(GT2)
                try:
                    GT1_new = GT_tuple_new[GT_tuple_old.index(GT1)]
                except ValueError:
                    GT1_new = '.'
                    sep = '/'
                try:
                    GT2_new = GT_tuple_new[GT_tuple_old.index(GT2)]
                except ValueError:
                    GT2_new = '.'
                    sep = '/'
                GT_new = '{}{}{}'.format(GT1_new, sep, GT2_new)
            if GT_new not in ('0/0','0|0','./.'):
                bool_print = True
            line_new += '\t{}'.format(GT_new)
            ## continue loop over samples
            continue
        if args.keep_monomorphics == True or bool_print:
            try:
                print(line_new, end='\n', file=fd_out, flush=True)
            except BrokenPipeError:
                sys.stderr.close()
                sys.exit()
                
    return


def print_header(fd_vcf, fd_out):

    for line in fd_vcf:
        ## print meta information lines
        if line[:2] == '##':
            print(line, end='', file=fd_out)
            continue
        ## print header line
        print(line, end='', file=fd_out)
        break

    return


def split_line_vcf(line):

    l = line.split('\t')
    REF = l[3]
    ALT = l[4]
    ## parse vcf
    CHROM = l[0]
    POS = int(l[1])
    ID = l[2]
    QUAL = l[5]
##            QUAL = '.'
    FILTER = l[6]
##            FITLER = '.'
    INFO = l[7]
##            INFO = '.'
##            FORMAT = l[8]
    FORMAT = 'GT'
##            GT = l[9].split(":",1)[0]

    return REF, ALT, CHROM, POS, ID, QUAL, FILTER, INFO, FORMAT, l


def trim_and_left_align(d_fai, fd_ref, CHROM, POS, REF, ALT):

    l_ALT = ALT.split(',')

    ## right trim and left align
    while all([REF[-1] == ALT[-1] for ALT in l_ALT]):
        ## right trim
        REF = REF[:-1]
        l_ALT = [ALT[:-1] for ALT in l_ALT]
        ## left align
        if len(REF) == 0 or any([len(ALT) == 0 for ALT in l_ALT]):
            POS -= 1
            nt = parse_ref(d_fai, fd_ref, CHROM, POS)
            REF = nt+REF
            l_ALT = [nt+ALT for ALT in l_ALT]

    ## left trim
    while all([REF[0] == ALT[0] for ALT in l_ALT]):
        if len(REF) == 1 or any([len(ALT) == 1 for ALT in l_ALT]):
            break
        POS += 1
        REF = REF[1:]
        l_ALT = [ALT[1:] for ALT in l_ALT]

    ALT = ','.join(l_ALT)

    return POS, REF, ALT


def parse_ref(d_fai, fd_ref, CHROM, POS, size=1):

    ## parse ref
    cnt = cnt_chars_per_line_excluding_newline = 60
    row1 = (POS-1)//cnt
    row2 = (POS-1+size)//cnt
    size += row2-row1
    col = (POS-1)%cnt
    byte_init = d_fai[CHROM]['start']
    offset = byte_init+(cnt+1)*row1+col
    fd_ref.seek(offset)
    read = fd_ref.read(size).replace('\n','')

    return read


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
    parser.add_argument(
        '--keep_mnps', help='Convert MNPs to SNPs by default',
        default=False, action='store_true')
    parser.add_argument(
        '--keep_multiallelics',
        help='Split multiallelics to biallelics by default',
        default=False, action='store_true')

    parser.add_argument(
        '--keep_monomorphics', help='Do not print monomorphic sites by default',
        default=False, action='store_true')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()
