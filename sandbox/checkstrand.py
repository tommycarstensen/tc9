#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute
## December 2013, January 2015

import argparse
import os
import datetime
import sys


def main():

    args = argparser()

    l_fam_sampleIDs = parse_fam(args.bfile + '.fam')

    ## Get count of samples.
    n_samples = len(l_fam_sampleIDs)
    ## Get count of SNPs. Only needed to keep track of completion percent.
    with open(args.bfile + '.bim', 'r') as bim:
        n_SNPs = len(bim.readlines())

    d_fai = read_fai(args.ref + '.fai')

    with open(args.bfile + '.bim', 'r') as bim, \
         open(args.ref, 'r') as ref:
        convert(
            args, bim, ref, d_fai,
            n_samples, n_SNPs)

    return

def convert(
    args, bim, ref, d_fai,
    n_samples, n_SNPs):

    ## data lines
    for i_SNP, line_bim in enumerate(bim):
        ## By default, the minor allele is coded A1
        ## and the major allele is coded A2
        CHROM, ID, _, POS, A1, A2 = line_bim.rstrip().split()
        if args.chrom and CHROM != args.chrom:
            continue
        if i_SNP % 100000 == 0:
            print('SNP {} of {} SNPs. CHROM={} POS={} time={}'.format(
                i_SNP, n_SNPs, CHROM, POS,
                datetime.datetime.now().strftime("%H:%M:%S")))

        ## Parse reference allele from reference sequence.
        REF = parse_ref(d_fai, ref, CHROM, int(POS))

        ## ALT = major
        ## Convert MAF to allele frequency for the ALT allele.
        if REF == A1:
            ALT = A2
        ## ALT = minor
        elif REF == A2:
            ALT = A1
        elif A1 == '0':
            ## 
            if not AF_A1 == 0:
                print('WARNING: REF={} CHROM={} POS={} ID={} A1={} A2={} AF={:.4f}'.format(
                    REF, CHROM, POS, ID, A1, A2, AF_A1), file=sys.stderr)
                ALT = A2
                continue
            ## monomorphic
            else:
                ALT = '.'
        ## Neither A1 nor A2 is identical to the reference allele.
        ## Should not happen. Print an error and continue.
        else:
            print(REF, A1, A2, CHROM, POS)
##            stop

    return


def read_fai(path_fai):

    assert os.path.isfile(path_fai)

    d = {}
    with open(path_fai) as file_fai:
        for line_fai in file_fai:
            l = line_fai.rstrip().split()
            chrom = l[0]
            byte_length = int(l[1])
            byte_start = int(l[2])
            bytes_per_line_excl_line_break = int(l[3])
            bytes_per_line_incl_line_break = int(l[4])
            d[chrom] = {'length': byte_length, 'start': byte_start}

    return d


def parse_ref(d_fai, fd_ref, CHROM, POS, size=1):

    ## parse ref
    cnt = cnt_chars_per_line_excluding_newline = 60
    row1 = (POS - 1) // cnt
    row2 = (POS - 1 + size) // cnt
    size += row2 - row1
    col = (POS - 1) % cnt
    byte_init = d_fai[CHROM]['start']
    offset = byte_init + (cnt + 1) * row1 + col
    fd_ref.seek(offset)
    read = fd_ref.read(size).replace('\n', '')

    return read


def parse_fam(fp_fam,):
    
    with open(fp_fam, 'r') as fam:
        l_sampleIDs = [line.rstrip().split()[0] for line in fam]
    with open(fp_fam, 'r') as fam:
        line = fam.readline()
        if line.split()[0] != line.rstrip().split()[1]:
            print('havent written code for when FID and IID are different')
            sys.exit()

    if len(l_sampleIDs) != len(set(l_sampleIDs)):
        print('duplicate sample IDs')
        print(l_sampleIDs)
        stop

    return l_sampleIDs


def argparser():

    parser = argparse.ArgumentParser()

    ## Required.
    parser.add_argument('--bfile', '--in', required = True)
    s_help = 'The A2 allele is saved as the reference by PLINK2, '
    s_help = 'but a reference sequence is needed to determine REF and ALT.'
    parser.add_argument(
        '--ref', required = True, default=None,
        help=s_help,
        )

    parser.add_argument('--chrom', required = False)

    args = namespace_args = parser.parse_args()

    ## Assert that input exists.
    assert all([
        os.path.isfile('{}.{}'.format(args.bfile, suffix))
        for suffix in ('bed','bim','fam')])
    assert os.path.isfile(args.ref)

    ## Do not allow compressed reference sequence for now.
    assert os.path.splitext(args.ref)[1] == '.fa'

    return args


if __name__ == '__main__':
    main()
