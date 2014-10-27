#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, February-March 2013, October-November 2014

import argparse
import fileinput
import os
import gzip
import contextlib
import bz2
import itertools


def main():

    args = argparser()

    l_path_vcf_predict = sort_nicely(args.predict)

    col_predict = sampleID2col(args.sampleID_predict, l_path_vcf_predict)
    col_gold = sampleID2col(args.sampleID_gold, args.gold)

    with contextlib.ExitStack() as stack:

        fd_predict = fileinput.FileInput(
            files=l_path_vcf_predict, openhook=hook_compressed_text)

        fd_vcf_gold = hook_compressed_text(args.gold, 'r')

        if args.recal:
            fd_recal = open(args.recal)
        else:
            fd_recal = None

        if args.bed:
            fd_bed = stack.enter_context(gzip.open(args.bed, 'rt'))
        else:
            fd_bed = None

        fd_out = open(args.out, 'w')

        loop_lines(
            fd_predict, fd_vcf_gold, fd_out, fd_bed,
            col_predict, col_gold,
            fd_recal, args)
        
    return


def sampleID2col(sampleID, l_path_vcf_predict):

    with fileinput.FileInput(
            files=l_path_vcf_predict, openhook=hook_compressed_text) as fd:
        for line in fd:
            if line[:2] == '##':
                continue
            break
    col = line.rstrip().split('\t').index(sampleID)

    return col


def split_line_bed(fd_bed):

    for line in fd_bed:
        chrom, pos1, pos2 = fd_bed.readline().split('\t')
        pos1 = int(pos1)
        pos2 = int(pos2)
        yield chrom, pos1, pos2

    return


def loop_lines(
    fd_predict, fd_gold, fd_out, fd_bed,
    col1, col2, fd_recal, args):

    l_chroms = [str(i) for i in range(1,23)]+['X','Y']

    l2, chrom2, pos2 = next(split_line_vcf(fd_gold))

    if fd_bed:
        chrom_bed, pos1_bed, pos2_bed = next(split_line_bed(fd_bed))

    nFalse = 0
    nTrue = 0

    for l1, chrom1, pos1 in split_line_vcf(fd_predict):
##    while True:
##        try:
##            l1, chrom1, pos1 = next(split_line_vcf(fd_predict))
##        except StopIteration:
##            break

        ## skip SNP/INDEL
        ref1 = l1[3]
        alt1 = l1[4]
        if args.mode != 'BOTH':
            bool_INDEL = is_INDEL(ref1, alt1)
    #        if pos1 == 726481:
    #            print(ref1,alt1)
    #            print(bool_INDEL)
    #            stop
            if args.mode == 'SNP' and bool_INDEL == True:
                continue
            if args.mode == 'INDEL' and bool_INDEL == False:
                continue

        ## skip non-call and REFREF-call
        GT1 = l1[col1].split(':')[0]
        if GT1 == './.' or GT1 == '0/0':
            if args.sort == 'VQSLOD':
                ## Keep recal in sync with vcf line reading.
                next(split_line_vcf(fd_recal))
            continue

        ## bed filtering
        if fd_bed:
            while l_chroms.index(chrom1) > l_chroms.index(chrom_bed):
                chrom_bed, pos1_bed, pos2_bed = next(split_line_bed(fd_bed))
            while pos1 > pos2_bed:
                chrom_bed, pos1_bed, pos2_bed = next(split_line_bed(fd_bed))
            if pos1_bed > pos1:
                if args.sort == 'VQSLOD':
                    next(split_line_vcf(fd_recal))
                continue
            pass
                
##        VQSLOD = parse_VQSLOD(fd_recal)
        if args.sort == 'VQSLOD':
            sort = VQSLOD = parse_VQSLOD(fd_recal, chrom1, pos1)
        else:
            sort = QUAL = l1[5] 

        bool_TF = None
        bool_break = False
        bool_read = False
        while True:
            if chrom1 != chrom2:
                if l_chroms.index(chrom1) > l_chroms.index(chrom2):
                    bool_read = True
                    pass
                else:
                    bool_TF = False
                    break
                pass
            elif pos1 == pos2:
                if bool_TF != None and l2[6] != 'PASS':
                    stop1
                    pass
                GT2 = l2[col2].split(':')[0].replace('|','/')
                ref2 = l2[3]
                alt2 = l2[4]
                gt1 = GT2gt(GT1, ref1, alt1)
                gt2 = GT2gt(GT2, ref2, alt2)
                if gt1 == gt2:
                    bool_TF = True
                elif gt1 == list(reversed(gt2)):
                    bool_TF = True
                else:
                    bool_TF = False
                break
##                ## compare genotypes
####                bool_TF = compareGT(GT1,ref1,alt1,GT2,ref2,alt2,bool_TF)
##                ## do not compare genotypes/haplotypes, but only positions/sites
##                bool_TF = True
##                pass
            elif pos1 > pos2:
                bool_read = True
                pass
            ## elif pos2 > pos1:
            else:
                if bool_TF == None:
                    bool_TF = False ## FP or TN
                else:
                    stop2
                break
            if bool_read:
                try:
                    l2, chrom2, pos2 = next(split_line_vcf(fd_gold))
                except StopIteration:
                    bool_break = True
                bool_read = False
            if bool_break:
                break
            continue
        ## last variant of gold vcf
        if bool_break:
            break

        assert bool_TF != None
        assert pos1 <= pos2
        
        QUAL = l1[5]

        if len(gt1[0]) == len(ref1) and len(gt1[1]) == len(ref1):
            bool_INDEL = False
            type_variant = 'SNP'
            if not len(ref1) == 1:  # assert not MNP
                print(ref1, alt1, ref2, alt2)
                print(gt1)
                print(gt2)
                stop
        else:
            bool_INDEL = True
            type_variant = 'INDEL'

        if bool_TF == False:
            nFalse +=1
        else:
            nTrue += 1
        line = '{chrom}\t{pos}\t{bool_TF}\t'.format(
            chrom=chrom1, pos=pos1, bool_TF=bool_TF)
        line += '{}\t{}\t{}\t'.format(QUAL, sort, GT1)
        line += '{}\t{}\t'.format(ref1, alt1)
        line += '\n'
        fd_out.write(line)
        print(line)

    print('False', nFalse)
    print('True', nTrue)

##sort -k5nr,5 INDEL.phasing_annotated.vcf.gz.2548897.out \
##| awk -v nTrue=623666 -v nFalse=1955934 \
##'{if($3=="False") {FP++} else {TP++}; print FP/nFalse, TP/nTrue}'

    return


def GT2gt(GT, ref, alt):

    l = [ref]+alt.split(',')
    gt = [l[int(i)] for i in GT.split('/')]

    if ref == 'ACAAAAAAAAAC':
        print(GT, ref, alt)
        print(gt[0])
        print(gt[1])
        stop

    return gt


def parse_VQSLOD(fd_recal,chrom1,pos1):

    l_recal, chrom_recal, pos_recal = next(split_line_vcf(fd_recal))
    if chrom_recal != chrom1:
        stop
    if pos_recal != pos1:
        print(chrom_recal,pos_recal,pos1)
        stop
    for s in l_recal[7].split(';'):
        l = s.split('=')
        if l[0] == 'VQSLOD':
            VQSLOD = float(l[1])
            break

    return VQSLOD


def compareGT(GT1,ref1,alt1,GT2,ref2,alt2,bool_TF):

    d = {0:[],1:[]}
    for i, (GT, alt, ref) in enumerate(
        [(GT1,alt1,ref1),(GT2,alt2,ref2)]):
##                    print(pos1,pos2,i,GT,alt,ref)
        l = [ref]+alt.split(',')
        try:
            d[i] = [l[int(j)] for j in GT.split('/')]
        except ValueError:
            print(s)
            print(l2)
            stop5
    if d[0] == d[1] or d[0] == list(reversed(d[1])):
        if bool_TF == False:
            print('\n')
            print(pos1,pos2)
            print(GT1,ref1,alt1)
            print(GT2,ref2,alt2)
            print(l1[col1])
            print(l2[col2])
            stop4
        bool_TF = True ## TP and FN
    else:
        if set(d[0]) != set(d[1]):
            if bool_TF == True:
                print('\n')
                print(pos1,pos2)
                print(GT1,ref1,alt1)
                print(GT2,ref2,alt2)
                print(l1[col1])
                print(l2[col2])
                stop5
            bool_TF = False
        else:
            print('\n')
            print(d)
            print(pos1,ref1,alt1,GT)
            print(l2[:9])
            print(l2[col2])
            print(l1[col1])
            print(list(reversed(d[1])))
            stop1
##                break
    return bool_TF


def is_INDEL(ref,alt):
    
    ## skip deletions
    if not ref in ['A','C','G','T',]:
        if ref in iter(','.join(tup) for tup in itertools.permutations('ACGT',2)):
            stoptmp
        bool_INDEL = True
    ## diallelic SNP
    elif alt in ['A','C','G','T',]:
        bool_INDEL = False
    ## triallelic SNP
    elif alt in iter(','.join(tup) for tup in itertools.permutations('ACGT',2)):
        bool_INDEL = False
    ## multiallelic SNP
    elif alt in iter(','.join(tup) for tup in itertools.permutations('ACGT',3)):
        bool_INDEL = False
    ## skip insertions
    else:
        bool_INDEL = True

    return bool_INDEL


def split_line_bgl(fd_bgl):

    for line in fd_bgl:
        if line[:6] == 'marker':
            continue
        l = line.rstrip().split()
        chrom, pos = l[0].split(':')
        pos = int(pos)
##        ref = l[3]
##        alt = l[4]
        yield l, chrom, pos

    return


def split_line_vcf(fd_vcf):

    for line in fd_vcf:
        if line == '':
            return
        if line[0] == '#':
            continue
        l = line.rstrip().split('\t')
        chrom = l[0]
        pos = int(l[1])
##        ref = l[3]
##        alt = l[4]
        yield l, chrom, pos

    return


def hook_compressed_text(filename, mode):

    ##http://stackoverflow.com/questions/21529163/python-gzipped-fileinput-returns-binary-string-instead-of-text-string/21529243

    ext = os.path.splitext(filename)[1]
    if ext == '.gz':
        f = gzip.open(filename, mode+'t')
    elif ext == '.bz2':
        f = bz2.open(filename, mode+'t')
    else:
        f = open(filename, mode)

    return f


def alphanum_key(s):
    import re
    ## http://dave.st.germa.in/blog/2007/12/11/exception-handling-slow/
    NUM_RE = re.compile('([0-9]+)')
    return [ int(c) if c.isdigit() else c for c in NUM_RE.split(s) ]


def sort_nicely(l):
    ## http://nedbatchelder.com/blog/200712/human_sorting.html
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)
    return l


def argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument('--gold', '--truth', required=True)

    parser.add_argument('--recal', required=False)

    parser.add_argument('--predict', required=True, nargs='+')

    parser.add_argument('--sampleID_predict', required=False)

    parser.add_argument(
        '--sampleID_gold', '--sampleID_truth', required=False)

    parser.add_argument('--mode', default='BOTH', choices=['BOTH','SNP','INDEL'])

    parser.add_argument(
        '--sort', default='QUAL', choices=['QUAL', 'VQSLOD'], required=True)

    parser.add_argument('--out')

    parser.add_argument('--bed', required=False)

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()
