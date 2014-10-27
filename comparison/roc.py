#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, February-March 2013

import argparse
import fileinput
import os
import gzip
import contextlib
import bz2
import itertools


def main():

    d_args = argparser()

    l_path_vcf_predict = sort_nicely(d_args['predict'])

    if d_args['tranches']:
        d_tranches = parse_tranches(d_args['tranches'])
    else:
        d_tranches = None

    col_predict = sampleID2col(d_args['sampleID_predict'], l_path_vcf_predict)
    col_gold = sampleID2col(d_args['sampleID_gold'], d_args['gold'])

    with contextlib.ExitStack() as stack:

        fd_predict = fileinput.FileInput(
            files=l_path_vcf_predict, openhook=hook_compressed_text)

        fd_vcf_gold = hook_compressed_text(d_args['gold'],'r')

        if d_args['recal']:
            fd_recal = open(d_args['recal'])
        else:
            fd_recal = None

        fd_out = open(d_args['out'],'w')

        loop_lines(
            fd_predict, fd_vcf_gold, fd_out,
            col_predict, col_gold, d_args['mode'],
            fd_recal, d_tranches, d_args)
        
    return


def parse_tranches(path_tranches):

    d_tranches = {}
    with open(path_tranches) as fd_tranches:
        for line in fd_tranches:
            if line[0] == '#':
                continue
            l = line.split(',')
            if l[0] == 'targetTruthSensitivity':
                index_targetTruthSensitivity = l.index('targetTruthSensitivity')
                index_minVQSLod = l.index('minVQSLod')
                continue
            targetTruthSensitivity = l[index_targetTruthSensitivity]
            minVQSLod = float(l[index_minVQSLod])
            d_tranches[targetTruthSensitivity] = minVQSLod

    return d_tranches


def sampleID2col(sampleID, l_path_vcf_predict):

    with fileinput.FileInput(
            files=l_path_vcf_predict, openhook=hook_compressed_text) as fd:
        line = fd.readline()
        if line[0] != '#':
            pass
        else:
            while line.split('\t')[0] != '#CHROM':
                line = fd.readline()
    try:
        if line[0] != '#':
            col = line.rstrip().split(' ').index(sampleID)
        else:
            col = line.rstrip().split('\t').index(sampleID)
    except:
        print(line)
        print(sampleID)
        stop

    return col


def loop_lines(
    fd_predict, fd_gold, fd_out, col1, col2, mode, fd_recal, d_tranches, d_args):

    predict = ['bgl','vcf'][1]  # never got it to work with bgl because no QUAL...

    l_chroms = [str(i) for i in range(1,23)]+['X','Y']

    l_FILTER_pass = ['PASS','PHI','PHI;PHR']
    l_FILTER_pass = ['.','PASS']

    ##FILTER=<ID=OC,Description="Coverage threshold exceeded">
    ##FILTER=<ID=RC,Description="RTG variant is a complex region">
    ##FILTER=<ID=RX,Description="RTG variant contained within hypercomplex region">
    ##FILTER=<ID=RCEQUIV,Description="RTG variant is equivalent to the previous variant">
    ##FILTER=<ID=OTHER,Description="Variant is invalid for unknown reasons">
    ##FILTER=<ID=PHI,Description="This variant has a phasing incompatibility">
    ##FILTER=<ID=PHR,Description="This variant should be ignored as it has been replaced by a variant with a repair for phasing incompatibility">
    ##FILTER=<ID=PME,Description="This variant should be ignored as the genotype ploidy of some of the samples did not match the expected ploidy">
    FILTER = None
##    while FILTER not in :
    while FILTER not in l_FILTER_pass:
        l2, chrom2, pos2 = next(split_line_vcf(fd_gold))
        FILTER = l2[6]

    nFalse = 0
    nTrue = 0

##    for l1, chrom1, pos1 in split_line_vcf(fd_predict):
    while True:
        if predict == 'vcf':
            try:
                l1, chrom1, pos1 = next(split_line_vcf(fd_predict))
            except StopIteration:
                break
        else:
            l1, chrom1, pos1 = next(split_line_bgl(fd_predict))

        ## skip SNP/INDEL
        if predict == 'vcf':
            ref1 = l1[3]
            alt1 = l1[4]
        else:
            ref1 = l1[1]
            alt1 = l1[2]
        bool_INDEL = is_INDEL(ref1, alt1)
#        if pos1 == 726481:
#            print(ref1,alt1)
#            print(bool_INDEL)
#            stop
        if mode == 'SNP' and bool_INDEL == True:
            continue
        if mode == 'INDEL' and bool_INDEL == False:
            continue

        ## skip non-call and reference-call
        ## http://www.1000genomes.org/node/101
        ## "The first sub-field must always be the genotype (GT)."
        GT1 = l1[col1].split(':')[0]
        if GT1 == './.' or GT1 == '0/0':
            if d_args['sort'] == 'VQSLOD':
                ## Keep recal in sync with vcf line reading.
                next(split_line_vcf(fd_recal))
            continue

##        VQSLOD = parse_VQSLOD(fd_recal)
        if d_args['sort'] == 'VQSLOD':
            sort = VQSLOD = parse_VQSLOD(fd_recal,chrom1,pos1)
        else:
            sort = QUAL = l1[5] 

        bool_TF = None
        bool_break = False
        while True:
            if chrom1 != chrom2:
                if l_chroms.index(chrom1) > l_chroms.index(chrom2):
                    pass
                else:
                    bool_TF = False
                    break
                pass
            elif pos1 == pos2:
                if bool_TF != None and l2[6] != 'PASS':
                    pass
                else:
                    GT2 = l2[col2].split(':')[0].replace('|','/')
                    ref2 = l2[3]
                    alt2 = l2[4]
                    ## compare genotypes
##                    bool_TF = compareGT(GT1,ref1,alt1,GT2,ref2,alt2,bool_TF)
                    ## do not compare genotypes/haplotypes, but only positions/sites
                    bool_TF = True
                    pass
                pass
            elif pos1 > pos2:
                pass
            ## elif pos2 > pos1:
            else:
                if bool_TF == None:
                    bool_TF = False ## FP or TN
                break
            FILTER = None
            while FILTER not in l_FILTER_pass:
                try:
                    l2, chrom2, pos2 = next(split_line_vcf(fd_gold))
                except StopIteration:
                    bool_break = True
                    break
                FILTER = l2[6]
            if bool_break:
                break
            continue
        if bool_break:
            break


        if bool_TF == None: stoptmpunexpected

        QUAL = l1[5]

        ## this is wrong
        if mode == 'INDEL':
            l_refalt = [ref1]+alt1.split(',')
            l_GT = GT1.split('/')
##            l_lengths = [len(l_refalt[int(i)])-len(ref1) for i in l_GT]
            if l_GT[0] == l_GT[1]:
                bool_het = False
            else:
                bool_het = True
        ## this is wrong
        else:
##            l_lengths = [0]
            if GT1 == '0/0' or GT1 == '1/1' or GT1 == '0|0' or GT1 == '1|1':
                bool_het = False
            else:
                bool_het = True
        l_lengths = [None]
        for lengthINDEL in l_lengths:
            if bool_TF == False:
                nFalse +=1
            else:
                nTrue += 1
            line = '{chrom}\t{pos}\t{bool_TF}\t{QUAL}'.format(
                chrom=chrom1, pos=pos1, bool_TF=bool_TF, QUAL=QUAL,)
##            line += '\t{sort}\t{length}\t{bool_het}\n'.format(
##                sort=sort, length=lengthINDEL, bool_het=bool_het)
            line += '\t{sort}\n'.format(
                sort=sort)
            fd_out.write(line)

    print('False', nFalse)
    print('True', nTrue)

##sort -k5nr,5 INDEL.phasing_annotated.vcf.gz.2548897.out \
##| awk -v nTrue=623666 -v nFalse=1955934 \
##'{if($3=="False") {FP++} else {TP++}; print FP/nFalse, TP/nTrue}'

    return


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

    parser.add_argument('--recal')

    parser.add_argument('--tranches')

    parser.add_argument('--predict', required=True, nargs='+')

    parser.add_argument('--sampleID_predict')

    parser.add_argument('--sampleID_gold', '--sampleID_truth')

    parser.add_argument('--mode', choices=['BOTH','SNP','INDEL'])

    parser.add_argument('--sort', default='QUAL', choices=['QUAL', 'VQSLOD'])

    parser.add_argument('--out')

##    print(vars(parser.parse_args()))
    d_args = {}
    for name,value in vars(parser.parse_args()).items():
        d_args[name] = value

    return d_args


if __name__ == '__main__':
    main()
