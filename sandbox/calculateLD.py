#!python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, February 2014

## http://en.wikipedia.org/wiki/Linkage_disequilibrium

## the script currently does not allow multiple chromosomes to be present in a single haps file
## currently just print above threshold to file - identical to LDselect.py

import argparse
import collections
import os
import sys
import fileinput
import contextlib
import gzip
from Bio import bgzf
import itertools


def main():

    d_args = argparser()

##    l_haps = d_args['haps']
##    threshold_min_LD = d_args['threshold_min_LD']
##    threshold_min_MAF = d_args['threshold_min_MAF']
##    path_out = d_args['out']
##    bool_verbose = d_args['verbose']

    with open_haps(d_args['haps']) as f:
        n = len(f.readline().rstrip().split())-5

    if d_args['sample_columns']:
        str_ = d_args['sample_columns']
        cols = list(itertools.chain.from_iterable(
            list(range(
                int(str_.split(',')[i].split('-')[0])-1,
                int(str_.split(',')[i].split('-')[-1])))
            for i in range(str_.count(',')+1)))
    else:
        cols = list(range(n))

    print('indexing', d_args['haps'])
    d_tell = {}
    pos_min = float('inf')
    pos_max = -float('inf')
##    d_args['pos1'] = int(d_args['pos1'])
##    d_args['pos2'] = int(d_args['pos2'])
    with open_haps(d_args['haps']) as f:
        for line in f:
            tell = f.tell()
##            key = tuple(f.readline().split(None,5)[:5])
##            pos = int(key[2])
            ID, pos, alleleA, alleleB, l_haps = parse_line(line, cols)
            ## skip INDELs (and non-biallelic SNPs)
            if alleleA in ['0','1'] and alleleB in ['0','1']:
                continue
            key = (pos,alleleA,alleleB)
            if pos < d_args['pos1']:
                continue
            if pos > d_args['pos2']+d_args['window']:
                break
            d_tell[key] = tell
            if pos < pos_min:
                pos_min = pos
            elif pos > pos_max:
                pos_max = pos
    print('indexed', d_args['haps'])

    with contextlib.ExitStack() as stack:

##        fd_out = stack.enter_context(
##            open(d_args['out'],'w'))

        fd_haps1 = stack.enter_context(open_haps(d_args['haps']))
        fd_haps2 = stack.enter_context(open_haps(d_args['haps']))
        if d_args['out'][-3:] == '.gz':
            fd_out = stack.enter_context(gzip.open(d_args['out'],'w'))
        else:
            fd_out = stack.enter_context(open(d_args['out'],'w'))

        pos_min = d_args['pos1']
        pos_max = d_args['pos2']
        window = d_args['window']

        print('looping', pos_min, '-', pos_max, d_args['haps'])
        for line1 in fd_haps1:
            ID1,pos1,A1,B1,l_haps1 = parse_line(line1, cols)
            ## skip INDELs (and non-biallelic SNPs)
            if A1 in ['0','1'] and B1 in ['0','1']:
                continue
            ## skip outside of range
            if pos1 < pos_min:
                continue
            print(pos1,A1,B1)
            if pos1 > pos_max:
                break
            MAF1 = (2*min(l_haps1.count('00'),l_haps1.count('11'))+l_haps1.count('01')+l_haps1.count('10'))/n
            if MAF1 > 0.5:
                print('MAF1',MAF1)
                stop8
            if MAF1 < d_args['min_MAF']:
                continue
            fd_haps2.seek(d_tell[pos1,A1,B1])
            for line2 in fd_haps2:
                ID2,pos2,A2,B2,l_haps2 = parse_line(line2, cols)
                ## skip INDELs (and non-biallelic SNPs)
                if A2 in ['0','1'] and B2 in ['0','1']:
                    continue
                ## skip outside of range
                if pos2 > pos1+window:
                    break
                if pos2 > pos_max+window:
                    print(pos1,pos2,pos_max+window,pos_max,window)
                    stoptmp9
                    break
                MAF2 = (2*min(l_haps2.count('00'),l_haps2.count('11'))+l_haps2.count('01')+l_haps2.count('10'))/n
                if MAF2 < d_args['min_MAF']:
                    continue
                if MAF2 > 0.5:
                    print('MAF2',MAF2)
                    stop8
                if pos1 == pos2:
                    if ID1 == ID2:
                        continue
                    else:
                        print(ID1,ID2,pos1,pos2)
                        stoptmp
                r2 = calc_LD(l_haps1,l_haps2)
                if r2 < d_args['min_LD']:
                    continue
                fd_out.write(
                    '{ID1} {ID2} {pos1} {pos2} {r2:.4f} {MAF1:.3f} {MAF2:.3f}\n'.format(
                        ID1=ID1,ID2=ID2,pos1=pos1,pos2=pos2,r2=r2,MAF1=MAF1, MAF2=MAF2))

##                ## monomorphic (or rare error l_haps1=['10', '01', '01'], l_haps2=['00', '10', '01'])
##                if not r2:
##                    continue

##                ## haplotypes in LD            
##                if(r2>=d_args['min_LD']):
##    ##                print('aaa', pos1, pos2, r2, p1, p2, q1, q2)
##    ##                print('aaa', pos1, pos2)
##                    continue

##                if q1 < d_args['min_MAF'] or 1-q1 < d_args['min_MAF']:
##                    stop3
##                    continue

    return


def main_old():

    (
        l_haps, path_out, bool_verbose, threshold_min_LD, threshold_min_MAF
        ) = argparser()

    with contextlib.ExitStack() as stack:

        fd_haps = stack.enter_context(
            fileinput.FileInput(files=l_haps,openhook=hook_compressed_text))
        fd_out = stack.enter_context(open(path_out,'w'))

        ## while monomorphic
        while True:
            line = fd_haps.readline()
            ID1,pos1,l_haps1 = parse_line(line)
            ## sample count
            n = len(l_haps1)
            ## not monomorphic
            if l_haps1 != n*['00'] and l_haps1 != n*['11']:
                ## not below MAF threshold
                MAF = (min(l_haps1.count('00'),l_haps1.count('11'))+l_haps1.count('01')+l_haps1.count('10'))/n
                print('MAF',MAF)
                if MAF > threshold_min_MAF:
                    print('MAF',MAF)
                    break

        ## loop after monomorphics
        for line in fd_haps:
            ID2,pos2,l_haps2 = parse_line(line)
            if pos2 and pos2%10000 == 0 and bool_verbose == True:
                print('position',pos2)

            r2 = calc_LD(l_haps1,l_haps2)

            ## monomorphic (or rare error l_haps1=['10', '01', '01'], l_haps2=['00', '10', '01'])
            if not r2:
                continue

            ## haplotypes in LD            
            if(r2>=threshold_min_LD):
##                print('aaa', pos1, pos2, r2, p1, p2, q1, q2)
##                print('aaa', pos1, pos2)
                continue

            if q1 < threshold_min_MAF or q2 < threshold_min_MAF:
                continue

##            if q1 > 0.49 and q2 < 0.01:
##                print(r2,q1,q2)
##                print(l_haps1)
##                print(l_haps2)
##                stop
            ## haplotypes not in LD
##            print('bbb', pos1, pos2, r2, p1, p2, q1, q2)
            fd_out.write('%s %s %s %s %s\n' %(ID2,pos2,ID1,pos1,r2))
            ID1 = ID2
            pos1 = pos2
            l_haps1 = l_haps2
                
    return


def open_haps(pathname):

    ##http://stackoverflow.com/questions/21529163/python-gzipped-fileinput-returns-binary-string-instead-of-text-string/21529243

    ext = os.path.splitext(pathname)[1]
    if ext == '.gz':
        stop_convert_to_bgz
        f = gzip.open(pathname, 'rt')
    elif ext == '.bgz':
        f = bgzf.BgzfReader(pathname)
    else:
        f = open(filename)

    return f


def calc_LD(l_haps1,l_haps2):

    if l_haps1 == l_haps2:
        r2 = 1
    else:
        d_cnt_haplotype = count_haplotypes(l_haps1,l_haps2,)
        if d_cnt_haplotype['11'] == 0 and d_cnt_haplotype['01'] == 0:
            r2 = p1 = p2 = q1 = q2 = None
            stoptmp7
        else:
            r2 = calc_r2(d_cnt_haplotype)

    return r2


def calc_r2(d_cnt_haplotype):

    ## nomenclature: http://en.wikipedia.org/wiki/Linkage_disequilibrium

##    ## unit tests
##    d_cnt_haplotype = {
##        '00':474,
##        '01':611,
##        '10':142,
##        '11':773,
##        }
####    r2 = 0.09239127328988425
####    d_cnt_hap = {'01': 41, '00': 205, '10': 440, '11': 88}
######    r2 = 0.0

    ## haplotype frequencies
    x11 = d_cnt_haplotype['00'] ## A1B1
    x12 = d_cnt_haplotype['01'] ## A1B2
    x21 = d_cnt_haplotype['10'] ## A2B1
    x22 = d_cnt_haplotype['11'] ## A2B2

##    ## allele frequencies
####    p1 = d_cnt_allele['p']['0']/(4*n)
####    p2 = d_cnt_allele['p']['1']/(4*n)
####    q1 = d_cnt_allele['q']['0']/(4*n)
####    q2 = d_cnt_allele['q']['1']/(4*n)
##    p1 = x11+x12
##    p2 = x21+x22
##    q1 = x11+x21
##    q2 = x12+x22

##    ## coefficient of linkage disequilibrium (Lewontin and Kojima 1960)
##    ## "coupling gametes" - "repulsion gametes"
##    ## actual/observed haplotype frequency - expected haplotype frequency if independence
##    ## D = x11-p1*q1 and hence x11=p1*q1+D and hence
##    D = x11*x22-x12*x21

##    if not D: ## e.g. l_haps1 ['10', '01', '01'] and l_haps2 ['00', '10', '01']
    ## SNP1: 00, SNP2: 11 not present
    if x11*x22 == x12*x21:
        r2 = 0
##    if not D: ## e.g. l_haps1 ['10', '01', '01'] and l_haps2 ['00', '10', '01']
####        print('x11,x12,x21,x22',x11,x12,x21,x22)
####        print(
####            d_cnt_haplotype['00'],
####            d_cnt_haplotype['01'],
####            d_cnt_haplotype['10'],
####            d_cnt_haplotype['11'],
####            )
####        print('D', D)
####        stop
##        r2 = None
    else:
        ## r2 measure of LD (Hill and Robertson 1968)
##        r2 = D**2/(p1*p2*q1*q2)
        r2 = (x11*x22-x12*x21)**2/((x11+x12)*(x21+x22)*(x11+x21)*(x12+x22))
    
    return r2


def count_haplotypes(l_haps1,l_haps2):

    ## haplotype count
    d_cnt_haplotype = {
        '00':0,
        '01':0,
        '10':0,
        '11':0,
        }
##    ## allele count
##    d_cnt_allele = {
##        'p':{'0':0,'1':0,},
##        'q':{'0':0,'1':0,},
##        }

    c = collections.Counter(zip(l_haps1,l_haps2,))
    for t in c.items():
        ## nomenclature: position/SNP 1 and 2, alleles/strands A and B
        allele1A = t[0][0][0]
        allele1B = t[0][0][-1]
        allele2A = t[0][1][0]
        allele2B = t[0][1][-1]
        haplotypeA = allele1A+allele2A
        haplotypeB = allele1B+allele2B
        d_cnt_haplotype[haplotypeA] += t[1]
        d_cnt_haplotype[haplotypeB] += t[1]
##        d_cnt_allele['p'][allele1A] += t[1]
##        d_cnt_allele['p'][allele1B] += t[1]
##        d_cnt_allele['q'][allele2A] += t[1]
##        d_cnt_allele['q'][allele2B] += t[1]

    return d_cnt_haplotype


def parse_line(line, cols):

    ## http://shapeit.fr/pages/m02_formats/hapssample.html

    l = line.split()
##    if l[0] in ['0','1']:
##        ID = None
##        pos = None
##        l_EUR = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,357,395,396,397,398,399,400,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,994,995,996,997,998,999,1000,1001,1002,1003,1004,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1034,1035,1036,1037,1038,1039,1040,1041,1042,1043,1044,1045,1046,1047,1048,1049,1050,1051,1052,1053,1054,1055,1056,1057,1058,1059,1060,1061,1062,1063,1064,1065,1066,1067,1068,1069,1070,1071,1072,1073,1074,1075,1076,1077,1078,1079,1080,1081,1082,1083,1084,1085,1086,1087,1088,1089,1090,1091]
##        l_haps = [''.join(l[2*i:2*i+2]) for i in l_EUR]
##    else:
##        ID = l[1]
##        pos = int(l[2])
##        l_haps = [''.join(l[i:i+2]) for i in range(5,len(l)-1,2)]
    ID = l[1]
    pos = int(l[2])
    alleleA = l[3]
    alleleB = l[4]
##    l_haps = [''.join(l[i:i+2]) for i in range(5,len(l)-1,2)]
    l_haps = [''.join(l[5+2*i:5+2*i+2]) for i in cols]

    return ID, pos, alleleA, alleleB, l_haps


def argparser():

    parser = argparse.ArgumentParser()

    s='SHAPEIT2 haps file; i.e. chrom rsID pos alleleA alleleB'
    s+=' (http://shapeit.fr/pages/m02_formats/hapssample.html);'
    s+= ' e.g. /lustre/scratch113/projects/agv/users/tc9/HiSeq/imputation_afrrefpan/'
    s+= 'out_SHAPEIT/Sotho_subsamp.postQC.autosomes.chrom22.haps'
    parser.add_argument('--haps', '--in', required=True, help=s)
    parser.add_argument('--out', required=True)
    parser.add_argument('--pos1', default=0, type=int)
    parser.add_argument('--pos2', default=float('inf'), type=int)
    parser.add_argument('--min_LD', type=float, default=.8)
    parser.add_argument('--min_MAF', type=float, default=.01)
    parser.add_argument('--verbose', action='store_true',)
    parser.add_argument('--sample_columns', help='e.g. 101-148,249-320')
    parser.add_argument('--window', type=int, default=250000)

    d_args = vars(parser.parse_args())

    return d_args


if __name__ == '__main__':
    main()
