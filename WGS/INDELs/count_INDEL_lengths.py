#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, 2013

import glob
import fileinput
import re
##import rpy
##from rpy import r
import sys
sys.path.append('/nfs/users/nfs_t/tc9/github/sandbox')
import gnuplot

def main():

    pop = 'union'
    coverage = '4x'
    ts_filter_level = 99.
    print(pop,coverage,ts_filter_level)

    ## parse the input tranches file describing where to cut the data
    fp_tranches = '../pipeline/%s%s/out_VariantRecalibrator/VariantRecalibrator.INDEL.tranches' %(pop,coverage)
    minVQSLOD = parse_minVQSLOD(ts_filter_level,fp_tranches)

    if pop == 'union':
        l_vcfs = glob.glob('../pipeline/%s%s/out_CombineVariants/*.vcf' %(pop,coverage,))
    else:
        l_vcfs = glob.glob('../pipeline/%s%s/out_UnifiedGenotyper/*.vcf' %(pop,coverage,))
    l_vcfs_sorted = sort_nicely(l_vcfs)
    l_vcfs_sorted = l_vcfs_sorted[:10]
    l_vcfs_sorted = ['/nfs/t149_1kg/phase1_v3/ALL.chr1.exome.1000GApr12.vcf']

    fp_recal = '../pipeline/%s%s/out_VariantRecalibrator/VariantRecalibrator.INDEL.recal' %(pop,coverage)
    print('fp_recal',fp_recal)
    fd_recal = open(fp_recal,'r')
    for line_recal in fd_recal:
        if line_recal[0] != '#':
            break

    d_lengths = loop_UG_out(l_vcfs_sorted,fd_recal,line_recal,minVQSLOD)

    plot_length_distribution(pop,coverage,d_lengths,ts_filter_level,)

    return


def plot_length_distribution(pop,coverage,d_lengths,ts_filter_level,):

    for key in d_lengths.keys():

        fn = 'lengths_%s%s_%4.1f_%s' %(pop,coverage,ts_filter_level,key)
        fd = open(fn,'w')
        fd.write(d_lengths[key])
        fd.close()

        gnuplot.histogram2(
            fn,
            x_step=1,
            x_max=12,
            xlabel='INDEL length',
            ylabel='INDEL count',
            color = 'blue',
            title = key,
            )

##    for key in d_lengths.keys():
##
##        x = d_lengths[key].split('\n')
##
##        r.png('lengths_%s_%s.png' %(pop,key))
##        r.hist(x, main='A histogram', xlab='x', col='lightblue')
##        r.dev_off()     

    return


def loop_UG_out(l_vcfs,fd_recal,line_recal,minVQSLOD):

    bool_exon = False

    if bool_exon == True:
        fn_exon = '/nfs/t149_influenza_exomes/working/Arrayid165_GRCh37_hs37d5_CTRplus_pos.txt'
        fd_exon = open(fn_exon)
        line_exon = fd_exon.readline()
        l_exon = line_exon.split()
        chrom_exon = l_exon[0]
        pos_exon = int(l_exon[1])

    d_lengths = {'PASS':{},'FAIL':{},}
    with fileinput.input(files=l_vcfs) as fd_vcf:
        for line_vcf in fd_vcf:
                    
            if fileinput.isfirstline():
                print(fileinput.filename())
                sys.stdout.flush()
            if line_vcf[0] == '#':
                continue
            l_vcf = line_vcf.split()
            ## skip SNPs
            if (
                l_vcf[3] in ['A','C','G','T',]
                and
                l_vcf[4] in [
                    ## UnifiedGenotyper
                    'A','C','G','T',
                    'A,C','A,G','A,T','C,G','C,T','G,T',
                    'A,C,G','A,C,T','A,G,T','C,G,T',
                    ## CombineVariants
                    '.',
                    'C,A','G,A','T,A','G,C','T,C','T,G',
                    'A,T,C','T,C,A','C,T,G','A,T,G','T,A,G','C,T,A','T,G,C','G,C,T',
                    ]
                ):
                continue
            ## CombineVariants
            if l_vcf[4] == '.': continue
            ## CombineVariants
            if len(l_vcf[3]) == 2 and len(l_vcf[4]) == 2:
                continue
            ## CombineVariants
            if len(l_vcf[3]) == len(l_vcf[4]) and ',' not in l_vcf[3] and ',' not in l_vcf[4]:
                print(2,l_vcf[3],l_vcf[4])
                continue
            ## CombineVariants
            bool_continue = True
            if ',' not in l_vcf[3]:
                for s in l_vcf[4].split(','):
                    if len(s) != len(l_vcf[3]):
                        bool_continue = False
                        break
            if bool_continue == True:
                print(3,l_vcf[3],l_vcf[4])
                continue

##            ## parse current VR recal line
##            l_recal = line_recal.split('\t')

            if bool_exon == True:
                pos_vcf = int(l_vcf[1])
                while True:
                    if chrom_exon == l_vcf[0]:
                        pass
                    else:
                        stoptmp
                    if pos_exon > pos_vcf:
                        bool_continue = True
                        break
                    elif pos_exon < pos_vcf:
                        pass
                    else:
                        bool_continue = False
                        break
                    line_exon = fd_exon.readline()
                    l_exon = line_exon.split()
                    chrom_exon = l_exon[0]
                    pos_exon = int(l_exon[1])
##
##            ## read next VR recal line
##            line_recal = fd_recal.readline()

            if bool_exon == True:
                if bool_continue == True:
                    continue
                line_exon = fd_exon.readline()
##
##            if l_recal[0] != l_vcf[0] or l_recal[1] != l_vcf[1]:
##                print('VRrecal')
##                print(l_recal)
##                print('UGvcf')
##                print(line_vcf)
##                stoptmp
##
##            ## parse INFO field/column
##            l_recal_INFO = l_recal[7].split(';')
##            ## parse VQSLOD value (would a regular expression be faster?)
##            VQSLOD = float(l_recal_INFO[1].split('=')[1])
##
##            if VQSLOD < minVQSLOD:
##                key = 'FAIL'
##            else:
##                key = 'PASS'

            key = 'PASS' ## tmp!!!

##            if VQSLOD == minVQSLOD:
##                print(VQSLOD)
##                print(minVQSLOD)
##                print(line_recal)
##                print(fileinput.filename())
##                stoptmp

            l_lengths1 = [len(s) for s in l_vcf[3].split(',')]
            l_lengths2 = [len(s) for s in l_vcf[4].split(',')]
##            if abs(len(l_lengths1)-len(l_lengths2)) >= 3:
##                print(l_vcf[3],l_vcf[4])
##                stop0
##            if len(l_lengths1) >= 1 and len(l_lengths2) >= 4:
##                print(l_vcf[3])
##                print(l_lengths2,l_vcf[4])
##                stop1
            if len(l_lengths1) >= 2:
                print(l_vcf[3])
                print(l_lengths2,l_vcf[4])
                stop2
            for len1 in l_lengths1:
                if len1 > 100: continue
                for len2 in l_lengths2:
                    if len2 > 100: continue
                    for i in range(9,len(l_vcf)):
                        l = l_vcf[i].split(':')[0].split('|')
                        if l_vcf[i][0] == '.':
                            print(l_vcf)
                            stop
                        ## DELETION
                        elif len1 != 1 and len2 == 1:
                            if l == ['1','1',]:
                                pass
                            elif l == ['0','0',] and len1 != 1 and len2 == 1:
                                if not abs(len2-len1) in d_lengths[key].keys():
                                    d_lengths[key][abs(len2-len1)] = 0
                                d_lengths[key][abs(len2-len1)] += 2
                            elif l == ['1','0',] and len1 != 1 and len2 == 1:
                                if not abs(len2-len1) in d_lengths[key].keys():
                                    d_lengths[key][abs(len2-len1)] = 0
                                d_lengths[key][abs(len2-len1)] += 1
                            elif l == ['0','1',] and len1 != 1 and len2 == 1:
                                if not abs(len2-len1) in d_lengths[key].keys():
                                    d_lengths[key][abs(len2-len1)] = 0
                                d_lengths[key][abs(len2-len1)] += 1
                        ## INSERTION
                        elif len1 == 1 and len2 != 1:
                            if l == ['0','0',]:
                                pass
                            elif l == ['1','0',] and len1 == 1 and len2 != 1:
                                if not abs(len2-len1) in d_lengths[key].keys():
                                    d_lengths[key][abs(len2-len1)] = 0
                                d_lengths[key][abs(len2-len1)] += 1
                            elif l == ['0','1',] and len1 == 1 and len2 != 1:
                                if not abs(len2-len1) in d_lengths[key].keys():
                                    d_lengths[key][abs(len2-len1)] = 0
                                d_lengths[key][abs(len2-len1)] += 1
                            elif l == ['1','1',] and len1 == 1 and len2 != 1:
                                if not abs(len2-len1) in d_lengths[key].keys():
                                    d_lengths[key][abs(len2-len1)] = 0
                                d_lengths[key][abs(len2-len1)] += 2
                        else:
##                            d_lengths[key] += '%i\n' %(abs(len2-len1))
##                        else:
                            print(len1,len2)
                            print(l)
                            print(l_vcf[i])
                            stop2
##                    d_lengths[key] += '%i\n' %(abs(len2-len1))

    if bool_exon == True:
        fd_exon.close()

    for x in range(max(d_lengths['PASS'].keys())):
        if x in d_lengths['PASS'].keys():
            print(x,d_lengths['PASS'][x])
    sys.exit()

    return d_lengths


def parse_minVQSLOD(ts_filter_level,fp_tranches):

    ts_filter_level
    fd = open(fp_tranches)
    lines = fd.readlines()
    fd.close()
    for line in lines:
        if line[0] == '#':
            continue
        l = line.split(',')
        if l[0] == 'targetTruthSensitivity':
            index = l.index('minVQSLod')
            continue
        targetTruthSensitivity = float(l[0])
        if targetTruthSensitivity == ts_filter_level:
            minVQSLOD = float(l[index])
            break

    return minVQSLOD


def alphanum_key(s):
    ## http://dave.st.germa.in/blog/2007/12/11/exception-handling-slow/
    NUM_RE = re.compile('([0-9]+)')
    return [ int(c) if c.isdigit() else c for c in NUM_RE.split(s) ]


def sort_nicely(l):
    ## http://nedbatchelder.com/blog/200712/human_sorting.html
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)
    return l


if __name__ == '__main__':
    main()
