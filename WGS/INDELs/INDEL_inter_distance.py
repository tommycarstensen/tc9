#!/bin/python3

## Tommy Carstensen (tc9)
## Wellcome Trust Sanger Institute, June 2013

import glob
import fileinput
from scipy import stats
import re
import matplotlib as mpl
## http://matplotlib.org/faq/howto_faq.html#matplotlib-in-a-web-application-server
mpl.use('Agg')
import matplotlib.pyplot as plt

def main():

    d_centromeres = get_centromeres()

    files_vcf = glob.glob('../pipeline/uganda4x/out_UnifiedGenotyper/*.vcf')
    files_vcf = sort_nicely(files_vcf)
    l_dist = []
    l_dist_even = []
    l_dist_odd = []
    with fileinput.input(files_vcf) as file_vcf:
        for chrom,pos,l_vcf in generate_line_vcf_INDEL(file_vcf):
            break
        chrom_prev = chrom
        pos_prev = pos
        for chrom,pos,l_vcf in generate_line_vcf_INDEL(file_vcf):
            if chrom == 'X': break
#            if chrom == 'Y': break
            if chrom != chrom_prev:
                pass
            elif pos_prev<d_centromeres[chrom]['acen']['max'] and pos>d_centromeres[chrom]['acen']['min']:
                pass
            else:
                bool_pass = False
                for posmin,posmax in d_centromeres[chrom]['gvar']:
                    if pos_prev<posmax and pos>posmin:
                        bool_pass = True
                        break
                if bool_pass == True:
                    pass
                else:
                    bool_even = False
                    bool_odd = False
                    alleleA = l_vcf[3]
                    if ',' in l_vcf[3]: stop
                    lenA = len(alleleA)
                    l_allelesB = l_vcf[4].split(',')
                    for alleleB in l_allelesB:
                        lenB = len(alleleB)
                        if (lenA-lenB)%2 == 0:
                            bool_even = True
                        else:
                            bool_odd = True
                    dist = pos-pos_prev
                    if dist < 2000:
                        l_dist += [dist]
                        if bool_odd == True:
                            l_dist_odd += [dist]
                        if bool_even == True:
                            l_dist_even += [dist]
                    else:
##                        print(dist,chrom_prev,chrom,pos_prev,pos)
                        pass
            chrom_prev = chrom
            pos_prev = pos
#            if len(l_dist) == 1000: break ## tmp!!!

    plt.hist(l_dist,200)
    plt.savefig('hist_INDEL_inter_distance.png')

    plt.hist(l_dist_odd,200)
    plt.savefig('hist_INDEL_inter_distance_odd.png')

    plt.hist(l_dist_even,200)
    plt.savefig('hist_INDEL_inter_distance_even.png')

    ## http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.describe.html#scipy.stats.describe
    size,tminmax,mean,var,skew,kurt = stats.describe(l_dist)
    print('size',size)
    print('tminmax',tminmax)
    print('mean',mean)
    print('var',var)
    print('skew',skew)
    print('kurt',kurt)

    return


def get_centromeres():

    d_centromeres = {str(chr):{'acen':[],'gvar':[],} for chr in list(range(1,23))+['X','Y',]}
    with open('cytoBand.txt') as f:
        for line in f:
            l = line.rstrip().split()
            k = l[4]
            if k not in ['acen','gvar',]: continue
            chr = l[0][3:]
            pos1 = int(l[1])
            pos2 = int(l[2])
            if k == 'acen':
                d_centromeres[chr][k] += [pos1,pos2]
            else:
                d_centromeres[chr][k] += [[pos1,pos2,]]
    for chr in d_centromeres.keys():
        for k in ['acen']:
            posmin = min(d_centromeres[chr][k])
            posmax = max(d_centromeres[chr][k])
            d_centromeres[chr][k] = {'min':posmin,'max':posmax}

    return d_centromeres


def generate_line_vcf_INDEL(f):

    set_nt = set(['A','C','G','T',])

    for line in f:
        if line[0] == '#': continue
        l = line.rstrip().split('\t')
        ## skip SNPs
        if (
            l[3] in ['A','C','G','T',]
            and
            l[4] not in [
                ## UnifiedGenotyper
                'A','C','G','T',
                'A,C','A,G','A,T','C,G','C,T','G,T',
                'A,C,G','A,C,T','A,G,T','C,G,T',
    ##                ## CombineVariants
    ##                'C,A','G,A','T,A','G,C','T,C','T,G',
##                'A,T,C','T,C,A','C,T,G','A,T,G','T,A,G','C,T,A','T,G,C','G,C,T',
    ##                ## UG GENOTYPE_GIVEN_ALLELES
    ##                'C,G,A','C,A,T',
    ##                'G,A,C','G,T,A','G,T,C','G,C,A',
    ##                'A,G,C',
    ##                'T,A,C','T,G,A',
                ]
            ):
            continue
        chrom = l[0]
        pos = int(l[1])

        yield chrom,pos,l


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
