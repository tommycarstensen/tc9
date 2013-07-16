#!/bin/python3

## Tommy Carstensen (tc9)
## Wellcome Trust Sanger Institute, May 2013

## built-ins
import fileinput
import re
import math
import sys
##sys.path.append('/nfs/users/nfs_t/tc9/github/sandbox')
##import gnuplot

def main():

    l_fp = ['../pipeline/uganda4x/out_ApplyRecalibration/ApplyRecalibration.recalibrated.filtered.SNP.vcf']
##    l_fp = ['../pipeline/uganda4x/out_BEAGLE/%i.gprobs' %(chrom) for chrom in range(1,23)]

    loop(l_fp)
    
    return


def return_samples(f):

    for line in f:

        ## VCF
        match = re.match('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t',line)
        if match:
            l_samples = line[len(match.group(0)):].rstrip().split('\t')
            break

##        ## BEAGLE
##        l_samples = line.rstrip().split()[3:-1:3]
##        break

    return l_samples


def loop(l_fp):

    set_nt = set(['A','C','G','T',])

    with fileinput.input(files=l_fp) as f:

        l_samples = return_samples(f)
        n_samples = len(l_samples)
        d_TiTv = {i:{
            'Ti':0,
            'Tv':0,
            } for i in range(n_samples)}

        d_allele2TiTv = {
            'GA':'Ti',
            'CT':'Ti',
            'AG':'Ti',
            'TC':'Ti',
            'CG':'Tv',
            'GC':'Tv',
            'GT':'Tv',
            'TG':'Tv',
            'TA':'Tv',
            'AT':'Tv',
            'CA':'Tv',
            'AC':'Tv',
            }

        j = 0
        for line in f:

            ## VCF
            l = line.rstrip().split('\t')
            if not l[3] in set_nt: continue ## deletion
            if not l[4] in set_nt: continue ## insertion or triallelic
            if l[6] != 'PASS': continue
            for i in range(n_samples):
                TiTv = d_allele2TiTv[l[3]+l[4]]
                gt = l[9+i].split(':')[0]
                if gt == '0/0':
                    continue
                elif gt == './.':
                    continue
                elif gt == '1/1':
                    d_TiTv[i][TiTv] += 2
                elif gt == '0/1':
                    d_TiTv[i][TiTv] += 1
                else:
                    print(gt)
                    stop

##            ## BEAGLE
##            l = line.rstrip().split(' ')
##            if not l[1] in set_nt: continue
##            if not l[2] in set_nt: continue
##            for i in range(n_samples):
##                TiTv = d_allele2TiTv[l[1]+l[2]]
##                if float(l[3+3*i]) > 0.9:
##                    continue
##                elif float(l[4+3*i]) > 0.9:
##                    d_TiTv[i][TiTv] += 1
##                elif float(l[5+3*i]) > 0.9:
##                    d_TiTv[i][TiTv] += 2
##                else:
##                    continue

            j += 1
            if j % 1000 == 0: print(j)
            if j > 1000000:
##            if j > 10000:
                break

    sumx = 0
    sumxx = 0
    for i in range(n_samples):
        TiTv = d_TiTv[i]['Ti']/d_TiTv[i]['Tv']
        sumx += TiTv
        sumxx += TiTv*TiTv
    SS = sumxx-(sumx**2)/n_samples
    var = SS/(n_samples-1) ## division with n if population, n-1 if sample
    stddev = math.sqrt(var)
    average = sumx/n_samples
    l = []
    for i in range(n_samples):
        TiTv = d_TiTv[i]['Ti']/d_TiTv[i]['Tv']
        Z = round((TiTv-average)/stddev,2)
        l += [(Z,i,d_TiTv[i]['Ti'],d_TiTv[i]['Tv'],TiTv,Z)]
    l.sort()
    for i in range(n_samples):
        print(l_samples[l[i][1]],l[i][2:])

    return


if __name__ == '__main__':

    main()
