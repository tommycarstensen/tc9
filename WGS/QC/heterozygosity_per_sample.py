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

##        ## BEAGLE
##        l_samples = line.rstrip().split()[3:-1:3]
##        break

        ## VCF
        match = re.match('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t',line)
        if match:
            break

    l_samples = line[len(match.group(0)):].rstrip().split('\t')

    return l_samples


def loop(l_fp):

    set_nt = set(['A','C','G','T',])

    with fileinput.input(files=l_fp) as f:

        l_samples = return_samples(f)
        n_samples = len(l_samples)
        d_het = {i:{
            'hom':0,
            'het':0,
            } for i in range(n_samples)}

        j = 0
        for line in f:

##            ## BEAGLE
##            l = line.rstrip().split(' ')
##            if not l[1] in set_nt: continue
##            if not l[2] in set_nt: continue
##            for i in range(n_samples):
##                if float(l[3+3*i]) > 0.9:
##                    d_het[i]['hom'] += 1
##                elif float(l[4+3*i]) > 0.9:
##                    d_het[i]['het'] += 1
##                elif float(l[5+3*i]) > 0.9:
##                    d_het[i]['hom'] += 1
##                else:
##                    continue

            ## VCF
            l = line.rstrip().split('\t')
            if not l[3] in set_nt: continue ## deletion
            if not l[4] in set_nt: continue ## insertion or triallelic
            if l[6] != 'PASS': continue
            for i in range(n_samples):
                gt = l[9+i].split(':')[0]
                if gt == '0/0':
                    d_het[i]['hom'] += 1
                elif gt == './.':
                    continue
                elif gt == '1/1':
                    d_het[i]['hom'] += 1
                elif gt == '0/1':
                    d_het[i]['het'] += 1
                else:
                    print(gt)
                    stop

            j += 1
            if j > 1000000:
                break

    sumx = 0
    sumxx = 0
    for i in range(n_samples):
        het = d_het[i]['het']/(d_het[i]['hom']+d_het[i]['het'])
        sumx += het
        sumxx += het*het
    SS = sumxx-(sumx**2)/n_samples
    var = SS/(n_samples-1) ## division with n if population, n-1 if sample
    stddev = math.sqrt(var)
    average = sumx/n_samples
    l = []
    for i in range(n_samples):
        het = round(d_het[i]['het']/(d_het[i]['hom']+d_het[i]['het']),4)
        Z = round((het-average)/stddev,2)
        l += [(Z,i,d_het[i]['het'],d_het[i]['hom'],het,Z)]
    l.sort()
    for i in range(n_samples):
        print(l_samples[l[i][1]],l[i][2:])

    return


if __name__ == '__main__':

    main()
