#!/bin/python3

## Tommy Carstensen (tc9)
## Wellcome Trust Sanger Institute, May 2013

## built-ins
import fileinput
import re
import math
import os
import sys


def main():

    path = '../pipeline/uganda4x/out_BEAGLE'
    affix = 'ethiopia'
    affix = sys.argv[-1]
    path = '../pipeline/%s4x/out_BEAGLE' %(affix)
    
    for chrom in range(1,23):
        if os.path.isfile('tped/%s%s.tped' %(affix,chrom)): continue ## tmp!!!
        fp_gprobs = os.path.join(path,'%s.gprobs' %(chrom))
        gprobs2tped(fp_gprobs,chrom,affix,)

    return


def gprobs2tped(fp_gprobs,chrom,affix,):

    set_nt = set(['A','C','G','T',])

    with open(fp_gprobs) as f_gprobs:

        ## skip header
        l_samples = f_gprobs.readline().strip().split()[3:-1:3]
        n_samples = len(l_samples)

        with open('tfam/%s%s.tfam' %(affix,chrom),'w') as f_tfam:

            lines_tfam = ['%s %s 0 0 -9 -9\n' %(ID,ID) for ID in l_samples]
            f_tfam.writelines(lines_tfam)

        with open('tped/%s%s.tped' %(affix,chrom),'w') as f_tped:

            for line_gprobs in f_gprobs:

                l_gprobs = line_gprobs.rstrip().split()
                chrom,pos = l_gprobs[0].split(':')
                if int(pos)%1000 == 0:
                    print(chrom,pos)
                alleleA = l_gprobs[1]
                alleleB = l_gprobs[2]
                if alleleA not in set_nt: continue
                if alleleB not in set_nt: continue

                line_tped = '%s %s:%s 0 %s' %(chrom,chrom,pos,pos)

                for i in range(3,3*(n_samples+1),3):
                    l_probs = [float(l_gprobs[i+j]) for j in range(3)]
                    if l_gprobs[i:i+3] == ['0.3333','0.3333','0.3333',]:
                        print(l_gprobs[i:i+3])
                        stop33
                    index = l_probs.index(max(l_probs))
                    if index == 0:
                        line_tped += ' %s %s' %(alleleA,alleleA)
                    elif index == 1:
                        line_tped += ' %s %s' %(alleleA,alleleB)
                    elif index == 2:
                        line_tped += ' %s %s' %(alleleB,alleleB)
                    else:
                        print(l_probs)
                        print(max(l_probs))
                        stop2

                line_tped += '\n'

                f_tped.write(line_tped)
    
    return


def generate_line_bgl(f):

    set_nt = set(['A','C','G','T',])

    for line in f:
        l = line.strip().split(' ')
        if not l[1] in set_nt: continue
        if not l[2] in set_nt: continue
        chrom,pos = l[0].split(':')
        pos = int(pos)

        yield chrom,pos,l


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


if __name__ == '__main__':

    main()
