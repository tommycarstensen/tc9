#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, January 2014

import os
import argparse
import contextlib
import fileinput
import sys


def main():

    '''assume that a vcf file with unphased genotypes is provided'''
    '''assume that samples are always in the same order if multiple vcf files'''
    '''assume that all samples are founders, when calculating LD and PCA'''
    '''optimize by carrying out ld pruning on fragments and merge binary files after pruning'''

    ## http://helix.nih.gov/Applications/README.eigenstrat

    ## For input in PACKEDPED format, snp file MUST be in genomewide order.
    ## For input in PACKEDPED format, genotype file MUST be in SNP-major order (the PLINK default: see PLINK documentation for details.)

    d_vars = init()

    vcf2bed(d_vars)

    bed2pca(d_vars['affix'], d_vars['EIGENSOFT'])

##    pca2png() ## let HGI decide what plotting software they want to use; gnuplot or matplotlib might not be their weapons of choice

    return


def vcf2bed(d_vars):

    '''
1) Estimate haplotype frequencies, 2) LD prune
(with no prior identification of founders) and
3) convert one or multiple vcf files to a single bed
'''

    ## split this long function into smaller sub routines...

    print('Estimation of haplotype frequencies and LD pruning')

    window = d_vars['window']
    step = d_vars['step']
    r2_max = d_vars['r2_max']
    itr_max = d_vars['itr_max']
    tol = d_vars['tol']

    ##
    ## vcf2fam
    ##
    with contextlib.ExitStack() as stack:
        fd_fam = stack.enter_context(
            open('EIGENSOFT/%s.pedind' % (d_vars['affix']), 'w'))
        fd_vcf = fileinput.FileInput(files=d_vars['vcf'])

        for line in fd_vcf:
            if line[0] != '#':
                break
            l = line.rstrip().split('\t')[9:]
        s = ''
        for sampleID in l:
            s += '%s %s 0 0 -9 -9\n' %(sampleID,sampleID)
        fd_fam.write(s)
    n = len(l)

    ##
    ## vcf2bed
    ##
    with contextlib.ExitStack() as stack:
        fd_bed = stack.enter_context(
            open('EIGENSOFT/%s.bed' % (d_vars['affix']), 'wb'))
        fd_bim = stack.enter_context(
            open('EIGENSOFT/%s.pedsnp' % (d_vars['affix']), 'w'))
        fd_vcf = fileinput.FileInput(files=d_vars['vcf'])

        magic_number = bytearray([108,27])
        mode = bytearray([1])
        fd_bed.write(magic_number+mode)

        l_l_dosages = []
        l_ID = []
        l_prune = []
        bool_break = False
        while True:

            ## parse vcf lines
            while len(l_l_dosages) < window:
                ## parse
                try:
                    chrom, pos, rsID, alleleA, alleleB, l_dosages = next(
                        parse_vcf(fd_vcf))
                except StopIteration:
                    bool_break = True
                    break
                if int(pos)%1000 == 0:
                    print(chrom,pos)
                    sys.stdout.flush()
                n_miss = l_dosages.count('.')
##                ## monomorphic or all heterozygous
##                if len(set(l_dosages) - set(['.'])) == 1:
##                    continue
                ## ignore sample call rates below 95%
                if n_miss/n > 0.05:
                    continue
                ## ignore MAF below 5%
                AF = (.5*l_dosages.count(1)+l_dosages.count(2))/(n-n_miss)
                if AF < 0.05 or AF > 0.95:
                    continue
                ## append
                l_l_dosages += [l_dosages]
                l_ID += [(chrom,pos,rsID,alleleA,alleleB)]

            ## calculate LD
            for i in range(len(l_l_dosages) - 1):
##                if i in l_prune:
##                    continue
                l_dosages1 = l_l_dosages[i]
                for j in range(i + 1, len(l_l_dosages)):
                    if j in l_prune:
                        continue
                    l_dosages2 = l_l_dosages[j]
                    ## different chromosome
                    if l_ID[i][0] != l_ID[j][0]:
                        break
                    rsq = esem_r(l_dosages1, l_dosages2, itr_max, tol)
                    if rsq >= r2_max:
                        l_prune += [j]
##                        print(i, j, rsq)
            
            ## append prune.in to bed/bim file
            barray = bytearray()
            s_bim = ''
            for i in range(step):
                b = ''
                if i in l_prune:
                    l_prune.remove(i)
                    continue
                genotype_data = ''
                chrom = l_ID[i][0]
                pos = l_ID[i][1]
                rsID = l_ID[i][2]
                alleleA = l_ID[i][3]
                alleleB = l_ID[i][4]
                s_bim += '%s\n' %('\t'.join((chrom,rsID,'0',pos,alleleA,alleleB)))
                for dosage in l_l_dosages[i]:
                    if dosage == 0:
                        b += '00'
                    elif dosage == 1:
                        b += '01'
                    elif dosage == 2:
                        b += '11'
                    else:
                        b += '10'
                    if len(b) == 8:
                        barray.append(int(b[::-1],2))
                        b = ''
                if b: ## i.e. if n % 4 != 0
                    barray.append(int(b[::-1].zfill(8),2))
            fd_bim.write(s_bim)
            fd_bed.write(bytes(barray))

            ## break after append
            if bool_break:
                break

            ## reset
            l_prune = [l_prune[i]-5 for i in range(len(l_prune))]
            l_l_dosages = l_l_dosages[step:]
            l_ID = l_ID[step:]

    return


def parse_vcf(fd_vcf):

    for line in fd_vcf:
        if line[0] == '#':
            continue
        l = line.rstrip().split('\t')
        chrom = l[0]
        pos = l[1]
        rsID = l[2]
        alleleA = l[3]
        alleleB = l[4]
        ## disallow INDELs and non-diallelic SNPs and monomorphic sites
        if not alleleB in ['A','C','G','T',]:
            continue
        l_dosages = []
        indexGT = l[8].split(':').index('GT')
        for i in range(9,len(l)):
            GT = l[i].split(':')[indexGT]
            if GT == './.':
                dosage = '.'
            else:
                dosage = sum(map(int, GT.split('/')))
            l_dosages += [dosage]
        yield chrom, pos, rsID, alleleA, alleleB, l_dosages

    return


def esem_r(Y, Z, itr_max, tol, h=[0.25, 0.25, 0.25, 0.25]):
    '''
# Use Excoffier-Slatkin EM algorithm to estimate r.
# Input:
#
# Y is vector of genotype values at 1st locus, coded as 0, 1, and 2
# to represent genotypes aa, aA, and AA.
#
# Z is the corresponding vector for 2nd locus.
#
# h is the initial vector of haplotype frequencies
#
# Function returns r.
'''
    h = esem(Y, Z, h, itr_max, tol)
    pA = h[0] + h[1]
    pB = h[0] + h[2]
    qA = 1.0 - pA
    qB = 1.0 - pB
    rsq = (h[0]*h[3] - h[1]*h[2])**2 / (pA*qA*pB*qB)
    return rsq


def esem(Y,Z, h, itr_max, tol):
    '''
# Excoffier-Slatkin EM algorithm for estimating haplotype frequencies.
# Input:
#
# Y is vector of genotype values at 1st locus, coded as 0, 1, and 2
# to represent genotypes aa, aA, and AA.
#
# Z is the corresponding vector for 2nd locus.
#
# h is the initial vector of haplotype frequencies
#
# Function returns h, a vector of 4 haplotype frequencies.
'''
    x = [[0,0,0],[0,0,0],[0,0,0]]
    for y,z in zip(Y,Z):
        if y == '.' or z == '.':
            continue
        try:
            x[y][z] += 1
        except:
            print(Z)
            print(y,z)
            stop
    for itr in range(itr_max):
        hh = esem_step(h, x)
        dh = 0.0
        for u,v in zip(h, hh):
            dh += abs(u-v)
        if dh <= tol:
            break
        h = hh
    if dh > tol:
        raise ConvergenceError
    return hh


def esem_step(h, x):

    '''
# Excoffier-Slatkin EM algorithm for estimating haplotype frequencies.
# This code implements the special case of the algorithm for two
# biallelic loci.  With two biallelic loci, there are 4 types of
# gamete, which I represent as follows:
#
#    AB  Ab  aB  ab
#     0   1   2   3
#
# Here A and a are the two alleles at the first locus and B and
# b are the alleles at the other locus.  The numbers below the
# gamete symbols are used below as indexes into arrays.
#
# Phenotypes at the 1st locus: AA, Aa, and aa are numbered 0, 1, and 2.
#
# Phenotypes at the 2nd locus: BB, Bb, and bb are numbered 0, 1, and 2.
#
# Input:
#
# h is a vector of 4 haplotype frequencies, indexed as shown above.
#
# x is a 3X3 matrix of phenotype counts.  The phenotypes are coded
# as explained above.  Thus, x[1][2] is the number of copies of
# the phenotype Aa/BB.
#
# n is the sample size and should equal the sum of x.
#
# Function returns a revised estimate of h after a single EM step.
'''

    # g is a 4X4 matrix of genotype frequencies. g[0][3]
    # is the frequency of the genotype that combines gamete 0 (AB)
    # with gamete 3 (ab).
    g = [[None,None,None,None],[None,None,None,None],
         [None,None,None,None],[None,None,None,None]]
    for i in range(4):
        g[i][i] = h[i]*h[i]
        for j in range(i):
            g[i][j] = 2*h[i]*h[j]

    # p is a 3X3 matrix of phenotype frequencies, recoded as
    # described for the input matrix x.
    p = [[None,None,None],[None,None,None],[None,None,None]]

    p[0][0] = g[0][0]
    p[0][1] = g[1][0]
    p[0][2] = g[1][1]

    p[1][0] = g[2][0]
    p[1][1] = g[3][0]+g[2][1]
    p[1][2] = g[3][1]

    p[2][0] = g[2][2]
    p[2][1] = g[3][2]
    p[2][2] = g[3][3]

    hh = [None, None, None, None]
    hh[0] = 2*x[0][0] + x[0][1] + x[1][0] + x[1][1]*g[3][0]/p[1][1]
    hh[1] = x[0][1] + 2*x[0][2] + x[1][1]*g[2][1]/p[1][1] + x[1][2]
    hh[2] = x[1][0] + x[1][1]*g[2][1]/p[1][1] + 2*x[2][0] + x[2][1]
    hh[3] = x[1][1]*g[3][0]/p[1][1] + x[1][2] + x[2][1] + 2*x[2][2]

    # haploid sample size
    n = float(sum(hh))

    # convert gamete counts counts to relative frequencies
    for i in range(4):
        hh[i] /= n

    return hh


def bed2pca(affix,path_EIGENSOFT):

    '''
run EIGENSOFT on a single input (bed or eigenstratgeno) file
(with no prior identification of founders)'''
    ## http://helix.nih.gov/Applications/README.eigenstrat

    ## For input in PACKEDPED format, snp file MUST be in genomewide order.
    ## For input in PACKEDPED format, genotype file MUST be in SNP-major order (the PLINK default: see PLINK documentation for details.)

    print('EIGENSOFT')

    s = ''
    s += 'genotypename: EIGENSOFT/%s.bed\n' %(affix)
    s += 'snpname: EIGENSOFT/%s.pedsnp\n' %(affix)
    s += 'indivname: EIGENSOFT/%s.pedind\n' %(affix)
    s += 'evecoutname: EIGENSOFT/%s.evec\n' %(affix)
    s += 'evaloutname: EIGENSOFT/%s.eval\n' %(affix)
    s += 'numoutlieriter: 0\n'
    s += 'outliersigmathresh: 6\n'
    s += 'outlieroutname: EIGENSOFT/%s.outlier\n' %(affix)
    s += 'qtmode: YES\n'
    s += 'snpweightoutname: EIGENSOFT/%s.snpweight\n' %(affix)

    with open('EIGENSOFT/%s.par' %(affix),'w') as fd:
        fd.write(s)

    os.system('%s -p EIGENSOFT/%s.par' %(path_EIGENSOFT,affix))

    return


def pca2png():

    '''plot EIGENSOFT PCA output'''

    return


def argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', '--in', required=True, nargs='+')
    parser.add_argument(
        '--affix', '--pca', '--out', required=True,
        help='affix of output and intermediate files',)
    parser.add_argument(
        '--EIGENSOFT', required=True,
        help='path to EIGENSOFT',
        default='/nfs/team149/Software/EIG4.2/bin/smartpca',
        )

    s = 'maximum change in sum of haplotype frequencies'
    s += ' between each EM step, faster if larger'
    parser.add_argument(
        '--tol', default=1e-4, type=float,
        help=s)
    parser.add_argument(
        '--itr_max', default=1000, type=int,
        help='maximum number of EM steps',)

    parser.add_argument(
        '--r2_max', default=0.2, type=float,
        help='r2 threshold (maximum) for LD pruning',)
    parser.add_argument(
        '--window', default=50, type=int,
        help='LD pruning window size, faster if smaller',)
    parser.add_argument(
        '--step', default=5, type=int,
        help='LD pruning step size, faster if larger',)

    namespace_args = parser.parse_args()
    d_vars = vars(namespace_args)

    if not os.path.isfile(d_vars['EIGENSOFT']):
        print('not found:', d_vars['EIGENSOFT'])
        sys.exit()

    if d_vars['window'] <= d_vars['step']:
        print('window size cannot be smaller than step size')
        sys.exit()

    if d_vars['r2_max'] >= 1 or d_vars['r2_max'] <= 0:
        print('r2 threshold outside of range 0-1')
        sys.exit()

    return d_vars


def mkdir():

    if not os.path.isdir('EIGENSOFT'):
        os.mkdir('EIGENSOFT')

    return


def init():

    d_vars = argparser()

    mkdir()

    return d_vars


if __name__ == '__main__':
    main()
