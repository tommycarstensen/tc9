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

    l_fp1, l_fp2 = parse_args()

    loop(l_fp1,l_fp2)
    
    return


def generate_line_vcf(f):

    set_nt = set(['A','C','G','T',])

    for line in f:
        if line[0] == '#': continue
        l = line.strip().split('\t')
        if not l[3] in set_nt:
            continue
        if not l[4] in set_nt:
            continue
        if not l[6] == 'PASS':
            continue
        chrom = l[0]
        pos = int(l[1])

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


def gt_cnt(n_samples,d_intersection,l1,k,):

    for i in range(n_samples):
        gt = l1[9+i].split(':')[0]
        if gt == '0/0':
            continue
        elif gt == '0/1':
            d_intersection[i][k] += 1
        elif gt == '1/1':
            d_intersection[i][k] += 1
        elif gt == './.':
            continue
        else:
            print(gt)
            stop

    return


def generate_line_bgl(f):

    set_nt = set(['A','C','G','T',])

    for line in f:
        l = line.strip().split(' ')
        if not l[1] in set_nt:
            continue
        if not l[2] in set_nt:
            continue
        chrom,pos = l[0].split(':')
        pos = int(pos)

        yield chrom,pos,l


def return_monomorphic(l1):

    j = 0
    k = 0
    for i in range(3,len(l1),3):
        j += l1[i:i+3].index(max(l1[i:i+3]))
        k += 1
    print(l1[0],l1[1],j,k)
    if j == 0 or j == 2*k:
        bool_monomorphic = True
    else:
        bool_monomorphic = False

    return bool_monomorphic


def loop(l_fp1,l_fp2):

    l_chrom_seq = [str(chrom) for chrom in range(1,23)]    

##    with fileinput.input(files=l_fp1) as file1, fileinput.input(files=l_fp2) as file2:
    with fileinput.input(files=l_fp1) as file1, fileinput.input(files=l_fp2) as file2:

##        l_samples = return_samples(file1)
##        n_samples = len(l_samples)
##
##        d_intersection = {i:{
##            'intersect':0,
##            'complement':0,
##            } for i in range(n_samples)}

        n_intersection = 0
        n_complement = 0
        n_monomorphic = 0

        generator1 = generate_line_vcf
        generator1 = generate_line_bgl
##        chrom1,pos1,l1 = next(generate_line_vcf(file1))
        chrom1,pos1,l1 = next(generator1(file1))
        chrom2,pos2,l2 = next(generate_line_vcf(file2))

        prev = 0
        while True:

            if chrom1 != chrom2:
                if l_chrom_seq.index(chrom2) > l_chrom_seq.index(chrom1):
                    bool_monomorphic = return_monomorphic(l1)
                    if bool_monomorphic == True:
                        n_monomorphic += 1
                    else:
                        n_complement += 1
                    try:
                        chrom1,pos1,l1 = next(generator1(file1))
                    except StopIteration:
                        break
                else:
                    try:
                        chrom2,pos2,l2 = next(generate_line_vcf(file2))
                    except StopIteration:
                        break
                continue
            else:
                if pos1 == pos2:
                    if l1[2] == '.':
                        stop1
                    if False:
##                    if l1[3] != l2[3] or l1[4] != l2[4]:
####                        print(l1[3],l1[4],l2[3],l2[4],)
                        pass
                    else:
                        bool_monomorphic = return_monomorphic(l1)
                        if bool_monomorphic == True:
                            n_monomorphic += 1
                            pass
                        else:
                            n_intersection += 1
##                        gt_cnt(n_samples,d_intersection,l1,'intersect')
##                    j += 1
##                    if j > 1000000:
##                        break
##                    if n_intersection+n_complement > 10000:
##                        break
                    try:
                        chrom1,pos1,l1 = next(generator1(file1))
                    except StopIteration:
                        break
                    try:
                        chrom2,pos2,l2 = next(generate_line_vcf(file2))
                    except StopIteration:
                        while True:
                            try:
                                chrom1,pos1,l1 = next(generator1(file1))
                            except StopIteration:
                                break
                            bool_monomorphic = return_monomorphic(l1)
                            if bool_monomorphic == True:
                                n_monomorphic += 1
                                pass
                            else:
                                n_complement += 1
                        break
                    prev = 1
                elif pos2 > pos1:
##                    if l1[2] == '.':
                    if True:
                        bool_monomorphic = return_monomorphic(l1)
                        if bool_monomorphic == True:
                            n_monomorphic += 1
                            pass
                        else:
                            n_complement += 1
##                        gt_cnt(n_samples,d_intersection,l1,'complement',)
                    try:
                        chrom1,pos1,l1 = next(generator1(file1))
                    except StopIteration:
                        break
##                    if n_intersection+n_complement > 10000:
##                        break
                    prev = 2
##                    j += 1
##                    if j > 1000000:
##                        break
                else:
                    try:
                        chrom2,pos2,l2 = next(generate_line_vcf(file2))
                    except StopIteration:
                        while True:
                            try:
                                chrom1,pos1,l1 = next(generator1(file1))
                            except StopIteration:
                                break
                            bool_monomorphic = return_monomorphic(l1)
                            if bool_monomorphic == True:
                                n_monomorphic += 1
                                pass
                            else:
                                n_complement += 1
                        break
                    prev = 3

    print('intersection',n_intersection)
    print('complement',n_complement)
    print('monomorphic',n_monomorphic)
    stop
    print('sampleID dbSNP% Z')

    sumx = 0
    sumxx = 0
    for i in range(n_samples):
        intersection = d_intersection[i]['intersect']/(d_intersection[i]['intersect']+d_intersection[i]['complement'])
##        print(i,d_intersection[i]['intersect'],intersection)
        sumx += intersection
        sumxx += intersection*intersection
    SS = sumxx-(sumx**2)/n_samples
    var = SS/(n_samples-1) ## division with n if population, n-1 if sample
    stddev = math.sqrt(var)
    average = sumx/n_samples

    l = []
    for i in range(n_samples):
        intersection = d_intersection[i]['intersect']/(d_intersection[i]['intersect']+d_intersection[i]['complement'])
        Z = (intersection-average)/stddev
        l += [[Z,[round(intersection,3),Z]]]

    l.sort()
    for i in range(n_samples):
        print(l_samples[i], l[i][1])
    stop

    return


def parse_args():

##    parser = argparse.ArgumentParser()
##
##    parser.add_argument(
##        '--fp1',
##        dest='l_fp1',nargs='+',
##        help='List of file paths 1 in sequential order',
##        metavar='FILE',default=None,
##        required = True,
##        )
##
##    parser.add_argument(
##        '--fp2',
##        dest='l_fp2',nargs='+',
##        help='List of file paths 2 in sequential order (e.g. chrom1 chrom2...).',
##        metavar='FILE',default=None,
##        required = True,
##        )
##
##    parser.add_argument(
##        '--affix',
##        dest='affix',
##        help='File prefix/suffix',
##        metavar='FILE',default=None,
##        required = True,
##        )
##
####    ## http://docs.python.org/2/library/functions.html#vars
####    d_args = {}
####    print(parser.parse_args())
####    print(dict(parser.parse_args()))
####    for k,v in vars(parser.parse_args()).items():
####        d_args[k] = v
####        continue
##
##    extension1 = l_fp1[0][l_fp1[0].rindex('.')+1:]
##    extension2 = l_fp2[0][l_fp2[0].rindex('.')+1:]

##    l_fp1 = ['../pipeline/uganda4x/out_ApplyRecalibration/ApplyRecalibration.recalibrated.filtered.SNP.vcf']
    l_fp1 = ['../pipeline/%s/out_BEAGLE/%i.gprobs' %(sys.argv[-1],chrom) for chrom in range(1,23)]
##    l_fp1 = ['../pipeline/%s/out_BEAGLE/%i.gprobs' %(sys.argv[-1],22)]
    l_fp2 = ['/lustre/scratch107/projects/uganda/users/tc9/in_GATK/dbsnp_135.b37.vcf']

    return l_fp1, l_fp2


if __name__ == '__main__':

    main()
