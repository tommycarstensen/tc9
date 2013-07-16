#!/bin/python3

## Tommy Carstensen (tc9)
## Wellcome Trust Sanger Institute, May 2013

## built-ins
import fileinput
import re
import math
##sys.path.append('/nfs/users/nfs_t/tc9/github/sandbox')
##import gnuplot

def main():

    (
        l_fp1, l_fp2,
        fp_samples1, fp_samples2,
        fp_samples1to2,
        format2) = parse_args()

    l_indexes1, l_indexes2, d_samples1to2 = samples2indexes(
        l_fp1,l_fp2,fp_samples1,fp_samples2,fp_samples1to2,format2,)

    loop(l_fp1,l_fp2,d_samples1to2,l_indexes1,l_indexes2,format2,)
    
    return


def samples2indexes(l_fp1,l_fp2,fp_samples1,fp_samples2,fp_samples1to2,format2,):

    l_samples1,l_samples2 = parse_samples(fp_samples1,fp_samples2,l_fp1,l_fp2,)

    l_indexes1 = []
    l_indexes2 = []

    with open(fp_samples1to2) as f:
        d_samples1to2 = {}
        for line in f:
            k,v = line.strip().split()
            d_samples1to2[k] = v

    for index1 in range(len(l_samples1)):
        sample1 = l_samples1[index1]
        try:
            sample2 = d_samples1to2[sample1]
        except:
            sample2 = d_samples1to2[sample1]
            continue
        try:
            if format2 == 'BEAGLE':
                ## marker alleleA alleleB
                index2 = 3+3*l_samples2.index(sample2)
            elif format2 == 'IMPUTE2':
                ## --- rsID position alleleA alleleB
                index2 = 5+3*l_samples2.index(sample2)
            else:
                stop
        except ValueError:
            continue
        except NameError:
            print(format2)
            stop
        except:
            stop
        l_indexes2 += [index2]
        l_indexes1 += [4+2*index1]

    return l_indexes1, l_indexes2, d_samples1to2


def parse_samples(fp_samples1,fp_samples2,l_fp1,l_fp2,):

    d_samples = {}
    l = [[fp_samples1,l_fp1,],[fp_samples2,l_fp2,],]
    for i in range(2):
        fp_samples = l[i][0]
        l_fp = l[i][1]
        if fp_samples == None:
            with open(l_fp[0]) as f:
                l_samples = f.readline().strip().split()[3:-1:3]
            print(l_samples)
            stop
        elif fp_samples[-5:] == '.tfam':
            l_samples = parse_samples_tfam(fp_samples)
        elif fp_samples[-4:] == '.bgl':
            with open(fp_samples) as f:
                l_samples = f.readline().strip().split()[3:-1:3]
        else:
            print(fp_samples)
            stop
        d_samples[i] = l_samples

    return d_samples[0],d_samples[1]


def parse_samples_tfam(fp_samples,):

    ## convert "plate_well_sample" format to "sample" format
    ## i.e. remove info about plate and well
    keyword = re.compile(r'(\d\d\d\d\d\d_[A-H]\d\d_)(.+\d\d\d\d\d\d\d)')
    l_samples = []
    with open(fp_samples) as lines:
        for line in lines:
            sampleID = line.split()[1]
            match = result = keyword.search(sampleID)
            if match:
                l_samples += [match.group(2)]
            else:
                l_samples += [sampleID]

    return l_samples


def generate_line_bgl(f):

    set_nt = set(['A','C','G','T',])

    for line in f:
        l = line.strip().split(' ')
        if not l[1] in set_nt: continue
        if not l[2] in set_nt: continue
        chrom,pos = l[0].split(':')
        pos = int(pos)

        yield chrom,pos,l


def generate_line_tped(f):

    for line in f:
        l = line.strip().split(' ')
        chrom = l[0]
        pos = int(l[3])

        yield chrom,pos,l


def correlation(
    n_samples,l1,l2,d_stats,l_indexes1,l_indexes2,i_alleleA2,i_alleleB2,):

    alleleA2 = l2[i_alleleA2]
    alleleB2 = l2[i_alleleB2]

    for i in range(n_samples):
        index1 = l_indexes1[i]
        alleleA1 = l1[index1]
        alleleB1 = l1[index1+1]
        index2 = l_indexes2[i]
        l_probs = l2[index2:index2+3]
        x = None
        y = None
        bool_discordant = False
        if alleleA1 == '0' and alleleB1 == '0':
##            d_stats[i]['missing'] += 1
            continue
        ## insertion/deletion
        elif alleleA2 == '0' and alleleB2 == '1':
            stop1
            continue
        ## homozygous, AA
        elif float(l_probs[0]) > 0.9:
            if (
                alleleA2 == alleleA1
                and
                alleleA2 == alleleB1
                ):
                x = 0
                pass
            else:
                ## homozygous
                if alleleA1 == alleleB1:
                    x = 2
                ## heterozygous
                else:
                    x = 1
                bool_discordant = True
        ## heterozygous, Aa or aA
        elif float(l_probs[1]) > 0.9:
            if any([
                all([
                    alleleA2 == alleleA1,
                    alleleB2 == alleleB1,
                    ]),
                all([
                    alleleB2 == alleleA1,
                    alleleA2 == alleleB1,
                    ]),
                ]):
                x = 1
                pass
            else:
                ## homozygous 1
                if (
                    alleleA1 == alleleA2
                    and
                    alleleB1 == alleleA2
                    ):
                    x = 0
                ## homozygous 2
                else:
                    x = 2
                bool_discordant = True
        ## homozygous, aa
        elif float(l_probs[2]) > 0.9:
            if (
                alleleB2 == alleleA1
                and
                alleleB2 == alleleB1
                ):
                x = 2
                pass
            else:
                ## homozygous
                if alleleA1 == alleleB1:
                    x = 0
                ## heterozygous
                else:
                    x = 1
                bool_discordant = True
        else:
            x = None

        if x != None:
            ## x is chip "dosage", y is sequence/imputation dosage
            y = 0
##            y += 0*float(l_probs_bgl[0])
            y += 1*float(l_probs[1])
            y += 2*float(l_probs[2])
            d_stats[i]['sumxy'] += x*y
            d_stats[i]['sumxx'] += x*x
            d_stats[i]['sumyy'] += y*y
            d_stats[i]['sumx'] += x
            d_stats[i]['sumy'] += y
            d_stats[i]['n'] += 1
            if bool_discordant == True:
                d_stats[i]['discordant'] += 1
        else:
            d_stats[i]['missing'] += 1

    return


def loop(l_fp1,l_fp2,d_samples1to2,l_indexes1,l_indexes2,format2,):

    n_samples = len(l_indexes1)## = len(l_indexes2)
    d_stats = {i:{
        'n':0,
        'sumx':0,
        'sumy':0,
        'sumxy':0,
        'sumxx':0,
        'sumyy':0,
        'discordant':0,
        'missing':0,
        } for i in range(n_samples)}

    extension1 = 'tped'
    extension2 = 'bgl'
    d_func = {
        'tped':generate_line_tped,
        'bgl':generate_line_bgl,
        }
    func1 = d_func[extension1]
    func2 = d_func[extension2]

    if format2 == 'BEAGLE':
        i_alleleA2 = 1
        i_alleleB2 = 2
    elif format2 == 'IMPUTE':
        i_alleleA2 = 3
        i_alleleB2 = 4
    else:
        print(format2)

    n_intersection = 0

    with fileinput.input(files=l_fp1) as file1, fileinput.input(files=l_fp2) as file2:

        chrom1,pos1,l1 = next(func1(file1))
        chrom2,pos2,l2 = next(func2(file2))

        while True:

            if chrom1 != chrom2:
                remembertodointegersorindexes
                print(chrom1,chrom2)
                stop
            else:
                if pos1 == pos2:
                    correlation(
                        n_samples,l1,l2,d_stats,l_indexes1,l_indexes2,
                        i_alleleA2,i_alleleB2,)
                    n_intersection += 1
                    if n_intersection == 10000: break ## tmp!!!
                    try:
                        chrom1,pos1,l1 = next(func1(file1))
                        chrom2,pos2,l2 = next(func2(file2))
                    except StopIteration:
                        break
                elif pos2 > pos1:
                    try:
                        chrom1,pos1,l1 = next(func1(file1))
                    except StopIteration:
                        break
                else:
                    try:
                        chrom2,pos2,l2 = next(func2(file2))
                    except StopIteration:
                        break

    print(n_intersection)
    print('sampleID r2 Z')
    sumx = 0
    sumxx = 0
    l_r2 = []
    for i in range(n_samples):
        n = d_stats[i]['n']
        nom = (d_stats[i]['sumxy']-d_stats[i]['sumx']*d_stats[i]['sumy']/n)
        den_sq = (
            (d_stats[i]['sumxx']-d_stats[i]['sumx']**2/n)*
            (d_stats[i]['sumyy']-d_stats[i]['sumy']**2/n))
        den = math.sqrt(abs(den_sq))
        r = nom/den
        r2 = r**2
##        print(i,r,r2)
        l_r2 += [r2]
        sumx += r2
        sumxx += r2*r2
    SS = sumxx-(sumx**2)/n_samples
    var = SS/(n_samples-1) ## division with n if population, n-1 if sample
    stddev = math.sqrt(var)
    average = sumx/n_samples
    l = []
    l_samples = parse_samples_tfam('../stats/uganda.tfam')
    for i in range(n_samples):
        r2 = l_r2[i]
        Z = (r2-average)/stddev
##        print(i,r,r2,Z)
        l += [[Z,[r2,Z]]]
    l.sort()
    for i in range(n_samples):
        print(l_samples[i], l[i][1])
    stop

    ##
    ## correlation / concordance
    ##
    fd = open('%s.correlation' %(affix),'w')
    fd.writelines(l_corr)
    fd.close()
    
    fd = open('%s.concordance' %(affix),'w')
    fd.writelines(l_conc)
    fd.close()

    cmd = 'cat %s.correlation' %(affix)
    cmd += ''' | awk '{if($2!="None") {sum+=$2;n++}} END{print sum/n}' '''
    r2_avg = float(os.popen(cmd).read())

    ##
    ##
    ##
    s = ''
    s += 'fp1: %s\n' %(fp1)
    s += 'fp2: %s\n' %(fp2)
    s += 'fp1: %i variants/SNPs\n' %(cnt_1_not_2+cnt_1_and_2)
    s += 'fp2: %i variants/SNPs\n' %(cnt_2_not_1+cnt_1_and_2)
    s += 'intersection: %i variants/SNPs\n' %(cnt_1_and_2)
    s += 'n_samples: %i\n' %(n_samples)
    if cnt1 != cnt_1_not_2+cnt_1_and_2:
        print(cnt1, cnt_1_not_2+cnt_1_and_2)
        print(cnt_1_not_2, cnt_1_and_2)
        stoptmp1
    if cnt2 != cnt_2_not_1+cnt_1_and_2:
        print(cnt2, cnt_2_not_1+cnt_1_and_2)
        print(cnt_2_not_1, cnt_1_and_2)
        stoptmp2
    if chrom1 != '22' or chrom2 not in['22','Y',]:
        print(chrom1, chrom2, pos1, pos2)
        stoptmp3
    ## discordant
    s += 'cnt_discordant_all_SNPs: %i\n' %(cnt_discordant_all_SNPs)
    ## missing
    s += 'cnt_missing_all_SNPs: %i\n' %(cnt_missing_all_SNPs)
    s += 'concordance: %.4f\n' %(
        1-(float(cnt_discordant_all_SNPs)/float(n_samples*cnt_1_and_2)))
    s += 'correlation: %.4f\n' %(r2_avg)            
    print(s)
    fd = open('%s.stats' %(affix),'a')
    fd.write(s)
    fd.close()
    for fn_SNPs,l_SNPs in [
        ['%s.1_not_2.SNPs' %(affix),l_1_not_2,],
        ['%s.1_and_2.SNPs' %(affix),l_1_and_2,],
        ]:
        s_SNPs = '\n'.join(l_SNPs)+'\n'
        fd = open(fn_SNPs,'w')
        fd.write(s_SNPs)
        fd.close()

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
##    parser.add_argument(
##        '--fp_samples1',
##        dest='fp_samples1',
##        help='File with sample IDs in sequence corresponding to sequence of genotype probabilities',
##        metavar='FILE',default=None,
##        required = False,
##        )
##
##    parser.add_argument(
##        '--fp_samples2',
##        dest='fp_samples2',
##        help='File with sample IDs in sequence corresponding to sequence of genotype probabilities',
##        metavar='FILE',default=None,
##        required = False,
##        )
##
##    parser.add_argument(
##        '--fp_samples1to2',
##        dest='fp_samples1to2',
##        help='Translation of sample IDs from file1 to file2',
##        metavar='FILE',default=None,
##        required = False,
##        )
##
##    parser.add_argument(
##        '--format1',
##        dest='format1',
##        metavar='STRING',default=None,
##        required = False,
##        )
##
##    parser.add_argument(
##        '--format2',
##        dest='format2',
##        metavar='STRING',default=None,
##        required = False,
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

    l_fp1 = ['../stats/uganda.tped']
##    fp2 = '../pipeline/uganda4x/out_ApplyRecalibration/ApplyRecalibration.recalibrated.filtered.SNP.vcf'
    l_fp2 = ['../pipeline/uganda4x/out_ProduceBeagleInput/ProduceBeagleInput.%i.bgl' %(chrom) for chrom in range(1,23)]
    l_fp2 = ['../pipeline/uganda4x/out_BEAGLE/%i.gprobs' %(chrom) for chrom in range(1,23)]
    fp_samples1 = '../stats/uganda.tfam'
    fp_samples2 = '../pipeline/uganda4x/out_ProduceBeagleInput/ProduceBeagleInput.bgl'
    fp_samples1to2 = '../stats/uganda.dic'
    format2 = 'BEAGLE'

    return l_fp1, l_fp2, fp_samples1, fp_samples2, fp_samples1to2, format2


if __name__ == '__main__':

    main()
