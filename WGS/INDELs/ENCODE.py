#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, 2013

from ftplib import FTP
import os
import glob
import re
import fileinput
import sys
##import rpy
##from rpy import r
sys.path.append('/nfs/users/nfs_t/tc9/github/sandbox')
import gnuplot

def main():

#    parse_ftp()

    pop = 'uganda'
    coverage = '4x'
    dn = 'DNAse-seq'
    dn = sys.argv[-1]
    print(pop,coverage,dn,)

    l_vcfs = glob.glob('../pipeline/%s%s/out_UnifiedGenotyper/*.vcf' %(pop,coverage,))
    l_vcfs_sorted = sort_nicely(l_vcfs)
    l_vcfs_sorted = ['../pipeline/uganda4x/out_UnifiedGenotyper/UnifiedGenotyper.2.2.vcf']
####    l_vcfs_sorted = l_vcfs_sorted[:10]
    sep = '/'
#    l_vcfs_exome = glob.glob('/nfs/t149_1kg/phase1_v3/ALL.chr*.exome.1000GApr12.vcf')
#    l_vcfs_exome_sorted = sort_nicely(l_vcfs_exome)
#    l_vcfs = glob.glob('/nfs/t149_1kg/phase1_v3/ALL.chr*.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz')
#    l_vcfs_sorted = sort_nicely(l_vcfs)
#    sep = '|'

##    gnuplot.contour_plot(
##        path_dat='INDEL_%s.dat' %(dn),
##        bool_remove=False,
##        xlabel='dist_m_i_n (1000bp)',ylabel='INDEL length',
##        bool_log=False,
##        x1=0,x2=50,y1=1,y2=20,
##        )
##    stop

    d_lengths = loop_UG_out(l_vcfs_sorted,dn,sep)

    plot_contour(dn,d_lengths)

    s = ''
    for k1 in d_lengths.keys():
        for k2,v in d_lengths.items():
            s += '%s %s %s\n' %(k1,k2,v)
    fd = open('%s.dict','w')
    fd.write(s)
    fd.close()

    stats(d_lengths,dn,)

########    plot_length_distribution(pop,coverage,d_lengths,)

    return


def stats(d_lengths,dn,):

    for bool_skip in [False,True,]:

        even = []
        odd = []
        for dist_min in d_lengths.keys():
            for len_diff in d_lengths[dist_min].keys():
                if bool_skip == True:
                    if len_diff == 1: continue
                if len_diff % 2 == 0:
                    even += d_lengths[dist_min][len_diff]*[dist_min]
                else:
                    odd += d_lengths[dist_min][len_diff]*[dist_min]

        import scipy
        from scipy import stats

        u,p = stats.mannwhitneyu(even,odd)
        fd = open('stats','a')
        fd.write('mannwhitneyu u %s p %s %s %s\n' %(u,p,dn,bool_skip))
        fd.close()
        
        z,p = stats.ranksums(even,odd)
        fd = open('stats','a')
        fd.write('ranksums z %s p %s %s %s\n' %(z,p,dn,bool_skip))
        fd.close()

        average_even = sum(even)/len(even)
        average_odd = sum(odd)/len(odd)
        fd = open('stats','a')
        fd.write('average even %s odd %s %s %s\n' %(average_even,average_odd,dn,bool_skip))
        fd.close()
    
    return


def plot_contour(dn,d_lengths):

    l_keys = list(d_lengths.keys())
##    for k in l_keys:
##        if k > 200:
##            del d_lengths[k]
    gnuplot.contour_plot(
        d_dat=d_lengths,fileprefix='INDEL_%s' %(dn),bool_remove=False,
        xlabel='dist_m_i_n (1000bp)',ylabel='INDEL length',
        bool_log=True,
        x2=50,
        y2=50,
        )

    return


def parse_ftp():

    d_paths = {
        'DNAse-seq':'pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/fdrPeaks',
        'FAIRE':'pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/faire_fseq_peaks',
        'Histone':'pub/databases/ensembl/encode/integration_data_jan2011/byDataType/peaks/jan2011/histone_macs/optimal/hub',
        }
    
    ftp = FTP('ftp.ebi.ac.uk')
    ftp.login()
    for k,path in d_paths.items():
##        print(ftp.dir(path))
        l_files = ftp.nlst(path)
##        print(ftp.mlsd(path))
        if not os.path.isdir(k):
            os.mkdir(k)
        for f in l_files:
            if f[-3:] != '.bb':
                continue
            path_out = os.path.join(k,os.path.basename(f))
            if os.path.isfile(path_out):
                continue
            ftp.retrbinary(
                'RETR %s' %(f),
                open(path_out,'wb').write,
                )
    ftp.quit()
    for dn in d_paths.keys():
        l_files = os.listdir(dn)
        for f in l_files:
            if f[-3:] != '.bb':
                continue
            path_in = os.path.join(dn,f)
            path_out = os.path.join(dn,f[:-3]+'.bed')
            if os.path.isfile(path_out):
                continue
            ## assume bigBedToBed binary is in cwd
            cmd = './bigBedToBed %s %s' %(path_in,path_out)
            os.system(cmd)

    return


def plot_length_distribution(pop,coverage,d_lengths,):

    for key in d_lengths.keys():

        fn = 'lengths_%s%s_%s' %(pop,coverage,key)
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


def generate_line_INDEL(fd_vcf):

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
                ## UG GENOTYPE_GIVEN_ALLELES
                'C,G,A','C,A,T',
                'G,A,C','G,T,A','G,T,C','G,C,A',
                'A,G,C',
                'T,A,C','T,G,A',
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

        yield line_vcf,l_vcf


def merge_intervals(dn_ENCODE):

    l_files = glob.glob(os.path.join(dn_ENCODE,'*.bed'))
##    l_files = l_files[:1] ## tmp!!!

    d_io = {}
    for f in l_files:
        d_io[f] = open(f)

    d_intervals_merged = {}
    l_intervals_next = []
    l_chroms = [str(chrom) for chrom in range(1,23)]+['X','Y',]
    l_chroms.sort()
##    for chrom in [1]+list(range(10,20))+[2]+list(range(20,23))+list(range(3,10)):
    for chrom in l_chroms:
        l_intervals = l_intervals_next
        l_intervals_next = []
        for f,io in d_io.items():
            print(chrom,f)
            for line in io:
                l = line.split()
                pos1 = int(l[1])
                pos2 = int(l[2])
                if chrom != l[0][3:]:
                    l_intervals_next += [[pos1,pos2,]]
                    break
                l_intervals += [[pos1,pos2,]]
        l_intervals.sort()

        d_intervals_merged[chrom] = [l_intervals[0]]
        for interval in l_intervals[1:]:
##            if interval[0] != d_intervals_merged[chrom][-1][0]:
##                d_intervals_merged[chrom] += [interval]
##            else:
            if d_intervals_merged[chrom][-1][1] < interval[0]:
                d_intervals_merged[chrom] += [interval]
            else:
                d_intervals_merged[chrom][-1][1] = interval[1]

##        break ## tmp!!!

    for f in l_files:
        d_io[f].close()

    del l_intervals

    return d_intervals_merged


def generate_chrom_pos(fd_vcf):

    for line_vcf in fd_vcf:

        if line_vcf[0] == '#':
            continue

        l_vcf = line_vcf[:40].split()
        chrom = l_vcf[0]
        pos = int(l_vcf[1])
        try:
            x = l_vcf[2]
        except:
            print(line_vcf)
            stoptmp

        yield chrom, pos


def skip_exome(generator_chrom_pos_skip,d_chroms,chrom,pos,chrom_skip,pos_skip,):

    bool_continue = False
    while True:
        if chrom_skip != chrom:
            if d_chroms[chrom_skip] < d_chroms[chrom]:
                chrom_skip, pos_skip = next(generator_chrom_pos_skip)
                continue
            else:
                bool_continue = False
                break
        if pos_skip < pos:
            chrom_skip, pos_skip = next(generator_chrom_pos_skip)
            continue
        elif pos_skip > pos:
            bool_continue = False
            break
        else:
            bool_continue = True
            break

    return bool_continue, chrom_skip, pos_skip


def return_dist_min(chrom,pos,d_intervals,d_interval_index,):

    try:
        interval1 = d_intervals[chrom][d_interval_index[chrom]]
        try:
            interval2 = d_intervals[chrom][d_interval_index[chrom]+1]
            while pos > interval2[1]:
                interval1 = interval2
                d_interval_index[chrom] += 1
                try:
                    interval2 = d_intervals[chrom][d_interval_index[chrom]+1]
##                    print('a',pos,interval2)
                except IndexError:
                    interval2 = interval1
                    break
        except IndexError:
            interval2 = interval1
    except IndexError:
        interval1 = interval2
    dist_min = round(abs(min(pos-interval1[1],interval2[0]-pos)),-3)

    return dist_min


def loop_UG_out(l_vcfs,dn_ENCODE,sep,):

    d_intervals = merge_intervals(dn_ENCODE)

    d_lengths = {} ## k=distfromhist or k1=length,k2=distfromhist
    d_interval_index = {str(chrom):0 for chrom in range(1,23)}
    chrom_skip = '1'
    pos_skip = 0
    d_chroms = {str(chrom):chrom for chrom in range(1,23)}
    d_chroms['X'] = 23
    d_chroms['Y'] = 24
    with fileinput.input(files=l_vcfs,openhook=fileinput.hook_compressed) as fd_vcf:

##    for vcf in l_vcfs:
##      with gzip.open(vcf) as fd_vcf:

        for line_vcf,l_vcf in generate_line_INDEL(fd_vcf):

            if l_vcf[0] == 'X': break

            chrom = l_vcf[0]
            pos = int(l_vcf[1])

            if pos % 1000 == 0:
                print(l_vcf[:2])

            dist_min = return_dist_min(chrom,pos,d_intervals,d_interval_index,)
            key = int(dist_min/1000)
            try:
                d_lengths[key]
            except KeyError:
                d_lengths[key] = {}
            if dist_min < 1000:
                print(dist_min,chrom,pos,l_vcf[3],l_vcf[4],)
            if dist_min > 2000000:
                print(l_vcf[:7])
                print(fileinput.filename())
                print(key,dist_min,pos,interval2)
                print(d_intervals[l_vcf[0]][-3])
                print(d_intervals[l_vcf[0]][-2])
                print(d_intervals[l_vcf[0]][-1])
                stoptmp

            l_lengths1 = [len(s) for s in l_vcf[3].split(',')]
            l_lengths2 = [len(s) for s in l_vcf[4].split(',')]
            l_lengths = l_lengths1+l_lengths2
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

            if len(l_lengths1) > 1 and len(l_lengths2) > 1:
                print(l_vcf[3])
                print(l_lengths2,l_vcf[4])
                stopstop

            for i in range(9,len(l_vcf)):
                if l_vcf[i] == './.':
                    continue
                if l_vcf[i] == '.|.':
                    print(l_vcf,i)
                    stop
                l = l_vcf[i].split(':')
                l_GT = [int(i) for i in l[0].split(sep)]
                len_ref = l_lengths[0]
                len1 = l_lengths[l_GT[0]]
                len2 = l_lengths[l_GT[1]]
                for len_indel in [len1,len2,]:
                    len_diff = abs(len_indel-len_ref)
                    if len_diff > 100:
                        continue
                    if len_diff == 0:
                        continue
                    try:
                        d_lengths[key][len_diff] += 1
                    except KeyError:
                        d_lengths[key][len_diff] = 1

##    for x in range(max(d_lengths['PASS'].keys())):
##        if x in d_lengths['PASS'].keys():
##            print(x,d_lengths['PASS'][x])
##    sys.exit()

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
