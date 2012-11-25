#!/software/bin/python-2.7.3

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import os
import sys
sys.path.append('/nfs/users/nfs_t/tc9/github/tc9/misc')
import gnuplot
sys.path.append('/nfs/users/nfs_t/tc9/github/tc9/math')
import statistics

##
## this script counts the number of passed variants
## and plots them as a function of AF and DP
##

def main():

    ## bsub -J"count$count" -o count$count.out -e count$count.err python ~/github/ms23/analysis/count_passed_variants.py $count
    ## bsub  -M500000 -R'select[mem>500] rusage[mem=500]' -J"count$count" -o count$count.out -e count$count.err python ~/github/ms23/analysis/count_passed_variants.py $count

##    fd = open('AF3.txt','r')
##    lines = fd.readlines()
##    fd.close()
##    l_AF3 = [float(s) for s in lines]
##    fd = open('AF13.txt','r')
##    lines = fd.readlines()
##    fd.close()
##    l_AF13 = [float(s) for s in lines]
##    fd = open('AF23.txt','r')
##    lines = fd.readlines()
##    fd.close()
##    l_AF23 = [float(s) for s in lines]
##    fd = open('AF123.txt','r')
##    lines = fd.readlines()
##    fd.close()
##    l_AF123 = [float(s) for s in lines]
##
##    import collections
##    AF13_multiset = collections.Counter(l_AF13)
##    AF23_multiset = collections.Counter(l_AF23)
##    AF123_multiset = collections.Counter(l_AF123)
##    AF3_multiset = collections.Counter(l_AF3)
##    print len(l_AF3)
##    l_AF3 = list((AF3_multiset - AF13_multiset).elements())
##    print len(l_AF3), len(l_AF13)
##    AF3_multiset = collections.Counter(l_AF3)
##    l_AF3 = list((AF3_multiset - AF23_multiset).elements())
##    print len(l_AF3), len(l_AF23)
##    AF3_multiset = collections.Counter(l_AF3)
##    l_AF3 = list((AF3_multiset - AF123_multiset).elements())
##    print len(l_AF3), len(l_AF123)
##    stop
##
####    print 'a'
####    for x in l_AF123:
####        l_AF3.remove(x)
####    print 'b'
####    for x in l_AF12:
####        l_AF3.remove(x)
####    print 'c'
####    for x in l_AF13:
####        l_AF3.remove(x)

    gnuplot.histogram2(
        'AF3',title='MAF distribution - 2.5M chip array',l_data=l_AF3,
        x_min=0,x_max=0.5,x_step=0.01,tic_step=0.05,xlabel='MAF',ylabel='SNP count',)
    gnuplot.histogram2(
        'AF13',title='MAF distribution - 2.5M chip array and HGI SNPs',l_data=l_AF13,
        x_min=0,x_max=0.5,x_step=0.01,tic_step=0.05,xlabel='MAF',ylabel='SNP count',)
    gnuplot.histogram2(
        'AF23',title='MAF distribution - 2.5M chip array and GATK SNPs',l_data=l_AF23,
        x_min=0,x_max=0.5,x_step=0.01,tic_step=0.05,xlabel='MAF',ylabel='SNP count',)
    gnuplot.histogram2(
        'AF123',title='MAF distribution - 2.5M chip array and HGI and GATK SNPs',l_data=l_AF123,
        x_min=0,x_max=0.5,x_step=0.01,tic_step=0.05,xlabel='MAF',ylabel='SNP count',)

    stop

    if sys.argv[-1] == '5':
        ## 5) compare mp15 steps
        fp1 = 'out_mp15/beagle/03.merged.vcf'
        fp2 = 'out_mp15/impute2/$CHROMOSOME'
        fp3 = '../omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.map'
        count_unique_and_intersect_vqsr(
            fp1,fp2,fp3,
            'mp15_BEAGLE_vs_IMPUTE2',
            bool_combined1 = True,
            bool_combined2 = False,
            )
        return

    elif sys.argv[-1] == '4':
        ## 4) compare mp15 steps
        fp1 = 'out_mp15/vqsr/$CHROMOSOME.vqsr.filt.vcf'
        fp2 = 'out_mp15/beagle/03.merged.vcf'
        fp3 = '../omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.map'
        count_unique_and_intersect_vqsr(
            fp1,fp2,fp3,
            'mp15_VQSR_vs_BEAGLE',
            bool_combined1 = False,
            bool_combined2 = True,
            )
        return

    elif sys.argv[-1] == '3':
        ## 3) compare mp15 steps
##        fp1 = 'out_mp15/pre-vqsr/$CHROMOSOME.vcf'
        fp1 = 'out_mp15/pre-vqsr/$CHROMOSOME.vqsr.vcf'
        fp2 = 'out_mp15/vqsr/$CHROMOSOME.vqsr.filt.vcf'
        fp3 = '../omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.map'
        count_unique_and_intersect_vqsr(
            fp1,fp2,fp3,
            'mp15_pre-VQSR_vs_post-VQSR',
            bool_combined2 = False,
            bool_combined1 = False,
            bool_ignore_FILTER1 = True,
            )
        return

    elif sys.argv[-1] == '2':
        ## 2) compare tc9 steps
        fp1 = 'out_GATK/join/CombineVariants.vcf'
        fp2 = 'out_GATK/join/ApplyRecalibration.recalibrated.filtered.vcf'
        fp3 = '../omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.map'
        count_unique_and_intersect_vqsr(
            fp1,fp2,fp3,
            'tc9_pre-VQSR_vs_post-VQSR',
            bool_combined1 = True,
            bool_combined2 = True,
            )
        return

    elif sys.argv[-1] == '1':
        ## 1) compare post-VQSR
        fp2 = 'out_GATK/join/ApplyRecalibration.recalibrated.filtered.vcf'
        fp1 = 'out_mp15/vqsr/$CHROMOSOME.vqsr.filt.vcf'
        fp3 = '../omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.map'
        count_unique_and_intersect_vqsr(
            fp1,fp2,fp3,
            'post-VQSR_tc9_vs_mp15',
            bool_combined2 = True,
            bool_combined1 = False,
            )
        return

    elif sys.argv[-1] == '0':
        ## 1) compare pre-VQSR
##        fp1 = 'out_mp15/pre-vqsr/$CHROMOSOME.vcf'
        fp1 = 'out_mp15/pre-vqsr/$CHROMOSOME.vqsr.vcf'
        fp2 = 'out_GATK/join/CombineVariants.vcf'
        fp3 = '../omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.map'
        count_unique_and_intersect_vqsr(
            fp1,fp2,fp3,
            'pre-VQSR_tc9_vs_mp15',
            bool_combined2 = True,
            bool_combined1 = False,
            bool_ignore_FILTER1 = True,
            )
        return

##    t1 = time.time()
##    singlevcf_vs_multiplevcfs()
##    t2 = time.time()
##    print 'time', t2-t1

    count_unique_and_intersect_impute2()
    stop

    count_and_plot()

    return


def count_unique_and_intersect_impute2():

    fp1 = 'out_mp15/impute2/04.merged.vcf'
    fp2 = 'out_IMPUTE2/join/IMPUTE2.$CHROMOSOME.gen'

    fd1 = open(fp1,'r')
    fd2 = open(fp2.replace('$CHROMOSOME','1',),'r')

    for line in fd2:
        POS2 = int(line.split()[2])
        break

    count1 = 0
    count2 = 1
    count_intersect = 0

    chromosome = None
    for line1 in fd1:

        l1, CHROM1, POS1, REF1, ALT1, FILTER1, bool_continue = parse_line_vcf(line1)
        if bool_continue == True:
            continue

        count1 += 1

        if CHROM1 != chromosome:
            print chromosome
            for line in fd2:
                count2 += 1
            chromosome = CHROM1
            fd2 = open(fp2.replace('$CHROMOSOME',chromosome,),'r')
            for line in fd2:
                POS2 = int(line.split()[2])
                count2 += 1
                break

        if POS1 == POS2:
            count_intersect += 1
            for line in fd2:
                POS2 = int(line.split()[2])
                count2 += 1
                break
        ## loop over lines1
        elif POS1 < POS2:
            continue
        else:
            for line in fd2:
                POS2 = int(line.split()[2])
                count2 += 1
                if POS1 == POS2:
                    count_intersect += 1
                    for line in fd2:
                        POS2 = int(line.split()[2])
                        count2 += 1
                        break
                ## loop over lines2
                elif POS2 < POS1:
                    continue
                ## loop over lines1
                else:
                    break
                break

    print count1
    print count2
    print count_intersect
    stop

    return


def count_unique_and_intersect_vqsr(
    fp_vcf_individual,fp_vcf_combined,fp_map,
    suffix,
    bool_combined1 = False,
    bool_combined2 = True,
    bool_ignore_FILTER1 = False,
    ):

    '''this function assumes that markers in the VCFs are sorted'''

    ## set file paths
    fp2_template = fp2 = fp_vcf_combined
    fp1_template = fp1 = fp_vcf_individual
    fp3 = fp_map

    print fp1
    print fp2
    print fp3

    ## set list of chromosomes (genotype array only contains autosomal SNPs)
    l_chromosomes = [str(i) for i in range(1,22+1,)]+['X','Y',]

    ## set initial chromosome
    chromosome1 = l_chromosomes[0]
    chromosome2 = l_chromosomes[0]
##    chromosome1 = chromosome2 = '22'

    fp1 = fp1.replace('$CHROMOSOME',chromosome1,)
    fp2 = fp2.replace('$CHROMOSOME',chromosome2,)

    ## set booleans before loop
    bool_read1 = True
    bool_read2 = True
    bool_read3 = True
    bool_EOF1 = False
    bool_EOF2 = False
    bool_EOF3 = False

    ## set counters before loop
    count_intersect12 = 0
    count_intersect13 = 0
    count_intersect23 = 0
    count1 = 0
    count2 = 0
    count3 = 0
    count_intersect123 = 0

    l_AF13 = []
    l_AF23 = []
    l_AF3 = []
    l_AF123 = []

    ## open files before loop
    fd1 = open(fp1,'r')
    fd2 = open(fp2,'r')
    fd3 = open(fp3,'r')
    fd3b = open('../omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.tped','r')
    
####    fd1.seek(1500000000000)
####    s = fd1.readline()
####    s = fd1.readline()
##    fd2.seek(65000000000)
##    s = fd2.readline()
##    s = fd2.readline()
##    fd3.seek(50000000)
##    s = fd3.readline()
##    s = fd3.readline()

    i = 0

    while True:

        i += 1

        if i % 10000 == 0:
            print CHROM1, CHROM2, CHROM3, POS1, POS2, POS3, '|', bool_read1, bool_read2, bool_read3, '|', bool_EOF1, bool_EOF2, bool_EOF3, i/10000
####            if CHROM1 != None and CHROM2 != None and CHROM3 != None:
####                if abs(l_chromosomes.index(CHROM1)-l_chromosomes.index(CHROM2)) > 1:
####                    print CHROM1, CHROM2
####                    stop1
####                if abs(l_chromosomes.index(CHROM2)-l_chromosomes.index(CHROM3)) > 1:
####                    print CHROM2, CHROM3
####                    stop2
##        if bool_read2 == True and count2-(count_intersect12+count_intersect23+count_intersect123) > 0:
##            print
##            print 2, fp2
##            print CHROM2, '***', POS2, '***'
##            print count2-(count_intersect12+count_intersect23+count_intersect123)
##            print count2, count_intersect12, count_intersect23, count_intersect123
##            print fp1, chromosome1
##            print fp2, chromosome2
##            stop2
####        if bool_read1 == True and count1-(count_intersect12+count_intersect23+count_intersect123) > 0:
####            print
####            print 1, fp1
####            print CHROM1, '***', POS1, '***'
####            print count1-(count_intersect12+count_intersect13+count_intersect123)
####            print count1, count_intersect12, count_intersect13, count_intersect123
####            print fp1, chromosome1
####            print fp2, chromosome2
####            stop1

        if bool_read1 == True:
            if bool_combined1 == True:
                CHROM1, POS1, count1, bool_EOF1 = loop_single_vcf(fd1,count1,)
            else:
                (
                    CHROM1, POS1, count1, fd1, bool_EOF1,
                    chromosome1,
                    ) = loop_multiple_vcf(
                        fd1,count1,fp1_template,chromosome1,l_chromosomes,
                        bool_ignore_FILTER=bool_ignore_FILTER1,
                        )
            bool_read1 = False

        if bool_read2 == True:
            if bool_combined2 == True:
                CHROM2, POS2, count2, bool_EOF2 = loop_single_vcf(fd2,count2,)
            else:
                (
                    CHROM2, POS2, count2, fd2, bool_EOF2,
                    chromosome2,
                    ) = loop_multiple_vcf(
                        fd2,count2,fp2_template,chromosome2,l_chromosomes,
                        )
            bool_read2 = False

        if bool_read3 == True:
            line3 = fd3.readline()
            if line3 == '':
                bool_EOF3 = True
                CHROM3 = None
                POS3 = None
            else:
                l3 = line3.split()
                CHROM3 = l3[0]
                POS3 = int(l3[3])
                count3 += 1

                ## tmp AF
                line3b = fd3b.readline()
                l = line3b.split()[4:]
                AF = l.count(l[0])/184.
                if AF > 0.5:
                    AF = 1-AF
                l_AF3 += [AF]

            bool_read3 = False


        if bool_EOF1 == True and bool_EOF2 == True and bool_EOF3 == True:
            break

##        print POS1, POS2, POS3, CHROM1, CHROM2, CHROM3
##        stop

        ## doing nested if statements is the fastest method of comparison
        ## looping over lines simultaneously to avoid reading all markers into memory

        ##
        ## triple intersection
        ##
##        if POS1 == POS2 == POS3 and CHROM1 == CHROM2 == CHROM3:
        if POS1 == POS2 == POS3:
            if POS1 % 1000 == 0:
                print CHROM1, '%6i' %(POS1/1000), '|',
                print '%8i' %(count_intersect12), '%8i' %(count_intersect13), '%8i' %(count_intersect23), '|',
                print '%8i' %(count_intersect123), '|',
                print '%8i' %(count1), '%8i' %(count2), '%8i' %(count3)
##            count_intersection12 += 1
##            count_intersection13 += 1
##            count_intersection23 += 1
            count_intersect123 += 1
            l_AF123 += [AF]
            bool_read1 = True
            bool_read2 = True
            bool_read3 = True
        ##
        ## double intersection
        ##
        else:
            ## it is faster to do nesting of logical statements
            ## when comparing long integers
            if POS1 == POS2 != None:
                if bool_EOF3 == True:
                    count_intersect12 += 1
                    bool_read1 = True
                    bool_read2 = True
                elif CHROM1 == CHROM3 and POS1 < POS3:
                    count_intersect12 += 1
                    bool_read1 = True
                    bool_read2 = True
                elif CHROM1 == CHROM3:
                    bool_read3 = True
                else:
                    if l_chromosomes.index(CHROM3) < l_chromosomes.index(CHROM1):
                        bool_read3 = True
                    else:
                        count_intersect12 += 1
                        bool_read1 = True
                        bool_read2 = True
            elif POS1 == POS3 != None:
                if bool_EOF2 == True:
                    count_intersect13 += 1
                    l_AF13 += [AF]
                    bool_read1 = True
                    bool_read3 = True
                elif CHROM1 == CHROM2 and POS1 < POS2:
                    count_intersect13 += 1
                    l_AF13 += [AF]
                    bool_read1 = True
                    bool_read3 = True
                elif CHROM1 == CHROM2:
                    bool_read2 = True
                else:
                    print CHROM1, CHROM2
                    stop2
            elif POS2 == POS3 != None:
                if bool_EOF1 == True:
                    count_intersect23 += 1
                    l_AF23 += [AF]
                    bool_read2 = True
                    bool_read3 = True
                elif CHROM1 == CHROM2 and POS2 < POS1:
                    count_intersect23 += 1
                    l_AF23 += [AF]
                    bool_read2 = True
                    bool_read3 = True
                elif CHROM1 == CHROM2:
                    bool_read1 = True
                else:
                    if l_chromosomes.index(CHROM1) < l_chromosomes.index(CHROM2):
                        bool_read1 = True
                    else:
                        count_intersect23 += 1
                        bool_read2 = True
                        bool_read3 = True
                        stop3tmp_wegethereornot
            ##
            ## no intersection
            ##
            else:
                ## different chromosomes
                if (
                    (CHROM1 != CHROM2)
                    or
                    (bool_EOF3 == False and CHROM2 != CHROM3)
                    ):
                    l_indexes= []
                    if bool_EOF1 == False:
                        index1 = l_chromosomes.index(CHROM1)
                        l_indexes += [index1]
                    if bool_EOF2 == False:
                        index2 = l_chromosomes.index(CHROM2)
                        l_indexes += [index2]
                    if bool_EOF3 == False:
                        index3 = l_chromosomes.index(CHROM3)
                        l_indexes += [index3]
                    min_index = min(l_indexes)
                    if bool_EOF1 == False and index1 == min_index:
                        bool_read1 = True
                    if bool_EOF2 == False and index2 == min_index:
                        bool_read2 = True
                    if bool_EOF3 == False and index3 == min_index:
                        bool_read3 = True
                ## same chromosome
                else:
                    ## either read 1 or 3
                    if CHROM1 == CHROM2 and POS1 < POS2:
                        if bool_EOF3 == True or POS1 < POS3:
                            bool_read1 = True
                        else:
                            bool_read3 = True
                    elif bool_EOF3 == True:
                        bool_read2 = True
                    ## either read 2 or 3
                    elif CHROM2 == CHROM3:
                        if bool_EOF3 == True or POS2 < POS3:
                            bool_read2 = True
                        else:
                            bool_read3 = True
                    else:
                        print CHROM1, CHROM2, CHROM3, POS1, POS2, POS3
                        stop

    print count1
    print count2
    print count3
    print count_intersect12
    print count_intersect13
    print count_intersect23
    print count_intersect123
    print
    print count1-count_intersect12-count_intersect13-count_intersect123
    print count2-count_intersect12-count_intersect23-count_intersect123
    print count3-count_intersect13-count_intersect23-count_intersect123
    print
    print fp1
    print fp2
    print fp3

    print 'AF3', sum(l_AF3)/len(l_AF3)
    print 'AF13', sum(l_AF13)/len(l_AF13)
    print 'AF23', sum(l_AF23)/len(l_AF23)
    print 'AF123', sum(l_AF123)/len(l_AF123)

    gnuplot.histogram2(
        'AF3',title='MAF distribution - 2.5M chip array',l_data=l_AF3,
        x_min=0,x_max=0.5,x_step=0.01,tic_step=0.05,xlabel='MAF',ylabel='SNP count',)
    gnuplot.histogram2(
        'AF13',title='MAF distribution - 2.5M chip array and HGI SNPs',l_data=l_AF13,
        x_min=0,x_max=0.5,x_step=0.01,tic_step=0.05,xlabel='MAF',ylabel='SNP count',)
    gnuplot.histogram2(
        'AF23',title='MAF distribution - 2.5M chip array and GATK SNPs',l_data=l_AF23,
        x_min=0,x_max=0.5,x_step=0.01,tic_step=0.05,xlabel='MAF',ylabel='SNP count',)
    gnuplot.histogram2(
        'AF123',title='MAF distribution - 2.5M chip array and HGI and GATK SNPs',l_data=l_AF123,
        x_min=0,x_max=0.5,x_step=0.01,tic_step=0.05,xlabel='MAF',ylabel='SNP count',)

    l_AF3 = [str(f) for f in l_AF3]
    l_AF13 = [str(f) for f in l_AF13]
    l_AF23 = [str(f) for f in l_AF23]
    l_AF123 = [str(f) for f in l_AF123]

    fd = open('AF3.txt','w')
    fd.write('\n'.join(l_AF3))
    fd.close()
    fd = open('AF13.txt','w')
    fd.write('\n'.join(l_AF13))
    fd.close()
    fd = open('AF23.txt','w')
    fd.write('\n'.join(l_AF23))
    fd.close()
    fd = open('AF123.txt','w')
    fd.write('\n'.join(l_AF123))
    fd.close()

    gnuplot.venn3(
        i1 = count1-count_intersect12-count_intersect13-count_intersect123,
        i2 = count2-count_intersect12-count_intersect23-count_intersect123,
        i3 = count3-count_intersect13-count_intersect23-count_intersect123,
        i4 = count_intersect12,
        i5 = count_intersect13,
        i6 = count_intersect23,
        i7 = count_intersect123,
        text1 = '%s' %(fp1),
        text2 = '%s' %(fp2),
        text3 = '%s' %(fp3),
        suffix = suffix,
        )

    return


def find_smaller(l_chromosomes, CHROM1, CHROM2, POS1, POS2,):

    if l_chromosomes.index(CHROM1) < l_chromosomes.index(CHROM2):
        smaller = 2
    else:
        smaller = 1

    return smaller


def loop_single_vcf(fd2,count2):

    CHROM2 = None
    POS2 = None
    bool_EOF2 = False

    while True:
        line2 = fd2.readline()
##            for line2 in fd2:
        ## EOF
        if line2=='':
            bool_EOF2 = True
            break
        (
            l2, CHROM2, POS2, REF2, ALT2, FILTER2, bool_continue,
            ) = parse_line_vcf(line2,)
        if bool_continue == True:
            continue
        else:
            count2 += 1
            break

    return CHROM2, POS2, count2, bool_EOF2


def loop_multiple_vcf(
    fd1,count1,fp1_template,chromosome,l_chromosomes,
    bool_ignore_FILTER=False,
    ):

    CHROM1 = None
    POS1 = None
    bool_EOF1 = False

    while True:
        line1 = fd1.readline()
##            for line1 in fd1:
        ## EOF
        if line1=='':
            index = l_chromosomes.index(chromosome)+1
            if index == len(l_chromosomes):
                bool_EOF1 = True
                break
            else:
                chromosome = l_chromosomes[index]
            fp1 = fp1_template.replace(
                '$CHROMOSOME',chromosome,
                )
            if not os.path.isfile(fp1):
                bool_EOF1 = True
                break
            ## tmp write time and file being opened to tmp file
            ## instead of waiting for stdout to be written to file
            ## when using the cluster
            import time
            fd = open('tmp.txt','a')
            fd.write('%s %s\n' %(str(time.gmtime()),str(fp1)))
            fd.close()
            print fp1
            fd1 = open(fp1,'r')
            continue
        l1, CHROM1, POS1, REF1, ALT1, FILTER1, bool_continue = parse_line_vcf(
            line1, bool_ignore_FILTER=bool_ignore_FILTER,
            )
        if bool_continue == True:
            continue
        else:
            count1 += 1
            bool_read1 = False
            break

    return CHROM1, POS1, count1, fd1, bool_EOF1, chromosome


def singlevcf_vs_multiplevcfs():

    '''count variant call differences between mpileup and unifiedgenotyper'''

    fp2 = 'out_GATK/join/ApplyRecalibration.recalibrated.filtered.vcf'
    fp1 = 'out_mp15/vqsr/$CHROMOSOME.vqsr.filt.vcf'

    fd2 = open(fp2,'r')

    ##
    ## parse first line of second file
    ##
    for line2 in fd2:
        l2, CHROM2, POS2, REF2, ALT2, FILTER2, bool_continue = parse_line_vcf(line2)
        if bool_continue == True:
            continue
        break

##    count1 = 0
##    count2 = 0
##    count_intersect = 0
##    count_unique1 = 0
##    count_unique2 = 0
##    for chromosome in [str(i) for i in range(1,22+1,)]+['X','Y',]:
##        print chromosome
##        l_pos2 = [POS2]
##        for line2 in fd2:
##            l2, CHROM2, POS2, REF2, ALT2, FILTER2, bool_continue = parse_line_vcf(line2)
##            if bool_continue == True: continue
##            if CHROM2 != chromosome: break
##            l_pos2 += [POS2]
##        l_pos1 = []
##        fd1 = open('out_mp15/vqsr/%s.vqsr.filt.vcf' %(chromosome),'r')    
##        for line2 in fd1:
##            l2, CHROM2, POS2, REF2, ALT2, FILTER2, bool_continue = parse_line_vcf(line2)
##            if bool_continue == True: continue
##            l_pos1 += [POS2]
##        count1 += len(l_pos1)
##        count2 += len(l_pos2)
##        set1 = set(l_pos1)
##        del l_pos1
##        set2 = set(l_pos2)
##        del l_pos2
##        count_intersect += len(set1&set2)
##        count_unique1 += len(set1-set2)
##        count_unique2 += len(set2-set1)
##        del set1
##        del set2
##    print count1
##    print count2
##    print count_intersect
##    print count_unique1, count1-count_intersect
##    print count_unique2, count2-count_intersect
##    stop

    d_QUAL = {}
    count_intersect = 0
    count1 = 0
    count2 = 1
    l_chromosomes = [str(i) for i in range(1,22+1,)]+['X','Y',]
    for chromosome in l_chromosomes:
        d_QUAL[chromosome] = []

        fd1 = open('out_mp15/vqsr/%s.vqsr.filt.vcf' %(chromosome),'r')
        for line1 in fd1:
            l1, CHROM1, POS1, REF1, ALT1, FILTER1, bool_continue = parse_line_vcf(line1)
            if bool_continue == True:
                continue

            count1 += 1

##            if CHROM1 == '2' or CHROM2 == '2':
##                break

            ## end of vcf2
            if CHROM2 == None:
                continue

            if CHROM1 != CHROM2:
                ## loop over lines2
                if l_chromosomes.index(CHROM1) > l_chromosomes.index(CHROM2):
                    for line2 in fd2:
                        l2, CHROM2, POS2, REF2, ALT2, FILTER2, bool_continue = parse_line_vcf(line2)
                        if bool_continue == True:
                            continue
                        count2 += 1
                        if CHROM1 != CHROM2:
                            continue
    ##                    bool_chromosome_diff = True
                        break
                ## loop over lines1
                else:
    ##                bool_chromosome_diff = True
                    continue

            if POS1 == POS2:
                if l1[5] != '999' and l2[5] != '999':
                    d_QUAL[CHROM1] += ['%s %s\n' %(l1[5],l2[5],)]
                count_intersect += 1
                for line2 in fd2:
                    l2, CHROM2, POS2, REF2, ALT2, FILTER2, bool_continue = parse_line_vcf(line2)
                    if bool_continue == True:
                        continue
                    count2 += 1
                    break
                continue

            ## loop over lines1
            elif POS2 > POS1:
                if CHROM1 != CHROM2:
                    print 'b', CHROM1, CHROM2
                    stop
                continue

            ## loop over lines2
##            elif POS1 > POS2:
            else:
                for line2 in fd2:
                    l2, CHROM2, POS2, REF2, ALT2, FILTER2, bool_continue = parse_line_vcf(line2)
                    if bool_continue == True:
                        continue
                    count2 += 1
                    if CHROM1 != CHROM2:
                        print 'c', CHROM1, CHROM2
                        stop
                    else:
                        ## loop over lines1
                        if POS2 > POS1:
                            break
                        elif POS1 == POS2:
                            if l1[5] != '999' and l2[5] != '999':
                                d_QUAL[CHROM1] += ['%s %s\n' %(l1[5],l2[5],)]
                            count_intersect += 1
                            for line2 in fd2:
                                l2, CHROM2, POS2, REF2, ALT2, FILTER2, bool_continue = parse_line_vcf(line2)
                                if bool_continue == True:
                                    continue
                                count2 += 1
                                break
                            break
                        else:
                            continue

        fd1.close()

    fd2.close()

    print count_intersect
    print count1
    print count2

    fd = open('QUAL.gnuplotdata','w')
    for chromosome in d_QUAL.keys():
        print chromosome, len(d_QUAL[chromosome])
        fd.writelines(d_QUAL[chromosome])
    fd.close()
    gnuplot.scatter_plot_2d(
        'QUAL',regression=True,
        xlabel='QUAL UnifiedGenotyper',
        ylabel='QUAL mpileup',
        )

    l_QUAL1 = []
    l_QUAL2 = []
    for chromosome in d_QUAL.keys():
        for line in d_QUAL[chromosome]:
            l = line.split()
            QUAL1 = float(l[0])
            QUAL2 = float(l[1])
            l_QUAL1 += [QUAL1]
            l_QUAL2 += [QUAL2]
    instance = statistics.tests()
    r = instance.correlation(l_QUAL1,l_QUAL2,)
    print r

    stop

    return


def parse_line_vcf(
    line, bool_ignore_FILTER=False,
    ):

    l = CHROM = POS = REF = ALT = FILTER = bool_continue = None
    if line[0] == '#':
        bool_continue = True
    else:
        l = line.strip().split()
        FILTER = l[6]
        REF = l[3]
        ALT = l[4]
        if bool_ignore_FILTER == False and FILTER not in ['.','PASS',]:
            bool_continue = True
        elif len(REF) != 1:
            bool_continue = True
        elif len(ALT) != 1:
            bool_continue = True
        else:
            CHROM = l[0]
            POS = int(l[1])
            bool_continue = False

    return l, CHROM, POS, REF, ALT, FILTER, bool_continue


def count_and_plot():

    chromosome = '1'
    import time
    t1 = time.time()
    for l_fp_in in [
        ['out_GATK/join/CombineVariants.vcf'],
        ['out_GATK/join/ApplyRecalibration.recalibrated.filtered.vcf'],
        ['SelectVariants_discordance1.vcf'],
        ['SelectVariants_discordance2.vcf'],
        ['%s.vqsr.filt.vcf' %(chromosome) for chromosome in range(1,23)+['X','Y',]],
        ['SelectVariants_concordance.vcf'],
        ]:

        import time
        t1 = time.time()

        ##
        ## prepare scatter lists
        ##
        l_gnuplot_MAF = []
        l_gnuplot_DP = []
        l_gnuplot_CR = []

        ##
        ## prepare contour dic
        ##
        d_contour = {}
        for AF in xrange(100+1):
            d_contour[AF*0.01] = {}
            for DP in xrange(150+1):
                d_contour[AF*0.01][DP*10.] = 0
        d_contour_CR = {}
        for AF in xrange(100+1):
            d_contour_CR[round(AF*0.01,2)] = {}
            for CR in xrange(100+1):
                d_contour_CR[round(AF*0.01,2)][CR] = 0

        for fp_in in l_fp_in:
            print fp_in
            fd = open(fp_in,'r')

            print fp_in
            for line in fd:
                if line[0] == '#':
                    continue
##                if line.count('./.')+line.count('0/0')+line.count('0/1')+line.count('1/1') != 100:
##                    print line
##                    stop
                CHROM, d_INFO, bool_continue = parse_line(line)
                if bool_continue == True:
                    continue

                CR = 100-line.count('./.')

                DP = int(d_INFO['DP'])
                try:
                    AF = float(d_INFO['AF'])
                except:
                    d_INFO['AF']
                    AF = 'N/A'

                if AF < 0.5:
                    MAF = AF
                else:
                    MAF = 1-AF

                ##
                ## append to list
                ##
                l_gnuplot_DP += [DP]
                l_gnuplot_CR += [CR]
                if fp_in != 'mp15_vqsr.vcf':
                    if AF == 'N/A':
                        stop
                    l_gnuplot_MAF += [MAF]
                    if DP < 1500:
                        d_contour[0.01*round(MAF/0.01,0)][10.*round(DP/10.,0)] += 1

                ##
                ## append to dic
                ##
                d_contour_CR[round(MAF,2)][CR] += 1

                if chromosome != CHROM:
                    t2 = time.time()
                    print fp_in, '%-2s' %(chromosome), '%2is' %(int(t2-t1))
                    chromosome = CHROM
                    t1 = t2

##            if CHROM == '2':
##                break
##            if POS[-1] == '0' and POS[-2] == '0' and POS[-3] == '0' and POS[-4] == '0':
##                print '%2s %9s %6s %4s' %(CHROM, POS, AF, DP,), fp_in
##                break
##            if POS[-1] == '0' and POS[-2] == '0' and POS[-3] == '0' and POS[-4] == '0':
##                print '%2s %9s %6s %4s' %(CHROM, POS, AF, DP,), fp_in

        title = fp_in.replace('_','').replace('out_GATK/','').replace('.vcf','')
        suffix = fp_in.replace('out_GATK','').replace('/','').replace('.vcf','')
        
        gnuplot.histogram2(
            'DP_%s' %(suffix),
            l_data=l_gnuplot_DP,
            x_step=10,x_min=0,x_max=1000,tic_step=100,
            xlabel='DP from VCF',
            ylabel='SNP count',
            title= title,
            )
        gnuplot.histogram2(
            'CR_%s' %(suffix),
            l_data=l_gnuplot_CR,
            x_step=1,x_min=0,x_max=100,tic_step=10,
            xlabel='SNP Call Rate',
            ylabel='SNP count',
            title= title,
            )
        if fp_in != 'mp15_vqsr.vcf':
            gnuplot.histogram2(
                'AF_%s' %(suffix),
                l_data=l_gnuplot_MAF,
                x_min=0,x_max=.5,tic_step=0.05,x_step=0.01,
                xlabel='AF from VCF',
                ylabel='SNP count',
                title = title,
                )

            lines = []
            for AF in xrange(50+1):
                for DP in xrange(150+1):
                    lines += ['%s %s %s\n' %(AF*0.01,DP*10.,d_contour[AF*0.01][DP*10.],)]
                lines += ['\n']
            gnuplot.contour_plot(
                'AFvDP_%s' %(suffix),
                lines,
                title = title,
                xlabel = 'AF from VCF',
                ylabel = 'DP from VCF',
                zlabel = 'count',
                )

            lines = []
            for AF in xrange(50+1):
                for CR in xrange(100+1):
                    lines += ['%s %s %s\n' %(AF*0.01,CR,d_contour_CR[round(AF*0.01,2)][CR],)]
                lines += ['\n']
            gnuplot.contour_plot(
                'AFvCR_%s' %(suffix),
                lines,
                title = title,
                xlabel = 'AF from VCF',
                ylabel = 'Call Rate from VCF',
                zlabel = 'count',
                )

##        t2 = time.time()
##        print t2-t1
##        stop

    return


def parse_line(line):

    bool_continue = False
    CHROM = None
    d_INFO = None
    while True:

        l = line.split('\t')
        CHROM = l[0]
        POS = l[1]
        ID = l[2]
        REF = l[3]
        ALT = l[4]
        QUAL = l[5]
        FILTER = l[6]
        INFO = l[7]
        FORMAT = l[8]
        samples = l[9:]
        import time
        t1 = time.time()
        if FILTER == '.':
            print line
            print fp_in
            stop
        if FILTER != 'PASS':
            bool_continue = True
            break
        ## avoid InDels and other stuff...
        if len(REF) != 1:
            bool_continue = True
            break
        if len(ALT) != 1:
            bool_continue = True
            break

    ##            ## this should be much faster than doing a regex
    ##            ## 456 seconds (fails when variable length)
    ##            index1 = INFO.index(';AF')+1+3
    ##            index2 = index1+INFO[index1:].index(';')
    ##            AF = INFO[index1:index2]

    ##            ## 477 seconds (fails when annotation not identical for each line)
    ##            if bool_INFO_indexed == False:
    ##                bool_INFO_indexed = True
    ##                d_INFO_indexes = {}
    ##                l = INFO.split(';')
    ##                for index in range(len(l)):
    ##                    s_key_and_val = l[index]
    ##                    s_key = s_key_and_val[:s_key_and_val.index('=')]
    ##                    d_INFO_indexes[s_key] = index

        ## 798 seconds (check if regex is faster than this...)
        l_INFO = INFO.split(';')
        d_INFO = {}
        for index in xrange(len(l_INFO)):
            s_key_and_val = l_INFO[index]
            if not '=' in s_key_and_val:
                continue
            index_equal = s_key_and_val.index('=')
            s_key = s_key_and_val[:index_equal]
            s_val = s_key_and_val[index_equal+1:]
            d_INFO[s_key] = s_val

        break

    return CHROM, d_INFO, bool_continue


if __name__ == '__main__':
    main()
