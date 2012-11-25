#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import os, math
import sys
sys.path.append('/nfs/users/nfs_t/tc9/github/ms23/GATK_pipeline')
import GATK_pipeline
sys.path.append('/nfs/users/nfs_t/tc9/github/tc9/math')
import statistics
instance_statistics = statistics.tests()

def main():

    instance = GATK_pipeline.main()
    d_chromosome_lengths = instance.parse_chromosome_ranges()
    d_centromere_ranges = parse_centromere_ranges()

    ped = 'omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped'
    sepjoin = 'sep'

    transpose(ped)

    IMPUTE2_tped(ped,d_chromosome_lengths,sepjoin,d_centromere_ranges,)

    BEAGLE_tped(ped,d_chromosome_lengths,sepjoin,d_centromere_ranges,)

    return


def parse_centromere_ranges():

    d_centromere_ranges = {}
    fd = open('centromeres_hg19.txt','r')
    lines = fd.readlines()
    fd.close()
    for line in lines:
        l = line.split()
        chrom = l[0][3:]
        init = int(l[1])
        term = int(l[2])
        if not chrom in d_centromere_ranges.keys():
            d_centromere_ranges[chrom] = []
        d_centromere_ranges[chrom] += [init,term,]

    return d_centromere_ranges


def parse_line_gen(line):

    l = line.strip().split()
    marker = l[0]
    pos = int(l[2])
    alleleA = l[3]
    alleleB = l[4]

    return l, marker, pos, alleleA, alleleB


def parse_line_gprobs(line_gprobs):

    l_gprobs = line_gprobs.strip().split()
    marker_gprobs = l_gprobs[0]
    pos_gprobs = int(marker_gprobs.split(':')[1])
    alleleA_gprobs = l_gprobs[1]
    alleleB_gprobs = l_gprobs[2]

    return l_gprobs, marker_gprobs, pos_gprobs, alleleA_gprobs, alleleB_gprobs


def BEAGLE_ped(ped):

    l_r2_allelic = []

    ##
    ## parse gprobs header and first line
    ##
    chromosome = '1'
    fp_gprobs = 'out_BEAGLE/sep/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.gprobs' %(
        chromosome,chromosome,
        )
    fd_gprobs = open(fp_gprobs,'r')
    ##
    ## parse header (parse samples in gprobs)
    ##
    for line in fd_gprobs:
        l = line.split()
        l_samples_gprobs = l_individuals = [l[i] for i in xrange(3,303,3)]
        break
    ##
    ## parse first line
    ##
    for line_gprobs in fd_gprobs:
        break
    (
        l_gprobs, marker_gprobs, pos_gprobs,
        alleleA_gprobs, alleleB_gprobs
        ) = parse_line_gprobs(line_gprobs)

    ##
    ## find relevant line numbers in ped
    ## instead of doing lookup in samples from gprobs all the time
    ##
    fd_ped = open('%s.ped' %(ped),'r')
    l_indexes_samples = []
    for line_ped in fd_ped:
        sample_ped = line_ped[0:10]
        l_indexes_samples += [l_samples_gprobs.index(sample_ped)]
    del l_samples_gprobs

    ##
    ## loop over markers from map file
    ##
    row_map = -1
    fd_map = open('%s.map' %(ped),'r')
    for line_map in fd_map:
        l_map = line_map.strip().split()
        chromosome_map = l_map[0]
        marker_map = l_map[1]
        pos_map = int(l_map[3])
        row_map += 1
        if chromosome != chromosome_map:
            fd_gprobs.close()
            chromosome = chromosome_map
            fp_gprobs = 'out_BEAGLE/sep/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.gprobs' %(
                chromosome,chromosome,
                )
            fd_gprobs = open(fp_gprobs,'r')
            for i in xrange(2):
                for line_gprobs in fd_gprobs:
                    break
            (
                l_gprobs, marker_gprobs, pos_gprobs,
                alleleA_gprobs, alleleB_gprobs,
                ) = parse_line_gprobs(line_gprobs)
        ## keep looping over markers from map file
        if pos_gprobs > pos_map:
            continue
        ## loop over markers from gprobs file
        elif pos_map > pos_gprobs:
            for line_gprobs in fd_gprobs:
                (
                    l_gprobs, marker_gprobs, pos_gprobs,
                    alleleA_gprobs, alleleB_gprobs,
                    ) = parse_line_gprobs(line_gprobs)
                ## keep looping over markers from map file
                if pos_gprobs > pos_map:
                    break
                ## keep looping over markers from gprobs file
                elif pos_map > pos_gprobs:
                    continue
                ## identical marker found in gprobs file
                else:
##                    import time
##                    t1 = time.time()
                    r2 = allelic_R2(ped,row_map,l_indexes_samples,l_gprobs,alleleA_gprobs,alleleB_gprobs,)
##                    t2 = time.time()
##                    print t2-t1
##                    print (t2-t1)*1622321/3600.
##                    print r2, marker_map, marker_gprobs
##                    stop
                    l_r2_allelic += ['%s %s\n' %(marker_map,r2,)]
        ## identical marker found in map file
        else:
            r2 = allelic_R2(ped,row_map,l_indexes_samples,l_gprobs,alleleA_gprobs,alleleB_gprobs,)
            print r2, marker_map, marker_gprobs
            l_r2_allelic += ['%s %s\n' %(marker_map,r2,)]
    fd_gprobs.close()

    fd = open('allelicR2_BEAGLEped.txt','w')
    fd.writelines(l_r2_allelic)
    fd.close()

    return


def allelic_R2_tped(
    l_tped,l_gprobs,
    n_samples,l_indexes_samples,
    alleleA_gprobs,alleleB_gprobs,
    cols_skip,
    ):

    sum_xy = 0
    sum_xx = 0
    sum_yy = 0
    sum_x = 0
    sum_y = 0
    n = 0
    count_0 = 0
    count_1 = 0
    count_2 = 0
    bool_break = False ## break if discordance (IMPUTE2, not BEAGLE)
    for i in xrange(n_samples):
        alleleA_ped = l_tped[4+2*i]
        alleleB_ped = l_tped[4+2*i+1]
        index = l_indexes_samples[i]
        if alleleA_ped == '0' and alleleB_ped == '0':
            continue
        if alleleA_ped == '0' or alleleB_ped == '0':
            print l_tped
            stop
        if alleleA_ped == alleleA_gprobs:
            if alleleB_ped == alleleA_gprobs:
                y = dosage_genotype = 0.
                count_0 += 1
            elif alleleB_ped == alleleB_gprobs:
                y = dosage_genotype = 1.
                count_1 += 1
            else:
                print alleleA_ped, alleleB_ped, alleleA_gprobs, alleleB_gprobs
                stop
        elif alleleA_ped == alleleB_gprobs:
            if alleleB_ped == alleleA_gprobs:
                y = dosage_genotype = 1.
                count_1 += 1
            elif alleleB_ped == alleleB_gprobs:
                y = dosage_genotype = 2.
                count_2 += 1
            else:
                print alleleA_ped, alleleB_ped, alleleA_gprobs, alleleB_gprobs
                stop
        else:
            print 'discordance',
            print alleleA_ped, alleleB_ped,
            print alleleA_gprobs, alleleB_gprobs, l_gprobs[0]
            bool_break = True
            break

        ##
        ## dosage
        ##
        l = l_gprobs[cols_skip+3*index:cols_skip+3*index+3]
##                        x = dosage_imputation = 0*l[0]+l[1]+2*l[2]
        x = dosage_imputation = float(l[1])+2*float(l[2])

        ##
        ## sums
        ##
        sum_xy += x*y
        sum_xx += x*x
        sum_yy += y*y
        sum_x += x
        sum_y += y
        n += 1

    if n > 1 and bool_break == False:
        denominator = math.sqrt((sum_xx-sum_x**2/n)*(sum_yy-sum_y**2/n))
        if denominator == 0:
            r = 1.
        else:
            nominator = (sum_xy-sum_x*sum_y/n)
            r = nominator/denominator
        r2 = r**2
##        ##
##        ## p
##        ##
##        DF = n-2
##        if r > 0.999999:
##            p_correlation = 0
##        else:
##            try:
##                t_correlation = r*math.sqrt((n-2)/(1-r**2))
##            except:
##                print r, r**2
##                stop
##            p_correlation = instance_statistics.tdist(t_correlation,DF)
##        ##
##        ##
##        ##
##        if r < 0.01:
####            print l_tped
####            print l_gprobs
####            print r
####            print nominator, denominator
##            print count_0, count_1, count_2, r
    else:
        r2 = -1

    return r2


def allelic_R2(ped,row_map,l_indexes_samples,l_gprobs,alleleA_gprobs,alleleB_gprobs,):

    fd_ped = open('%s.ped' %(ped),'r')
    sum_xy = 0
    sum_xx = 0
    sum_yy = 0
    sum_x = 0
    sum_y = 0
    n = 0
    bool_break = False
    i_line_ped = -1
    for line_ped in fd_ped:
        i_line_ped += 1
        sample_ped = line_ped[0:10]
        column_ped = 31+4*row_map
        alleleA_ped = line_ped[column_ped]
        alleleB_ped = line_ped[column_ped+2]
        if alleleA_ped == '0' and alleleB_ped == '0':
            continue
        if alleleA_ped == alleleA_gprobs:
            if alleleB_ped == alleleA_gprobs:
                y = dosage_genotype = 0.
            elif alleleB_ped == alleleB_gprobs:
                y = dosage_genotype = 1.
            else:
                print alleleA_ped, alleleB_ped, alleleA_gprobs, alleleB_gprobs
                stop
        elif alleleA_ped == alleleB_gprobs:
            if alleleB_ped == alleleA_gprobs:
                y = dosage_genotype = 1.
            elif alleleB_ped == alleleB_gprobs:
                y = dosage_genotype = 2.
            else:
                print alleleA_ped, alleleB_ped, alleleA_gprobs, alleleB_gprobs
                stop
        else:
            print alleleA_ped, alleleB_ped
            print alleleA_gprobs, alleleB_gprobs
            print sample_ped
            stop
            bool_break = True
            break

        ##
        ## dosage
        ##
        index = l_indexes_samples[i_line_ped]
        l = l_gprobs[3+3*index:3+3*index+3]
##                        x = dosage_imputation = 0*l[0]+l[1]+2*l[2]
        x = dosage_imputation = float(l[1])+2*float(l[2])

        sum_xy += x*y
        sum_xx += x*x
        sum_yy += y*y
        sum_x += x
        sum_y += y
        n += 1

    fd_ped.close()

    if bool_break == False and n > 1:
        r = (sum_xy-sum_x*sum_y/n)/math.sqrt((sum_xx-sum_x**2/n)*(sum_yy-sum_y**2/n))
        r2 = r**2
        print n,
    else:
        print alleleA_ped
        print alleleB_ped
        print n
        print marker_map
        print marker_gprobs
        print pos_map
        print pos_gprobs
        stop

    return r2


def transpose(ped):

    cmd = 'plink \\\n'
    cmd += '--noweb \\\n'
    cmd += '--file %s \\\n' %(ped)
    cmd += '--recode \\\n'
    cmd += '--transpose \\\n'
    cmd += '--out %s \\\n' %(ped)

    if not os.path.isfile('%s.tped' %(ped)):
        print cmd
        os.system(cmd)

    if not os.path.isfile('%s.tped' %(ped)):
        stop

    return


def IMPUTE2_tped(ped,d_chromosome_lengths,sepjoin,d_centromere_ranges,):

    i_first_line = 1

    ##
    ## parse header (parse samples in gprobs)
    ##
    chromosome = '1'
    fd_gprobs = open(
        'out_BEAGLE/sep/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.gprobs' %(
            chromosome,chromosome,
            ),'r')
    for line in fd_gprobs:
        l = line.split()
        l_samples_gprobs = l_individuals = [l[i] for i in xrange(3,303,3)]
        break
    fd_gprobs.close()

    ##
    ## find relevant columns in tped (parse samples in tfam)
    ##
    fd_tfam = open('%s.tfam' %(ped),'r')
    l_indexes_samples = []
    for line_tfam in fd_tfam:
        sample_tfam = line_tfam[0:10]
        l_indexes_samples += [l_samples_gprobs.index(sample_tfam)]
    fd_tfam.close()
    del l_samples_gprobs

    n_samples = len(l_indexes_samples)

    ##
    ## parse first data line of gen file
    ##
    fd_gen = open('out_IMPUTE2/sep/IMPUTE2.%s.gen' %(chromosome),'r')
    for line_gen in fd_gen:
        break
    (
        l_gen, marker_gen, pos_gen,
        alleleA_gen, alleleB_gen
        ) = parse_line_gen(line_gen)

    ##
    ## loop over markers in tped file
    ##
    l_r2_allelic = []
    pos_add = 0
    print 'chromosome', chromosome
    fd_tped = open('%s.tped' %(ped),'r')
    for line_tped in fd_tped:
        l_tped = line_tped.strip().split()
        chromosome_tped = l_tped[0]
        marker_tped = l_tped[1]
        pos_tped = int(l_tped[3])
        if chromosome != chromosome_tped:
            print 'chromosome', chromosome_tped
            fd_gen.close()
            pos_add += d_chromosome_lengths[chromosome]
            chromosome = chromosome_tped
            fd_gen = open('out_IMPUTE2/sep/IMPUTE2.%s.gen' %(chromosome),'r')
            ## parse first line of gen file
            for i in xrange(i_first_line):
                for line_gen in fd_gen: break
            (
                l_gen, marker_gen, pos_gen,
                alleleA_gen, alleleB_gen,
                ) = parse_line_gen(line_gen)
        if pos_tped < pos_gen:
            continue
        elif pos_tped > pos_gen:
            for line_gen in fd_gen:
                (
                    l_gen, marker_gen, pos_gen,
                    alleleA_gen, alleleB_gen
                    ) = parse_line_gen(line_gen)
                if pos_tped > pos_gen:
                    continue
                elif pos_gen > pos_tped:
                    break
                else:
                    r2 = allelic_R2_tped(
                        l_tped,l_gen,
                        n_samples,l_indexes_samples,
                        alleleA_gen,alleleB_gen,
                        5,
                        )
##                    print r2, marker_tped, marker_gen
##                    l_r2_allelic += ['%s %s\n' %(marker_gen,r2,)]
                    l_r2_allelic += [
                        '%s %s\n' %(
                            pos_gen+pos_add,r2,
                            )
##                        r2,
##                        '%s %s\n' %(
##                            r2,p,
##                            )
                        ]
        else:
            r2 = allelic_R2_tped(
                l_tped,l_gen,
                n_samples,l_indexes_samples,
                alleleA_gen,alleleB_gen,
                5,
                )
##            print r2, marker_tped, marker_gprobs
##            l_r2_allelic += ['%s %s\n' %(marker_gprobs,r2,)]
            l_r2_allelic += [
                '%s %s\n' %(
                    pos_gen+pos_add,r2,
                    )
##                r2,
##                '%s %s\n' %(
##                    r2,p,
##                    )
                ]
    fd_gen.close()
    fd_tped.close()

    fd = open('allelicR2_%s_IMPUTE2.gnuplotdata' %(sepjoin),'w')
    fd.writelines(l_r2_allelic)
    fd.close()

##    print sum(l_r2_allelic)/len(l_r2_allelic)

    gnuplot('allelicR2_%s_IMPUTE2' %(sepjoin),d_chromosome_lengths,d_centromere_ranges,)

    return


def IMPUTE2_low_mem(n_individuals):

    sum_xy = 0
    sum_xx = 0
    sum_yy = 0
    sum_x = 0
    sum_y = 0

    fd = fd_map_omni = open('omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.map','r')
    fd_ped_omni = open('omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.ped','r')
    print len(fd_map_omni.readlines())
    print len(fd_ped_omni.readlines())
    stop

##    i = 0
##    d = {}
##    for line_map_omni in fd_map_omni:
##        if line_map_omni[:2] != '18':
##            continue
##        pos_omni = line_map_omni.split()[3]
##        d[pos_omni] = [i,line_map_omni,]
##        i += 1

    fd = open('IMPUTE2.xx.gen','r')
    bool_continue = False
    for line in fd:
        l = line.split()
        pos_wgs = l[2]
        if bool_continue == True:
            print pos_wgs, pos_omni
            if int(pos_wgs) > int(pos_omni):
                bool_continue = False
            elif pos_wgs == pos_omni:
                print line
                print line_map_omni
                stop0
                bool_continue = False
            ## wgs still lagging
            else:
                continue

        for line_map_omni in fd_map_omni:
            pos_omni = line_map_omni.split()[3]
            if line_map_omni[:2] != 'xx':
                continue
            print pos_wgs, pos_omni
            if int(pos_wgs) < int(pos_omni):
                bool_continue = True
                break
            stop
            if pos_wgs == pos_omni:
                print line
                print line_map_omni
                stop1

        if bool_continue == True:
            continue
        if d.has_key(pos_wgs):
            print d[pos_wgs]
        else:
            continue
        stop2
        for i in xrange(5,3*n_individuals+5,3,):
##            dosage_imputation = 0*l[i]+1*l[i+1]+2*l[i+2]
            x = dosage_imputation = l[i+1]+2*l[i+2]
            y = dosage_genotype
            sum_xy += x*y
            sum_xx += x*x
            sum_yy += y*y
            sum_x += x
            sum_y += y
        stop
    fd.close()

    r = (sum_xy-sum_x*sum_y/n)/math.sqrt((sum_xx-sum_x**2/n)*(sum_yy-sum_y**2/n))
    r2 = r**2

    return r2


def omni2dic():

    print 'read map'

    fd = fd_map_omni = open('omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.map','r')

    i = 0
    d = {}
    l_markers = []
    for line_map_omni in fd_map_omni:
        l = line_map_omni.strip().split()
        chrom_omni = l[0]
        pos_omni = l[3]
        marker_omni = l[1]
        l_markers += [marker_omni]
##        d[pos_omni] = [i,line_map_omni,]
        d[marker_omni] = i
        i += 1

    return d, l_markers


def BEAGLE_count():

    ##
    d_map, l_map = omni2dic()

    ##
    ## loop chromosomes
    ##
    print 'read gprobs'
    l_wgs_markers_all = []
    wgs_not_omni = 0
    wgs = 0
    wgs_omni_intersect = 0
    for chromosome in xrange(1,22+1):
        print chromosome
        l_wgs_markers = []
        ##
        ## read gprobs
        ##
        fd = open(
            'out_BEAGLE/sep/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.gprobs' %(
                chromosome,chromosome,
                ),'r')
        ##
        ## header
        ##
        for line in fd:
            l = line.split()
            l_samples = l_individuals = l[3:]
            break
        for line in fd:
            l = line.split()
            marker = l[0]
            l_wgs_markers += [marker]
            l_wgs_markers_all += [marker]
        fd.close()
        set_omni_markers = set(l_map)
        set_wgs_markers = set(l_wgs_markers)
        wgs += len(set_wgs_markers)
        wgs_not_omni += len(set_wgs_markers-set_omni_markers)
        wgs_omni_intersect += len(set_wgs_markers&set_omni_markers)
    print 'wgs', wgs
    print 'wgs not omni', wgs_not_omni
    print 'wgs and omni overlap', wgs_omni_intersect

##    ##
##    ## initiate dictionary
##    ##
##    d_ped = {}
##    for sample in l_samples:
##        d_ped[sample] = {}
##    ## read ped
##    fd_ped_omni = open('omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.ped','r')
##    ## count markers
##    n_map = len(l_map)
##    ## loop ped lines
##    for line_ped in fd_ped_omni:
##        sample = line_ped[0:10]
##        for row_map in xrange(n_map):
##            marker = l_map[row_map]
##            column_ped = 31+4*row_map
##            alleleA = line_ped[column_ped]
##            alleleB = line_ped[column_ped+2]
##            d_ped[sample][marker] = [alleleA,alleleB,]

    set_omni_markers = set(l_map)
    set_wgs_markers = set(l_wgs_markers_all)
    print 'omni', len(set_omni_markers)
    print 'omni not wgs', len(set_omni_markers-set_wgs_markers)
    stop

    n_samples = len(l_samples)
    sum_xy = 0
    sum_xx = 0
    sum_yy = 0
    sum_x = 0
    sum_y = 0
    ## loop over markers
    for line in fd:
        l = line.split()
        marker = l[0]
        alleleA_imputation = l[1]
        alleleB_imputation = l[2]
        if d_map.has_key(marker):
            print d_map[marker]
            print l[:6]
            for i_sample in xrange(n_samples):
##                x = dosage_imputation = 0*l[6+i_sample+0]+l[6+i_sample+1]+2*l[6+i_sample+1]
                x = dosage_imputation = l[6+i_sample+1]+2*l[6+i_sample+1]
                sample = l_samples[i_sample]
                alleles = d_ped[sample][marker]
                print alleleA_imputation
                print alleleB_imputation
                alleleA_omni = alleles[0]
                alleleB_omni = alleles[1]
                print alleleA_omni
                print alleleB_omni
                ## homozygous
                if alleleA_omni == alleleA_imputation and alleleB_omni == alleleA_imputation:
                    y = dosage_genotype = 0
                ## homozygous
                elif alleleA_omni == alleleB_imputation and alleleB_omni == alleleB_imputation:
                    y = dosage_genotype = 2
                ## homozygous
                else:
                    print alleleA_omni, alleleA_imputation
                    print alleleB_omni, alleleB_imputation
                    stop
                stop
            stop

            x = dosage_imputation = l[i+1]+2*l[i+2]
            y = dosage_genotype
            sum_xy += x*y
            sum_xx += x*x
            sum_yy += y*y
            sum_x += x
            sum_y += y
    fd.close()
    stop

    r = (sum_xy-sum_x*sum_y/n)/math.sqrt((sum_xx-sum_x**2/n)*(sum_yy-sum_y**2/n))
    r2 = r**2
        
    stop
    return r2


def BEAGLE_tped(ped,d_chromosome_lengths,sepjoin,d_centromere_ranges,):

    ##
    ## parse header (parse samples in gprobs)
    ##
    chromosome = '1'
    fd_gprobs = open(
        'out_BEAGLE/sep/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.gprobs' %(
            chromosome,chromosome,
            ),'r')
    for line in fd_gprobs:
        l = line.split()
        l_samples_gprobs = l_individuals = [l[i] for i in xrange(3,303,3)]
        break

    ##
    ## find relevant columns in tped (parse samples in tfam)
    ##
    fd_tfam = open('%s.tfam' %(ped),'r')
    l_indexes_samples = []
    for line_tfam in fd_tfam:
        sample_tfam = line_tfam[0:10]
        l_indexes_samples += [l_samples_gprobs.index(sample_tfam)]
    fd_tfam.close()
    del l_samples_gprobs

    n_samples = len(l_indexes_samples)

    ##
    ## parse first data line of gprobs file
    ##
    for line_gprobs in fd_gprobs:
        break
    (
        l_gprobs, marker_gprobs, pos_gprobs,
        alleleA_gprobs, alleleB_gprobs
        ) = parse_line_gprobs(line_gprobs)

    ##
    ## loop over markers in tped file
    ##
    l_r2_allelic = []
    pos_add = 0
    print 'chromosome', chromosome
    fd_tped = open('%s.tped' %(ped),'r')
    for line_tped in fd_tped:
        l_tped = line_tped.strip().split()
        chromosome_tped = l_tped[0]
        marker_tped = l_tped[1]
        pos_tped = int(l_tped[3])
        if chromosome != chromosome_tped:
            print 'chromosome', chromosome_tped
            fd_gprobs.close()
            pos_add += d_chromosome_lengths[chromosome]
            chromosome = chromosome_tped
            fd_gprobs = open(
                'out_BEAGLE/sep/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.gprobs' %(
                    chromosome,chromosome,
                    ),'r')
            ## parse first line of gprobs file
            for i in xrange(2):
                for line_gprobs in fd_gprobs: break
            (
                l_gprobs, marker_gprobs, pos_gprobs,
                alleleA_gprobs, alleleB_gprobs
                ) = parse_line_gprobs(line_gprobs)
        if pos_tped < pos_gprobs:
            continue
        elif pos_tped > pos_gprobs:
            for line_gprobs in fd_gprobs:
                (
                    l_gprobs, marker_gprobs, pos_gprobs,
                    alleleA_gprobs, alleleB_gprobs
                    ) = parse_line_gprobs(line_gprobs)
                if pos_tped > pos_gprobs:
                    continue
                elif pos_gprobs > pos_tped:
                    break
                else:
                    r2 = allelic_R2_tped(
                        l_tped,l_gprobs,
                        n_samples,l_indexes_samples,
                        alleleA_gprobs,alleleB_gprobs,
                        3,
                        )
##                    print r2, marker_tped, marker_gprobs
##                    l_r2_allelic += ['%s %s\n' %(marker_gprobs,r2,)]
                    l_r2_allelic += [
                        '%s %s\n' %(
                            pos_gprobs+pos_add,r2,
                            )
##                        r2,
##                        '%s %s\n' %(
##                            r2,p,
##                            )
                        ]
        else:
            r2 = allelic_R2_tped(
                l_tped,l_gprobs,
                n_samples,l_indexes_samples,
                alleleA_gprobs,alleleB_gprobs,
                3,
                )
##            print r2, marker_tped, marker_gprobs
##            l_r2_allelic += ['%s %s\n' %(marker_gprobs,r2,)]
            l_r2_allelic += [
                '%s %s\n' %(
                    pos_gprobs+pos_add,r2,
                    )
##                r2,
##                '%s %s\n' %(
##                    r2,p,
##                    )
                ]
    fd_gprobs.close()
    fd_tped.close()

    fd = open('allelicR2_%s_BEAGLE.gnuplotdata' %(sepjoin),'w')
    fd.writelines(l_r2_allelic)
    fd.close()

##    print sum(l_r2_allelic)/len(l_r2_allelic)

    gnuplot('allelicR2_%s_BEAGLE' %(sepjoin),d_chromosome_lengths,d_centromere_ranges,)

    return


def gnuplot(prefix,d_chromosome_lengths,d_centromere_ranges,):

    lines = []
    lines += ['set terminal postscript enhanced eps\n']
    lines += ['set output "%s.ps"\n' %(prefix)]
    lines += ['set size 8,8\n']
    lines += ['set autoscale fix\n']
    lines += ['set encoding iso_8859_1\n']
    lines += ['set xlabel "{/=48 position}"\n']
    lines += ['set ylabel "{/=48 allelic R^2}"\n']
    lines += ['set title "{/=48 %s}"\n' %(prefix.replace('_',' '))]
    pos = 0
    for chromosome in range(1,22+1,):#+['X','Y',]:
        chromosome = str(chromosome)
        ## centromeres
        lines += ['set arrow from %s.0,0 to %s.0,1.1 nohead lc 1\n' %(
            pos+min(d_centromere_ranges[chromosome]),
            pos+min(d_centromere_ranges[chromosome]),
            )]
        lines += ['set arrow from %s.0,0 to %s.0,1.1 nohead lc 1\n' %(
            pos+max(d_centromere_ranges[chromosome]),
            pos+max(d_centromere_ranges[chromosome]),
            )]
        ## telomeres
        pos += d_chromosome_lengths[chromosome]
        lines += ['set arrow from %s.0,0 to %s.0,1.1 nohead lc 0\n' %(pos,pos,)]
    s_plot = 'plot [0:%i.0][0:1] "%s.gnuplotdata" ps 1 t ""\n' %(pos,prefix,)
    lines += [s_plot]

    fd = open('%s.gnuplotsettings' %(prefix),'w')
    fd.writelines(lines)
    fd.close()

    os.system('gnuplot %s.gnuplotsettings' %(prefix))

##    cmd = 'convert -rotate 90 %s.ps %s.png' %(prefix,prefix,)
    cmd = 'convert %s.ps %s.png' %(prefix,prefix,)
    print cmd
    os.system(cmd)
    if os.path.isfile('%s.png' %(prefix)):
        os.remove('%s.ps' %(prefix))
        os.remove('%s.gnuplotsettings' %(prefix))

    return


if __name__ == '__main__':
    main()
