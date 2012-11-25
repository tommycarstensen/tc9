#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import sys, os, time
import PEDgenotype_vs_VCFMAF

##
## this script generates concordance tables
##

##
## todo: 2012-08-16: also do this as a function of MAF
##

def main():

##    prefix = sys.argv[-1]
##    if not prefix in ['ApplyRecalibration','IMPUTE2','BEAGLE',]:
##        print 'Invalid prefix', sys.argv[-1]
##        raise Exception

    sepjoin = 'join'

    print 'parse SNP differences'
    print

    s_table = ''
    for prefix in ['ApplyRecalibration','BEAGLE','IMPUTE2',]:

        ## reset counts
        i_minor = 0
        i_major = 0
        i_both = 0
        ## loop over chromosomes
        for chromosome in range(1,22+1,):#+['X','Y',]:

            ##
            ## parse allele frequencies
            ##
##            print 'parse allele frequencies', chromosome
##            fp_in_vcf = 'out_GATK/sep/CombineVariants.%s.vcf' %(chromosome)
##            d_AF = PEDgenotype_vs_VCFMAF.loop_vcf_lines(fp_in_vcf)

            fp_in = 'out_plink/%s/plink_merge_%s_%s.diff' %(
                sepjoin,prefix,chromosome,
                )
            bool_continue = check_file_in(fp_in)
            if bool_continue == True: continue

            print 'loop diff lines', chromosome
            i_minor,i_major,i_both = loop_diff_lines(
                fp_in,i_minor,i_major,i_both,
##                d_AF,
                )

        if i_both == 0:
            continue
        s_table += '%s\n' %(prefix)
        s_table += 'major %s%% %s\n' %(
            round(100*float(i_major)/sum([i_minor,i_major,i_both]),1), i_major,
            )
        s_table += 'minor %s%% %s\n' %(
            round(100*float(i_minor)/sum([i_minor,i_major,i_both]),1), i_minor,
            )
        s_table += 'both %s%% %s\n' %(
            round(100*float(i_both)/sum([i_minor,i_major,i_both]),1), i_both,
            )
        s_table += '--------\n'
    print s_table
    fd = open('concordance_%s.table' %(sepjoin),'w')
    fd.write(s_table)
    fd.close()

    return


def loop_diff_lines(
    fp_in,i_minor,i_major,i_both,
##    d_AF,
    ):

    fd = open(fp_in)

    ## skip header
    for line in fd:
        break
    ## loop over lines
    for line in fd:
        l = line.strip().split()
        POS = l[0].split(':')[1]
        try:
            NEW = l[3].split('/')
        except:
            print line
            print fp_in
            print 'time', time.time()-os.path.getmtime(fp_in)
            stop
        OLD = l[4].split('/')
        ## e.g. G/G A/G
        if NEW[1] == OLD[1] and NEW[0] != OLD[0]:
            i_major += 1
        ## e.g. T/T T/C
        elif NEW[0] == OLD[0] and NEW[1] != OLD[1]:
            i_minor += 1
        ## e.g. T/T A/A
        elif NEW[0] != OLD[0] and NEW[1] != OLD[1]:
            i_both += 1
        else:
            print line
            stop

    fd.close()

    return i_minor,i_major,i_both


def check_file_in(fp_in):

    if not os.path.isfile(fp_in):
        print 'missing', fp_in
        bool_continue = True
    elif time.time()-os.path.getmtime(fp_in) < 60:
        print 'modified within the last 60 seconds', fp_in
        bool_continue = True
    else:
        bool_continue = False

    return bool_continue


if __name__ == '__main__':
    main()
