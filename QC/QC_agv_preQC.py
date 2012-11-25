#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import sys, os, time
sys.path.append('/nfs/users/nfs_t/tc9/lib/python2.7/site-packages')
import xlrd

l_populations = [
    u'Muganda', u'KIKUYU', u'Mandinka', u'ZULU', u'KALENJIN', u'Fula', u'Murundi', u'Munyarwanda', u'Sotho', u'Wolloff', u'Jola', u'GA-ADANGBE',
    u'IBO',
    'Ethiopia', ## SOMALI,AMHARA,OROMO
    ## 1000g
    'YRI','TSI','PUR','PEL','MXL','MKK','LWK','KHV','JPT','IBS','IBO','GIH','GBR','FIN','CLM','CHS','CHD','CHB','CEU','CDX','ASW','ACB',
    ]

def main():

    ## 0) check that he have the liftOver executable and the map.chain file
    init()

    bfile = 'omni2.5-4_20120904_agv_gtu' ## requires 2797MB

    ##
    ## update from build 36 to build 37
    ## and update sexes at the same time
    ##

    ## 1a) create liftOver input BED file from PLINK bim file
    create_input_liftOver(bfile,)

    ## 1b) run liftOver
    run_liftOver(bfile,)

    ## 1c) create PLINK input file from liftOver output BED file
    create_input_PLINK_update(bfile,)

    ## 1d) update quad SNP positions from build 36 to build 37
    ## i.e. run PLINK --update-map --update-sex
    PLINK_update_map_and_sex(bfile,)

    ## 2) exclude duplicates *after* updating from build 36 to build 37
    ## otherwise SNP4-69432884, SNP4-69611990 and SNP6-9123955 will not be excluded
    ## omni2.5-8_agv_20120910_gtu.bim:4        rs4860941       0       68963570        A       G
    ## omni2.5-8_agv_20120910_gtu.bim:4        SNP4-69432884   0       69432884        A       G
    ## omni2.5-4_20120904_agv_gtu.liftOver.BED:chr4    69280974        69280975        rs4860941
    ## omni2.5-4_20120904_agv_gtu.liftOver.BED:chr4    69280974        69280975        SNP4-69432884
    PLINK_exclude_duplicates()

    ## 3a) identify common positions and build a rsID dictionary (quad2octo) for those SNPs
    identify_common_positions()

    ## 3b) identify but do not exclude allele mismatches yet
    ## e.g. fgrep -w 106129520 *ted.bim *gtu.bim mismatch.SNPs *duplicates*
    ## omni2.5-4_20120904_agv_gtu_build37_sexupdated.bim:1     SNP1-105931043  0       106129520       A       G
    ## omni2.5-8_agv_20120910_gtu.bim:1        kgp15817565     0       106129520       A       C
    ## e.g. fgrep -w 207131740 *ted.bim *gtu.bim mismatch.SNPs *duplicates* *.BED
    ## omni2.5-4_20120904_agv_gtu_build37_sexupdated_deduplicated.bim:1        rs17018316      0       207131740       A       C
    ## omni2.5-8_agv_20120910_gtu_deduplicated.bim:1   kgp15398577     0       207131740       0       G
    ## omni2.5-4_20120904_agv_gtu.bim:1        SNP1-205198363  0       205198363       A       G
    ## omni2.5-4_20120904_agv_gtu.bim:1        rs17018316      0       205198363       A       C
    identify_allele_mismatches()

    ## 3c) rename quad SNPs
    ## i.e. run PLINK --update-map --update-name
    PLINK_rename()

    ## 4a) DO NOT flip SNPs
    ## i.e. DO NOT run PLINK --flip
    ## instead do this as part of QC.py (preQC of QC.py should perhaps be moved to a generic preQC.py script...)
    ## 4b1) DO NOT extract/exclude SNPs in/not in the strand file
    ## instead do this as part of QC.py
    ## 4b2) exclude/extract allele mismatches/matches
    ## 4b3) extract/exclude positions that are/are not common between quad and octo
    ## 4c) split by population (previously QC_agv_preQC_split_populations.py)
    ## i.e. run PLINK --keep
    split_by_population()

    ## 5) merge selected quad and octo samples (previously QC_agv_preQC_merge_chips.py)
    ## i.e. run PLINK --bmerge
    PLINK_bmerge()

    ## 6) count SNPs at each step for Manj
    count_SNPs()

    return


def count_SNPs():

    fd = open('SNPstats.txt','w')
    fd.close()

    cmd = 'cat omni2.5-4_20120904_agv_gtu.bim | wc -l'
    append_stats(cmd)
    cmd = "cat omni2.5-4_20120904_agv_gtu.bim | awk '{if($1>=1 && $1<=23) print}' | wc -l"
    append_stats(cmd)
    cmd = '''cat omni2.5-4_20120904_agv_gtu.unlifted.BED | awk 'NR%2==0{if ($1!="chr24" && $1!="chr25" && $1!="chr26") print}' | wc -l'''
    append_stats(cmd)
    cmd = "cat omni2.5-4_20120904_agv_gtu_build37_sexupdated.bim | wc -l"
    append_stats(cmd)

    fd = open('SNPstats.txt','a')
    fd.write('\n')
    fd.close()

    for bfile in [
        'omni2.5-8_agv_20120910_gtu',
        'omni2.5-4_20120904_agv_gtu_build37_sexupdated',
        ]:
        for suffix in ['.duplicates.SNPs','.bim','_deduplicated.bim',]:
            cmd = 'cat %s%s | wc -l' %(bfile,suffix,)
            append_stats(cmd)

    fd = open('SNPstats.txt','a')
    fd.write('\n')
    fd.close()

    for bfile in [
        'omni2.5-8_agv_20120910_gtu',
        'omni2.5-4_20120904_agv_gtu_build37_sexupdated',
        ]:
        for suffix in ['_deduplicated.mismatch.SNPs','_deduplicated.common.SNPs',]:
            cmd = 'cat %s%s | wc -l' %(bfile,suffix,)
            append_stats(cmd)

    fd = open('SNPstats.txt','a')
    fd.write('\n')
    fd.close()

    cmd = "cat quad2octo.dic | wc -l"
    append_stats(cmd)
    cmd = "cat omni2.5-4_20120904_agv_gtu_build37_sexupdated_deduplicated_renamed.bim | wc -l"
    append_stats(cmd)

    return


def split_by_population():

    l_bfiles = [
        'omni2.5-8_agv_20120910_gtu_deduplicated',
        'omni2.5-4_20120904_agv_gtu_build37_sexupdated_deduplicated_renamed',
        ]

    l_out = []

    for bfile in l_bfiles:
        for population in l_populations:

            ## 1000g
            if (
                len(population) == 3
                and
                population == population.upper()
                and
                type(population) == str
                ):
                population_out = population
                pass
            ## agv
            else:
                population_out = '-'.join(
                    [str(s[0]).upper()+str(s[1:]).lower() for s in population.split('-')]
                    )
                pass

            out = '%s_%s' %(bfile,population_out,)
            if os.path.isfile('%s.bed' %(out)):
                continue

            if population == 'Ethiopia':
                l_populations_sub = [u'SOMALI', u'AMHARA', u'OROMO',]
            else:
                l_populations_sub = [population,]

            ## set file name of .keep file
            keep = '%s_%s.keep' %(bfile,population_out,)

            ## init file
            fd = open(keep,'w')
            fd.close()

            for population_sub in l_populations_sub:
                ## append to file
                if not os.path.isfile('samples/%s.samples' %(population_sub)):
                    print population_sub
                    stop
                cmd = 'fgrep -f samples/%s.samples %s.fam' %(population_sub,bfile,)
                cmd += ' | '
                cmd += "awk '{print $1,$2}' >> %s" %(keep,)
##                print cmd
                os.system(cmd)
                continue
            cmd = 'cat %s | wc -l' %(keep)
            count = int(os.popen(cmd).read())
            ## no samples in bfile
            if count == 0:
                os.remove(keep)
                continue

            cmd = 'plink \\\n'
            cmd += '--bfile %s \\\n' %(bfile)
            cmd += '--keep %s \\\n' %(keep)
            cmd += '--exclude %s.mismatch.SNPs \\\n' %(l_bfiles[0])
            cmd += '--extract %s.common.SNPs \\\n' %(l_bfiles[0])
            cmd += '--noweb \\\n'
            cmd += '--nonfounders \\\n'
            cmd += '--make-bed \\\n'
            cmd += '--allow-no-sex \\\n'
            cmd += '--out %s' %(out)

##            print count
##            print int(os.popen('cat %s.fam | wc -l' %(bfile)).read())
            bsub = "bsub -M3000000 -R'select[mem>3000] rusage[mem=3000]' "
##            bsub += '-P agv '
            bsub += '-G agv '
            bsub += '-q normal '
            bsub += "-J'%s' " %(population_out)
            bsub += '%s' %(cmd)

            os.system(bsub)

            l_out += [out]

            ## continue loop over populations
            continue

        ## continue loop over bfiles
        continue

    ##
    ## check that output has been generated
    ##
    for x in xrange(8*60/5):
        print x
        bool_break = True
        for out in l_out:
            for extension in ['bed','bim','fam',]:
                fn_out = '%s_%s.%s' %(bfile,population_out,extension,)
                if not os.path.isfile(fn_out):
                    bool_break = False
                    break
                if os.path.getsize(fn_out) == 0:
                    bool_break = False
                    break
                continue
            if bool_break == False:
                break
            continue
        if bool_break == True:
            break
        time.sleep(5*60)
        continue

    return


def append_stats(cmd):

    print cmd
    i = int(os.popen(cmd).read())
    fd = open('SNPstats.txt','a')
    fd.write('%s\n' %(cmd))
    fd.write('%i\n' %(i))
    fd.close()

    return


def extract_SNP_selection():

    ## include
    incl1 = 'omni2.5-8_agv_20120910_gtu_deduplicated.common.SNPs'
    excl1 = 'omni2.5-8_agv_20120910_gtu_deduplicated.mismatch.SNPs'
    incl2 = 'bim_strand_intersection'
    cmd += '--exclude %s.mismatch.SNPs \\\n' %(bfile)
    cmd += '--extract %s.common.SNPs \\\n' %(bfile)
    print 'common - mismatch > extract.SNPs'
    stop

    l_bfiles = [
        'omni2.5-8_agv_20120910_gtu_deduplicated',
        'omni2.5-4_20120904_agv_gtu_build37_sexupdated_deduplicated_renamed',
        ]

    ## 1) select SNPs to be flipped
    fn_strand = 'HumanOmni2.5-8v1_A-b37.strand'
    cmd = '''cat %s | awk '{if($5=="-") print $1}' | sort > flip.txt''' %(fn_strand)
    print cmd
    os.system(cmd)

    ## 2) only SNPs to be extracted
    ## 2a) sort .strand file
    cmd = '''cat %s | sort -k1,1 > %s.sorted''' %(fn_strand,fn_strand,)
    print cmd
    os.system(cmd)
    ## 2a) sort .bim file
    cmd = '''cat %s | sort -k2,2 > %s.sorted''' %(l_bfiles[0],l_bfiles[0],)
    print cmd
    os.system(cmd)
    ## 2c) join sorted build 37 .strand and .bim file
    cmd = 'join -1 1 -2 2 -o 0 %s.sorted %s.bim.sorted > common_bim_strand.SNPs' %(
        fn_strand,l_bfiles[1],)
    print cmd
    os.system(cmd)
    stop

    for bfile in l_bfiles:
        cmd = 'plink \\\n'
        cmd += '--bfile %s \\\n' %(bfile,)
        if bfile == l_bfiles[1]:
            cmd += '--make-bed --out %s_flipped \\\n' %(bfile,)
        elif bfile == l_bfiles[0]:
            cmd += '--make-bed --out %s_flipped \\\n' %(bfile,)
        cmd += '--flip flip.txt \\\n'
        cmd += '--exclude %s.mismatch.SNPs \\\n' %(bfile)
        cmd += '--extract extract.SNPs \\\n' %(bfile)
        cmd += '--noweb --nonfounders --allow-no-sex \\\n'

    ## xx) clean up
    for bfile in l_bfiles:
        os.remove('%s.bim.sorted' %(bfile))
        os.remove('%s.common.SNPs' %(bfile))
        os.remove('%s.mismatch.SNPs' %(bfile))
    os.remove('%s.sorted' %(fn_strand))
    os.remove('extract.SNPs')

    return


def PLINK_rename():

    l_bfiles = [
        'omni2.5-8_agv_20120910_gtu_deduplicated',
        'omni2.5-4_20120904_agv_gtu_build37_sexupdated_deduplicated',
        ]

    bfile = 'omni2.5-4_20120904_agv_gtu_build37_sexupdated_deduplicated'

    cmd = 'plink \\\n'
    cmd += '--bfile %s \\\n' %(bfile,)
    cmd += '--make-bed --out %s_renamed \\\n' %(bfile,)
##    cmd += '--exclude %s.mismatch.SNPs \\\n' %(bfile)
##    cmd += '--extract %s.common.SNPs \\\n' %(bfile)
    cmd += '--update-map quad2octo.dic --update-name \\\n'
    cmd += '--noweb --nonfounders --allow-no-sex \\\n'
    print cmd
    os.system(cmd)

    return


def PLINK_bmerge():

    l_bfiles = [
        'omni2.5-8_agv_20120910_gtu_deduplicated',
        'omni2.5-4_20120904_agv_gtu_build37_sexupdated_deduplicated_renamed',
        ]

    ##
    ## define dic
    ##
    d_bfiles = {
        'omni2.5-4_20120904_agv_gtu_build37_sexupdated_deduplicated_renamed':'quad',
        'omni2.5-8_agv_20120910_gtu_deduplicated':'octo',
        }

    ##
    ## make clean up dirs
    ##
    for extension in ['nof','hh','log',]:
        if not os.path.isdir(extension):
            os.mkdir(extension)
            pass

    ##
    ## loop over populations
    ##
    for population in l_populations:

        ## 1000g
        if (
            len(population) == 3
            and
            population == population.upper()
            and
            type(population) == str
            ):
            population_out = population
            pass
        ## agv
        else:
            population_out = '-'.join(
                [str(s[0]).upper()+str(s[1:]).lower() for s in population.split('-')]
                )
            pass

        bool_both_exist = True
        for bfile in l_bfiles:
            out = '%s_%s' %(bfile,population_out,)
            if not os.path.isfile('%s.bed' %(out)):
                bool_both_exist = False
                break
            continue
        ## 1) all data present on one chip
        if bool_both_exist == False:
            for bfile in l_bfiles:
                out = '%s_%s' %(bfile,population_out,)
                if not os.path.isfile('%s.bed' %(out)):
                    continue
                for suffix in ['bed','fam','bim',]:
                    os.system('mv %s.%s %s.%s' %(
                        out,suffix,population_out,suffix,
                        ))
                    continue
                break
            pass
        ## 2) keep separate
        elif population_out in ['Muganda', 'Munyarwanda',]:
            for bfile in l_bfiles:
                fn1 = '%s_%s' %(bfile,population_out,)
                fn2 = '%s_%s' %(population_out,d_bfiles[bfile])
                for suffix in ['bed','fam','bim',]:
                    os.system('mv %s.%s %s.%s' %(
                        fn1,suffix,fn2,suffix,
                        ))
                    continue
                continue
            pass
        ## 3) merge
        else:

            l_out = ['%s_%s' %(bfile,population_out,) for bfile in l_bfiles]

            ##
            ## merge quad and octo (--bmerge)
            ##
            l_cmd = ['plink']
            l_cmd += ['--bmerge']
            l_cmd += ['%s.bed' %(l_out[0])]
            l_cmd += ['%s.bim' %(l_out[0])]
            l_cmd += ['%s.fam' %(l_out[0])]
            l_cmd += ['--bfile %s' %(l_out[1])]
            l_cmd += ['--make-bed --out %s' %(population_out)]
            l_cmd += ['--allow-no-sex']
            l_cmd += ['--nonfounders']
            l_cmd += ['--noweb']
            cmd = ' \\\n'.join(l_cmd)

            bsub = "bsub -M3000000 -R'select[mem>3000] rusage[mem=3000]' "
            bsub += '-P agv '
            bsub += '-q normal '
            bsub += "-J'%s' " %(population_out)

            os.system('%s %s' %(bsub, cmd))

            pass

        ## continue loop over populations
        continue

    return


def identify_common_positions():

    l_bfiles = [
        'omni2.5-4_20120904_agv_gtu_build37_sexupdated_deduplicated',
        'omni2.5-8_agv_20120910_gtu_deduplicated',
        ]

    col_bim_chr = 1
    col_bim_name = 2
    col_bim_pos = 4

    ##
    ## 1) sort
    ##
    sort(
        l_bfiles,
        'sort -k%i,%i -k%i,%i' %(
            col_bim_chr,col_bim_chr, col_bim_pos,col_bim_pos,
            ),
        )

    ##
    ## 2) identify common positions
    ##
    ## 2a) concatenate chromosome and position
    for bfile in l_bfiles:
        cmd = 'cat %s.bim.sorted' %(bfile)
        cmd += ' | sort -k%i,%i -k%i,%i' %(
            col_bim_chr,col_bim_chr, col_bim_pos,col_bim_pos,
            )
        cmd += ''' | awk '{print $0,$%i":"$%i}' ''' %(col_bim_chr, col_bim_pos,)
        cmd += ' > %s.bim.sorted.concatenated''' %(bfile,)
        print cmd
        os.system(cmd)
    ## 2b) join on concatenated field/column
    cmd = 'join -1 %i -2 %i' %(7,7,)
    cmd += ' -o 1.%i,2.%i,1.%i,1.%i,2.%i,2.%i' %(
        col_bim_name,col_bim_name,5,6,5,6,
        )
    cmd += ' %s.bim.sorted.concatenated' %(l_bfiles[0])
    cmd += ' %s.bim.sorted.concatenated' %(l_bfiles[1])
    cmd += ' > %s.%s.concatenated' %(l_bfiles[0],l_bfiles[1],)
    print cmd
    os.system(cmd)

    ## 3) write quad2octo.dic
    cmd = 'cat %s.%s.concatenated' %(l_bfiles[0],l_bfiles[1],)
    cmd += " | awk '{print $1,$2}'"
    cmd += ' > quad2octo.dic'
    print cmd
    os.system(cmd)

    ## 4) write %s.common.SNPs
    for i_bfile in xrange(len(l_bfiles)):
        bfile = l_bfiles[i_bfile]
        cmd = 'cat %s.%s.concatenated' %(l_bfiles[0],l_bfiles[1],)
        cmd += " | awk '{print $%i}'" %(i_bfile+1)
        cmd += ' > %s.common.SNPs' %(bfile)
        print cmd
        os.system(cmd)

    ##
    ## 4) clean up
    ##
    for bfile in l_bfiles:
        os.remove('%s.bim.sorted.concatenated' %(bfile))

    return


def identify_allele_mismatches():

    l_bfiles = [
        'omni2.5-4_20120904_agv_gtu_build37_sexupdated_deduplicated',
        'omni2.5-8_agv_20120910_gtu_deduplicated',
        ]

    col_bim_chr = 1
    col_bim_name = 2
    col_bim_pos = 4

    ## 2c) 
    cmd = 'cat %s.%s.concatenated' %(l_bfiles[0],l_bfiles[1],)
    cmd += ''' | awk '{if (!('''
    ## alleleA1==alleleB1 && alleleA2==alleleB2
    cmd += '''(($3==0 || $5==0 || $3==$5) && ($4==0 || $6==0 || $4==$6))'''
    cmd += ' || '
    ## alleleA1==alleleB2 && alleleA2==alleleB1
    cmd += '''(($3==0 || $6==0 || $3==$6) && ($4==0 || $5==0 || $4==$5))'''
    cmd += ''')) print $1,$2}' '''
    cmd += ' > mismatch.SNPs'
    print cmd
    os.system(cmd)

    ##
    ## 3) print columns of mismatch.SNPs to individual mismatch.SNPs files
    ##
    for i in range(len(l_bfiles)):
        bfile = l_bfiles[i]
        cmd = "cat mismatch.SNPs | awk '{print $%i}' > %s.mismatch.SNPs" %(
            i+1, l_bfiles[i],
            )
        print cmd
        os.system(cmd)
        continue

    ##
    ## 4) clean up
    ##
    os.remove('mismatch.SNPs')
    os.remove('%s.%s.concatenated' %(l_bfiles[0],l_bfiles[1],))

    return


def PLINK_exclude_duplicates():

    l_bfiles = [
        'omni2.5-8_agv_20120910_gtu',
        'omni2.5-4_20120904_agv_gtu_build37_sexupdated',
        ]

    col_bim_chr = 1
    col_bim_name = 2
    col_bim_pos = 4

    ##
    ## 1) sort
    ##
    sort(
        l_bfiles,
        'sort -k%i,%i -k%i,%i' %(
            col_bim_chr,col_bim_chr, col_bim_pos,col_bim_pos,
            ),
        )

    ##
    ## 2) identify duplicates (chromosome + position)
    ##
    for bfile in l_bfiles:
        identify_duplicates(bfile,col_bim_chr,col_bim_pos,col_bim_name,)

    ##
    ## 3) exclude duplicates
    ##
    PLINK_exclude(l_bfiles)

    ##
    ## 4) clean up
    ##
    for bfile in l_bfiles:
        os.remove('%s.bim.sorted' %(bfile))

    fd = open('SNPstats.txt','a')
    fd.write()
    fd.close()

    return

    
def PLINK_exclude(l_bfiles):

    for bfile in l_bfiles:

        cmd = 'plink \\\n'
        cmd += '--bfile %s \\\n' %(bfile)
        cmd += '--exclude %s.duplicates.SNPs \\\n' %(bfile)
        cmd += '--make-bed --out %s_deduplicated \\\n' %(bfile)
        cmd += '--noweb \\\n'
        cmd += '--allow-no-sex \\\n'
        cmd += '--nonfounders \\\n'
        cmd += '\n'
        print cmd
        os.system(cmd)

    stop

    return


def sort(l_bfiles,sort,):

    for bfile in l_bfiles:
        l_cmds = []
        l_cmds += ['cat %s.bim' %(bfile)]
        l_cmds += [sort]
##        l_cmds += ['''awk '{print $%i":"$%i,$%i}' ''' %(
##            col_bim_chr, col_bim_pos, col_bim_name,
##            )]
        cmd = ' | '.join(l_cmds)
        cmd += ' > %s.bim.sorted' %(bfile,)
        print cmd
        os.system(cmd)

    return


def identify_duplicates(bfile,col_bim_chr,col_bim_pos,col_bim_name,):

    print 'identify duplicates %s' %(bfile)
    fd_in = open('%s.bim.sorted' %(bfile,),'r')
    fd_out = open('%s.duplicates.SNPs' %(bfile,),'w')
    ## read first line
    for line in fd_in:
        line_prev = line
        break
    ## read subsequent lines
    for line in fd_in:
        l = line.split()
        l_prev = line_prev.split()
        if l[col_bim_chr-1] == '0':
            continue
        ## not a duplicate
        if not (
            l_prev[col_bim_chr-1] == l[col_bim_chr-1]
            and
            l_prev[col_bim_pos-1] == l[col_bim_pos-1]
            ):
            line_prev = line
            continue
        if (
            l[4] not in ['0','A','C','G','T',]
            or
            l[5] not in ['0','A','C','G','T',]
            ):
            name_exclude = l[col_bim_name-1]
        ## both are rs (exclude current)
        elif l[col_bim_name-1][:2] == 'rs' and l_prev[col_bim_name-1][:2] == 'rs':
            name_exclude = l[col_bim_name-1]
        ## current line is rs (exclude previous)
        elif l[col_bim_name-1][:2] == 'rs':
            name_exclude = l_prev[col_bim_name-1]
        ## prev line is rs (exclude current)
        elif l_prev[col_bim_name-1][:2] == 'rs':
            name_exclude = l[col_bim_name-1]
        ## none are rs (exclude current)
        else:
            name_exclude = l_prev[col_bim_name-1]
        fd_out.write('%s\n' %(name_exclude))
        line_prev = line
        ## continue loop over lines
        continue
    fd_in.close()
    fd_out.close()

##    cmd = 'cat %s.duplicates.SNPs %s.duplicates.SNPs | sort -u > duplicates.SNPs'

    return


def PLINK_update_map_and_sex(bfile,):

    cmd = 'plink \\\n'
    cmd += '--bfile %s \\\n' %(bfile)
    cmd += '--update-map build37.txt \\\n'
    cmd += '--extract build37.SNPs \\\n'
    cmd += '--update-sex %s_update-sex.txt \\\n' %(bfile)
##    cmd += '--make-bed --out %s_build37' %(bfile)
    cmd += '--make-bed --out %s_build37_sexupdated \\\n' %(bfile)
    cmd += '--noweb \\\n'
    cmd += '--allow-no-sex \\\n'
    cmd += '--nonfounders \\\n'
    cmd += '\n'
    print cmd
    os.system(cmd)

    return


def create_input_PLINK_update(bfile):

    create_input_PLINK_updatemap(bfile)

    create_input_PLINK_updatesex(bfile)

    return


def create_input_PLINK_updatesex(bfile):

    ##
    ## write update-sex.txt
    ##
    workbook = '/lustre/scratch107/projects/agv/data/1000genomes_all2141_150812.xlsx'
    wb = xlrd.open_workbook(workbook)
    sh = wb.sheet_by_index(0)
    l_header = sh.row_values(0)
    col_sex = l_header.index(u'Gender')
    col_iid = l_header.index(u'Broad_ID')

    fd = open('%s_update-sex.txt' %(bfile),'w')
    for rownum in range(1,sh.nrows):
        row = sh.row_values(rownum)
        sex = row[col_sex]
        iid = row[col_iid]

        ## convert float to string
        if type(iid) == float:
            iid = str(int(iid))
            pass

        fid = iid

        if sex == 'male':
            sex = 1
        elif sex == 'female':
            sex = 2
        else:
            print 'sex', sex
            sys.exit(0)

        fd.write('%s\t%s\t%s\n' %(fid,iid,sex,))
    fd.close()

    return


def create_input_PLINK_updatemap(bfile,):

    cmd = 'cat %s.liftOver.BED' %(bfile)
    cmd += ' | '
    cmd += "awk '{print $4,$2}'"
##    cmd += "awk '{print $4,$3}'"
    cmd += ' > '
    cmd += 'build37.txt'
    os.system(cmd)

    ## also create input for --extract option...
    cmd = 'cat %s.liftOver.BED' %(bfile)
    cmd += ' | '
    cmd += "awk '{print $4,$4}'"
    cmd += ' > '
    cmd += 'build37.SNPs'
    os.system(cmd)

    return


def run_liftOver(bfile,):

    cmd = './liftOver \\\n'
    cmd += '%s.BED \\\n' %(bfile)
    cmd += 'hg18ToHg19.over.chain \\\n'
    ## output
    cmd += '%s.liftOver.BED \\\n' %(bfile)
    cmd += '%s.unlifted.BED \\\n' %(bfile)
    cmd += '\n'
    os.system(cmd)

    return


def create_input_liftOver(bfile,):

    fn_out = '%s.BED' %(bfile)

    ## skip this time consuming step if file was recently generated
    ## i.e. within the last hour
    if os.path.isfile(fn_out):
        if time.time()-os.path.getmtime(fn_out) < 60*60:
            return

    cmd = 'cat %s.bim' %(bfile)
    cmd += ' | '
    cmd += """awk '{sub(/23/,"X",$1);print "chr"$1,$4,$4+1,$2}'"""
##    cmd += """awk '{sub(/23/,"X",$1);print "chr"$1,$4-1,$4,$2}'"""
    cmd += ' > '
    cmd += '%s' %(fn_out)
    os.system(cmd)

    return


def init():

    if not os.path.isfile('liftOver'):
        os.system('wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver')
        os.system('chmod +x liftOver')
    fn = 'hg18ToHg19.over.chain'
    if not os.path.isfile(fn):
        os.system('wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/%s.gz' %(fn))
        os.system('gunzip %s.gz' %(fn))

    fn = 'HumanOmni2.5-8v1_A-b37-strand.zip'
    if not os.path.isfile('HumanOmni2.5-8v1_A-b37.strand'):
        os.system('wget http://www.well.ox.ac.uk/~wrayner/strand/%s' %(fn))
        os.system('unzip %s.zip' %(fn))
        os.system('dos2unix HumanOmni2.5-8v1_A-b37.strand')
        os.system('dos2unix HumanOmni2.5-8v1_A-b37.miss')

    return



if __name__ == '__main__':
    main()
