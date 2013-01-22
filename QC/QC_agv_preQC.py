#!/software/bin/python
# -*- coding: utf-8 -*-

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012-2013

## built-ins
import sys, os, time, inspect
from xml.dom import minidom
## add-ons
import numpy
sys.path.append('/nfs/users/nfs_t/tc9/lib/python2.7/site-packages')
import xlrd
sys.path.append('/nfs/users/nfs_t/tc9/github/sandbox')
import gnuplot
##sys.path.append('/nfs/users/nfs_t/tc9/github/sandbox/math')
##import statistics


hours = 24*90

col_bim_name = col_bim_rsID = 2
col_bim_chr = 1
col_bim_pos = 4
col_strand_name = col_strand_rsID = 1
col_strand_chr = 2
col_strand_pos = 3
col_strand_alleles = 6


def main():

##    for country in root.findall('country'):
##          rank = country.find('rank').text
##          name = country.get('name')
##          print name, rank

    ## 0a) check that he have the liftOver executable and the map.chain file
    (
        l_bfiles, l_strands, l_miss, l_multiple,
        l_populations, d_populations,
        ) = init()

    ## 0b) update quad sex and remove quad/octo samples without consent
    ## i.e. run PLINK --update-sex --remove
    l_bfiles = PLINK_update_sex_remove_consent(l_bfiles)

    ##
    ## strands
    ##

    ## 1) update positions in quad strand file from build 36 to build 37
    ## *before* deduplication based on position
    l_strands[0] = update_quad_strand_to_build37(l_strands[0])

    ## 2)  exclude SNPs in miss/multiple files from strand files
    ## *before* deduplication based on position
    l_strands = exclude_strand_miss_and_multiple(l_strands,l_miss,l_multiple,)

    ## 3) exclude chr+pos duplicates within strand files
    ## *before* join based on position
    ## *after* updating positions and excluding miss/multiple SNPs
    l_strands = exclude_strand_duplicates(l_strands)

    ## 4) join strand files on chr+pos+allele
    ## *after* deduplication
    l_strands_venn = list(l_strands)
    l_strands = exclude_strand_noncommon(l_strands)

    ## 5) flip strand files
    l_strands = flip_strands(l_strands)

    ## 6) rename SNPs in quad strand file and chromosomes in both strand files
    l_strands_not_renamed = list(l_strands)
    l_strands = rename_strand(l_strands,)


    ##
    ## genotype data files (bed,bim,fam)
    ##

    ## 7) remove rsIDs not in strand files from bim files
    ## N.B. removal before name update and flip, otherwise all SNPs not renamed/flipped and renamed/flipped suffix will be misleading
    ## N.B. also risk of introducing duplicates during renaming (e.g. rs671483, rs166815, rs2019, rs3176785, rs2297476) and there is no deduplication step!
    l_bfiles = PLINK_exclude_rsIDs_not_in_strand_file(l_bfiles, l_strands_not_renamed,)

    ## 8) update positions in quad bim file from build 36 to build 37
    ## i.e. run PLINK --update-map
    ## N.B. chromosomes are not changed by liftOver, so they do not need to be updated
    l_bfiles[0] = PLINK_update_map_pos(l_bfiles[0],)

    ## 9) rename quad SNPs (update rsIDs from build 36 to build 37)
    ## i.e. run PLINK --update-map --update-name
    l_bfiles[0] = PLINK_update_map_name(l_bfiles[0], l_strands_not_renamed,)

    ## 10) flip bim files according to build 36 and build 37 strand files
    ## *after* changing rsIDs if harmonised/renamed quad strand file used
    l_bfiles = PLINK_flip(l_bfiles,l_strands,)

    ## 11a) find SNPs common between strand/bim files
    ## 11b) find SNPs common between bim/bim files
    ## 11c) and extract
    l_bfiles_venn = list(l_bfiles)
    l_bfiles, l_strands = PLINK_extract_intersection(l_strands, l_bfiles)


    ##
    ## summarize
    ##

    ## ) count SNPs at each step
    count_SNPs(l_bfiles,l_strands,)

##    ## ) create 4 set Venn
##    venn(l_strands_venn,l_bfiles_venn,)


    ##
    ## split by population
    ##

    ## 12) split by population (previously QC_agv_preQC_split_populations.py)
    ## i.e. run PLINK --keep
    PLINK_split_by_population(l_bfiles,l_populations,d_populations,)

    ## 13) merge selected quad and octo samples (previously QC_agv_preQC_merge_chips.py)
    ## i.e. run PLINK --bmerge
    PLINK_bmerge(l_bfiles,l_populations,d_populations,)

    ## 14) run sample and SNP QC

    return


def PLINK_exclude_rsIDs_not_in_strand_file(l_bfiles,l_strands,):

    suffix = 'commonID'

    l_bfiles_out, bool_exists = check_file_out(l_bfiles,suffix)
    if bool_exists == True:
        return l_bfiles_out

    for i in xrange(2):

        bfile = l_bfiles[i]
        strand = l_strands[i]

        bfile_out = bfile+'_'+suffix

        ## 2) only SNPs common to .strand and .bim file to be extracted
        ## 2a) sort .strand file
        cmd = '''cat %s.strand | sort -k1,1 > %s.strand.sorted''' %(
            strand,strand,)
        execmd(cmd)
        ## 2a) sort .bim file
        cmd = '''cat %s.bim | sort -k2,2 > %s.bim.sorted''' %(bfile,bfile,)
        execmd(cmd)
        ## 2c) join sorted build 37 .strand and .bim file on rsID
        fn_extract = 'join.%s.%s.SNPs' %(strand,bfile,)
        cmd = 'join -1 1 -2 2 -o 0 %s.strand.sorted %s.bim.sorted' %(
            strand,bfile,)
        cmd += ' > %s' %(fn_extract)
        execmd(cmd)

        cmd = 'plink \\\n'
        cmd += '--bfile %s \\\n' %(bfile,)
        cmd += '--make-bed --out %s \\\n' %(bfile_out,)
        cmd += '--extract %s \\\n' %(fn_extract)
        cmd += '--noweb --nonfounders --allow-no-sex \\\n'
        execmd(cmd)
    
    return l_bfiles_out


def rename_strand(l_strands,):

    suffix = 'renamed'
    l_strands_out, bool_exists = check_file_out(
        l_strands,suffix,extension='strand',)
    if bool_exists == True:
        return l_strands_out

    ## 1) rename quad SNPs and chromosomes
    cmd = 'paste %s.strand %s.strand' %(l_strands[0],l_strands[1],)
    cmd += " | awk '{"
    cmd += 'sub(/X/,23,$2);sub(/Y/,24,$2);sub(/XY/,25,$2);sub(/MT/,26,$2);'
    cmd += 'print $7,$2,$3,$4,$5,$6'
    cmd += "}'"
    cmd += ' > %s.strand' %(l_strands_out[0])
    execmd(cmd)

    ## 2) rename octo chromosomes
    cmd = 'cat %s.strand' %(l_strands[1],)
    cmd += " | awk '{"
    cmd += 'sub(/X/,23,$2);sub(/Y/,24,$2);sub(/XY/,25,$2);sub(/MT/,26,$2);'
    cmd += 'print $0'
    cmd += "}'"
    cmd += ' > %s.strand' %(l_strands_out[1])
    execmd(cmd)

    return l_strands_out


def flip_strands(l_strands_in,):

    suffix = 'flipped'
    l_strands_out, bool_exists = check_file_out(
        l_strands_in,suffix,extension='strand',)
    if bool_exists == True:
        return l_strands_out

    for i in xrange(2):
        strand_in = l_strands_in[i]
        strand_out = l_strands_out[i]
        ## init awk
        cmd = "cat %s.strand | awk '{" %(strand_in)
        cmd += ' if($5=="+") {print $0}'
        ## AG
        cmd += ' else if($5=="-" && $6=="AG") {print $1,$2,$3,$4,$5,"TC"}'
        ## AC
        cmd += ' else if($5=="-" && $6=="AC") {print $1,$2,$3,$4,$5,"TG"}'
        ## AT and CG
        cmd += ' else {print $0}'
        ## term awk
        cmd += "}'"
        cmd += ' > %s.strand' %(strand_out)
        execmd(cmd)

    return l_strands_out


def venn(l_strands,l_bfiles,):

    if os.path.isfile('venn4_avg.png'):
        if time.time()-os.path.getmtime('venn4_avg.png')<hours*60*60:
            return

##    ##
##    ## 1) sort
##    ##
    for strand in l_strands:
        cmd = 'cat %s.strand' %(strand)
##        cmd += ''' | awk '{print $0,$%i":"$%i}' ''' %(
##            col_strand_chr,col_strand_pos,)
        cmd += ''' | awk '{'''
        cmd += 'sub(/X/,23,$2);sub(/Y/,24,$2);sub(/XY/,25,$2);sub(/MT/,26,$2);'
        cmd += 'print $%i":"$%i":"$%i, $%i":"$%i' %(
            col_strand_chr,col_strand_pos,col_strand_alleles,
            col_strand_chr,col_strand_pos)
        cmd += "}'"
        cmd += ' | sort -k1,1 '
        cmd += ' > %s.strand.sorted' %(strand)
        execmd(cmd)
        continue

##    ##
##    ## 2) comm
##    ##
    for strand in l_strands:
        cmd = "comm -12 %s.strand.sorted %s.strand.sorted | awk '{print $2}' | sort > %s.strand.comm12" %(
            l_strands[0],l_strands[1],strand,)
        execmd(cmd)
    cmd = "comm -23 %s.strand.sorted %s.strand.sorted | awk '{print $2}' | sort > %s.strand.comm23" %(
        l_strands[0],l_strands[1],l_strands[0],)
    execmd(cmd)
    cmd = "comm -23 %s.strand.sorted %s.strand.sorted | awk '{print $2}' | sort > %s.strand.comm23" %(
        l_strands[1],l_strands[0],l_strands[1],)
    execmd(cmd)
    for i in xrange(2):
        strand = l_strands[i]
        cmd = "cat %s.strand.comm12 | awk '{print $1,0}' > %s.strand.venn" %(strand,strand,)
        execmd(cmd)
        cmd = "cat %s.strand.comm23 | awk '{print $1,%i}' >> %s.strand.venn" %(strand,i+1,strand,)
        execmd(cmd)
        cmd = 'sort %s.strand.venn -o %s.strand.venn' %(strand,strand,)
        execmd(cmd)

    ##
    ## 3) comm
    ##
    for i in xrange(2):
        bfile = l_bfiles[i]
        strand = l_strands[i]
        cmd = "comm -12 %s.strand.comm12 %s.bim.positions | awk '{print $1,0}' > %s.bim.positions.venn" %(
            strand,bfile,bfile,)
        execmd(cmd)
        cmd = "comm -13 %s.strand.comm12 %s.bim.positions | awk '{print $1,%i}' >> %s.bim.positions.venn" %(
            strand,bfile,i+1,bfile,)
        execmd(cmd)
        cmd = 'sort %s.bim.positions.venn -o %s.bim.positions.venn' %(bfile,bfile,)
        execmd(cmd)

    ##
    ## 4) create 4 set Venn diagram
    ##
    fn1='%s.strand.venn' %(l_strands[0])
    fn2='%s.strand.venn' %(l_strands[1])
    fn3='%s.bim.positions.venn' %(l_bfiles[0])
    fn4='%s.bim.positions.venn' %(l_bfiles[1])
    text1='quad.strand'
    text2='octo.strand'
    text3='quad.bim'
    text4='octo.bim'
    gnuplot.venn4(
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        fn1=fn1,
        fn2=fn2,
        fn3=fn3,
        fn4=fn4,
        text1=text1,
        text2=text2,
        text3=text3,
        text4=text4,
        suffix='agv',
        )

    return


def exclude_strand_miss_and_multiple(l_strands_in,l_miss,l_multiple,):

    suffix = 'exclmissmult'

    l_strands_out, bool_exists = check_file_out(
        l_strands_in,suffix,extension='strand',)
    if bool_exists == True:
        return l_strands_out

    for i in xrange(2):

        strand_in = l_strands_in[i]
        strand_out = strand_in+'_'+suffix

        miss = l_miss[i]
        multiple = l_multiple[i]

        ##
        ## 1) sort strand by rsID
        ##
        cmd = "cat %s.strand | awk '{print $1}' | sort -k1,1 > %s.strand.SNPs" %(strand_in,strand_in,)
        execmd(cmd)

        ##
        ## 2a) parse miss rsIDs
        ##
        cmd = "cat %s | awk '{print $3}' > %s.SNPs" %(miss,miss,)
        execmd(cmd)

        ##
        ## 2b) parse multiple rsIDs
        ##
        cmd = "cat %s | awk '{print $1}' > %s.SNPs" %(multiple,multiple,)
        execmd(cmd)

        ##
        ## 2c) concatenate miss and multiple rsIDs
        ##
        cmd = 'cat %s.SNPs %s.SNPs | sort > %s.%s.SNPs' %(
            miss,multiple,miss,multiple,)
        execmd(cmd)

        ##
        ## 2d) comm
        ##
        cmd = 'comm -23 %s.strand.SNPs %s.%s.SNPs > %s.strand.comm.SNPs' %(
            strand_in,miss,multiple,strand_in,)
        execmd(cmd)

        ##
        ## 3) sort strand by rsID
        ##
        cmd = "cat %s.strand | sort -k1,1 > %s.strand.sorted" %(strand_in,strand_in,)
        execmd(cmd)

        ##
        ## 4) join
        ##
        cmd = 'join -1 1 -2 1 -o 2.1,2.2,2.3,2.4,2.5,2.6 '
        cmd += ' %s.strand.comm.SNPs %s.strand.sorted ' %(strand_in,strand_in,)
        cmd += ' > %s.strand' %(strand_out)
        execmd(cmd)

    return l_strands_out


def PLINK_extract_intersection(l_strands,l_bfiles):

    suffix = 'commonpos'

    l_bfiles_out, bool_exists = check_file_out(l_bfiles,suffix)
    l_strands_out = [strand+'_'+suffix for strand in l_strands]
    if bool_exists == True:
        return l_bfiles_out,l_strands_out

    identify_common_strand_bim_and_bim_bim(l_strands,l_bfiles,)

    for i in xrange(2):
        strand = l_strands[i]
        strand_out = strand+'_'+suffix
        bfile = l_bfiles[i]

        ##
        ## update harmonised strand files to only contain SNPs common to
        ## both strand files and both bim files
        ##
        cmd = 'sort -k1,1 %s.strand > %s.strand.sorted' %(strand,strand,)
        execmd(cmd)
        cmd = 'sort -k1,1 %s.common.SNPs > %s.common.SNPs.sorted' %(bfile,bfile,)
        execmd(cmd)
        cmd = 'join -1 1 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6'
        cmd += ' %s.strand.sorted %s.common.SNPs.sorted' %(strand,bfile,)
        cmd += ' | sort -n -k2,2 -k3,3'
        cmd += ' > %s.strand' %(strand_out)
        execmd(cmd)

        ## extract positions that are shared
        ## between both strand files and both bim files
        bfile_out = bfile+'_'+suffix
        cmd = 'plink \\\n'
        cmd += '--bfile %s \\\n' %(bfile)
        cmd += '--extract %s.common.SNPs \\\n' %(bfile)
        cmd += '--make-bed --out %s \\\n' %(bfile_out,)
        cmd += '--noweb --allow-no-sex --nonfounders \\\n'
        execmd(cmd)

    return l_bfiles_out, l_strands_out


def identify_common_strand_bim_and_bim_bim(l_strands,l_bfiles,):

    ## also include rsID in comparison, as different duplicates might have been removed
    ## e.g. SNP6-153357967/kgp17184978 (b37 pos 6:153316274) instead of rs6553229 in the strand files
    ## only kgp17184978 (b37 pos 6:153316274) present in the bim file
    ## no, instead exclude deduplicated bim SNPs from strand

    for i in xrange(2):
        bfile = l_bfiles[i]
        strand = l_strands[i]
        ## 1) sort strand
        cmd = 'cat %s.strand ' %(strand)
        cmd += " | awk '{"
        cmd += 'sub(/X/,23,$2);sub(/Y/,24,$2);sub(/XY/,25,$2);sub(/MT/,26,$2);'
##        cmd += 'print $2":"$3'
        cmd += 'print $2":"$3":"$1'
        cmd += "}'"
        cmd += ' | sort '
        cmd += ' > %s.strand.positions' %(strand,)
        execmd(cmd)
        ## 2) sort bim
##        cmd = '''cat %s.bim | awk '{print $1":"$4}' | sort > %s.bim.positions''' %(bfile,bfile,)
        cmd = '''cat %s.bim | awk '{print $1":"$4":"$2}' | sort > %s.bim.positions''' %(bfile,bfile,)
        execmd(cmd)
        ## 3) comm between bim and strand
        cmd = 'comm -12 %s.bim.positions %s.strand.positions ' %(bfile,strand,)
        cmd += ' > %s.bim.%s.strand.common.positions' %(bfile,strand,)
        execmd(cmd)
    ## 4) comm between common bim/strand
    cmd = 'comm -12 '
    cmd += ' %s.bim.%s.strand.common.positions ' %(l_bfiles[0],l_strands[0],)
    cmd += ' %s.bim.%s.strand.common.positions ' %(l_bfiles[1],l_strands[1],)
    cmd += ' > bim.strand.common.positions'
    execmd(cmd)

    for i in xrange(2):
        bfile = l_bfiles[i]
        ## 5) sort bim
##        cmd = '''cat %s.bim | awk '{print $0,$1":"$4}' | sort -k7,7 ''' %(bfile,)
        cmd = '''cat %s.bim | awk '{print $0,$1":"$4":"$2}' | sort -k7,7 ''' %(bfile,)
        cmd += ' > %s.bim.sorted' %(bfile,)
        execmd(cmd)
        ## 6) join sorted bim and sorted positions
        cmd = 'join -1 7 -2 1 -o 1.2,1.5,1.6 %s.bim.sorted bim.strand.common.positions ' %(bfile)
        cmd += ' > %s.bim.joined' %(bfile)
        execmd(cmd)

    i1 = int(os.popen('cat %s.bim.joined | wc -l' %(l_bfiles[0])).read())
    i2 = int(os.popen('cat %s.bim.joined | wc -l' %(l_bfiles[1])).read())
    if i1 != i2:
        print 'use join instead of paste'
        sys.exit(0)

    ##
    ## 7) additional exclusion of bim/bim allele mismatches (i.e. kgp22784170)
    ##
    cmd = 'paste %s.bim.joined %s.bim.joined ' %(l_bfiles[0],l_bfiles[1],)
    ## init awk
    cmd += " | awk '{"
    ##
    cmd += "if($2==0||$5==0||(($2==$5&&$3==$6)||($2==$6&&$3==$5)))"
##    cmd += "if((($2==0||$5==0)&&$3==$6)||(($2==$5&&$3==$6)||($2==$6&&$3==$5)))"
##        cmd += 'print $%i > "%s.common.SNPs"' %(3*i+1,bfile,)
    cmd += 'print $1,$4'
    ## term awk
    cmd += "}'"
    cmd += ' > common.SNPs.pasted'
    execmd(cmd)

    for i in xrange(2):
        bfile = l_bfiles[i]
        cmd = "cat common.SNPs.pasted | awk '{print $%i}' > %s.common.SNPs" %(
            i+1,bfile,)
        execmd(cmd)

##    kgp22784170
##    join -1 7 -2 1 -o 1.2,1.5,1.6 omni2.5-4_20120904_agv_gtu_sexupdated_flipped_build37_exclmissmult_deduplicated.bim.sorted bim.strand.common.positions  > out4
##    join -1 7 -2 1 -o 1.2,1.5,1.6 omni2.5-8_agv_20120910_gtu_flipped_exclmissmult_deduplicated.bim.sorted bim.strand.common.positions  > out8
##    paste out4 out8 | awk '{if($2!=0&&$5!=0&&!(($2==$5&&$3==$6)||($2==$6&&$3==$5))) print}'
##    paste out4 out8 | awk '{if(!(($3==$6)||($3==$5)||($2==$6))) print}' | wc -l
##        paste out4 out8 | awk '{if(!(($2==0&&$3==0)||($5==0&&$6==0))&&!(($3==$6)||($3==$5)||($2==$6))) print}'
##        paste out4 out8 | awk '{if(!(($2==$5&&$3==$6)||($2==$6&&$3==$5))&&!(!(($2==0&&$3==0)||($5==0&&$6==0))&&!(($3==$6)||($3==$5)||($2==$6)))) print}'

    ## count after step 4
    l_fn = []
    for i in xrange(2):
        bfile = l_bfiles[i]
        strand = l_strands[i]
        l_fn += [
            '%s.bim.positions' %(bfile,),
            '%s.bim.%s.strand.common.positions' %(bfile,strand,),
            ]
    l_fn += ['bim.strand.common.positions',]
    for fn in l_fn:
        for setname,condition in [
            ['all','',],
            ['autosomal','if($1>=1&&$1<=22)',],
            ['X','if($1==23)',],
            ]:
            cmd = '''cat %s | awk 'BEGIN{FS=":"}{%s print}' | wc -l ''' %(fn,condition,)
##            print cmd
            print setname, fn, int(os.popen(cmd).read())

    return


def exclude_strand_noncommon(l_strands_in):

    ## echo A ; cat HumanOmni2.5M-b36_build37_exclmissmult_deduplicated_common.strand | awk '{print $2,$3,$5,$6}' > out1 ; echo B ; cat HumanOmni2.5-8v1_A-b37_exclmissmult_deduplicated_common.strand | awk '{print $2,$3,$5,$6}' > out2 ; echo C ; comm -12 out1 out2 | awk '{print $1,$2,$4}' > out3 ; echo D ; cat HumanOmni2.5M-b36_build37_exclmissmult_deduplicated_common.strand | awk '{print $2,$3,$6}' > out1 ; echo E ; cat HumanOmni2.5-8v1_A-b37_exclmissmult_deduplicated_common.strand | awk '{print $2,$3,$6}' > out2 ; echo F ; comm -12 out1 out2 > out4 ; wc -l out3 out4 ; comm -13 out3 out4 > out5 ; cat out5 | awk '{if($3!="AT" && $3!="CG") print}'
    ## 7 142163044 AG +/- (142163046 in octo bim)
    ## 7 142214205 AG +/- (142214207 in octo bim)

    l_strands_out, bool_exists = check_file_out(
        l_strands_in,'common',extension='strand',)
    if bool_exists == True:
        return l_strands_out

    col_name = 1
    col_chr = 2
    col_pos = 3
    col_flip = 5
    col_alleles = 6

    ##
    ## 1) sort
    ##
    for strand_in in l_strands_in:
        ## loop
        cmd = 'cat %s.strand' %(strand_in)

##        ## print chromosome:position
##        cmd += ''' | awk '{print $0,$%i":"$%i}' ''' %(
##            col_chr,col_pos,)
##        ## print chromosome:position:allele:flip (flip is +/-)
##        cmd += ''' | awk '{print $0,$%i":"$%i":"$%i":"$%i}' ''' %(
##            col_chr,col_pos,col_alleles,col_flip,)
##        ## print chromosome:position:allele
##        cmd += ''' | awk '{print $0,$%i":"$%i":"$%i}' ''' %(
##            col_chr,col_pos,col_alleles,)
        ## conditional print because
        ## AT and CG can't be distinguished from their flipped counterparts
        ## TA and GC
        ## N.B. This code needs to be changed if flipping of strand files is to precede this step!!!
        cmd += " | awk '{"
        ## print chromosome:position:allele:flip
        cmd += ' if($6!="AT"&&$6!="CG") {print $0,$%i":"$%i":"$%i":"$%i}' %(
            col_chr,col_pos,col_alleles,col_flip,)
        ## print chromosome:position:allele
        cmd += ' else {print $0,$%i":"$%i":"$%i}' %(
            col_chr,col_pos,col_alleles,)
        cmd += "}'"

        ## sort
        cmd += ' | sort -k7,7 '
        ## redirect stream
        cmd += ' > %s.strand.sorted' %(strand_in)
        execmd(cmd)
        continue

    ##
    ## 2) join
    ##
    for i in xrange(2):
        
        j = i+1
        cmd = 'join -1 %i -2 %i' %(7,7,)
##        cmd += ' -o 1.%i,2.%i' %(col_name,col_name,)
        cmd += ' -o %i.1,%i.2,%i.3,%i.4,%i.5,%i.6' %(j,j,j,j,j,j,)
        cmd += ' %s.strand.sorted' %(l_strands_in[0])
        cmd += ' %s.strand.sorted' %(l_strands_in[1])
        cmd += ' > %s_common.strand' %(l_strands_in[i],)
##        cmd += ' | sort > tmp1_chrpos'
##        cmd += ' | sort > tmp2_chrposall'
##        cmd += ' | sort > tmp3_chrposallfli'
    ##    print os.popen(cmd).read()
        execmd(cmd)

        print os.popen('cat %s_common.strand | wc -l' %(l_strands_in[i])).read() ## tmp
    stop ## tmp

##    cat HumanOmni2.5-8v1_A-b37_deduplicated_common.strand | awk '{print $2":"$3}' | sort > tmp.strand
##    cat omni2.5-4_20120904_agv_gtu_build37.bim | awk '{print $1":"$4}' | sort > tmp4.bim; cat omni2.5-8_agv_20120910_gtu.bim | awk '{print $1":"$4}' | sort > tmp8.bim
##    comm -12 tmp8.bim tmp.strand > tmp8merged.bim; wc -l tmp8merged.bim
##    comm -12 tmp4.bim tmp.strand > tmp4merged.bim; wc -l tmp4merged.bim
##    comm -12 tmp4merged.bim tmp8merged.bim | wc -l

##    comm -13 tmp3_chrposallfli tmp2_chrposall | sort -k1,1 > tmp4
##    comm -13 tmp3_chrposallfli tmp2_chrposall | sort -k2,2 > tmp8
##    sort -k2,2 omni2.5-4_20120904_agv_gtu_build37.bim > tmp4.bim
##    sort -k2,2 omni2.5-8_agv_20120910_gtu.bim > tmp8.bim
##    join -1 1 -2 2 -o 1.1,1.2,2.5,2.6 tmp4 tmp4.bim | sort -k1,1 > tmp4.bim.joined
##    join -1 2 -2 2 -o 1.1,1.2,2.5,2.6 tmp8 tmp8.bim | sort -k1,1 > tmp8.bim.joined
##    sort -k1,1 HumanOmni2.5M-b36.strand > tmp4.strand
##    sort -k1,1 HumanOmni2.5-8v1_A-b37.strand > tmp8.strand
##    join -1 1 -2 1 -o 1.1,1.2,2.5,2.6 tmp4 tmp4.strand | sort -k1,1 > tmp4.strand.joined
##    join -1 2 -2 1 -o 1.1,1.2,2.5,2.6 tmp8 tmp8.strand | sort -k1,1 > tmp8.strand.joined
##    paste tmp4.bim.joined tmp8.bim.joined tmp4.strand.joined tmp8.strand.joined | awk '{if( $11!=$15 && ( ($3=="A"&&$4=="T"&&$7=="T"&&$8=="A") || ($3=="C"&&$4=="G"&&$7=="G"&&$8=="C") || ($3=="G"&&$4=="C"&&$7=="C"&&$8=="G") || ($3=="T"&&$4=="A"&&$7=="A"&&$8=="T") ) ) print $1,$2,$3,$4,$7,$8,$11,$12,$15,$16}'

##    ## overlap between 17900 discordant post QC SNPs and 17042/17044 AT-/AT+/CG-/CG+ SNPs
##    cat HumanOmni2.5M-b36_build37_exclmissmult_deduplicated.strand | awk '{print $0,$2":"$3":"$6}' | sort -k7,7 > out1; cat HumanOmni2.5-8v1_A-b37_exclmissmult_deduplicated.strand | awk '{print $0,$2":"$3":"$6}' | sort -k7,7 > out2; join -1 7 -2 7 -o 1.1,2.1,1.5,1.6,2.5,2.6 out1 out2 | awk '{if(($4=="AT"||$4=="CG")&&$3!=$5) print $1,$2}' | awk '{print $1; print $2}' | sort > out3; sort ../discordant_SNPQC.SNPs -o out4 ; comm -12 out3 out4 | wc -l

    return l_strands_out


def exclude_strand_duplicates(l_strands_in):

    ## otherwise SNP4-69432884, SNP4-69611990 and SNP6-9123955 will not be excluded
    ## omni2.5-8_agv_20120910_gtu.bim:4        rs4860941       0       68963570        A       G
    ## omni2.5-8_agv_20120910_gtu.bim:4        SNP4-69432884   0       69432884        A       G
    ## omni2.5-4_20120904_agv_gtu.liftOver.BED:chr4    69280974        69280975        rs4860941
    ## omni2.5-4_20120904_agv_gtu.liftOver.BED:chr4    69280974        69280975        SNP4-69432884

    col_name = 1
    col_chr = 2
    col_pos = 3
    col_alleles = 6

    suffix = 'deduplicated'

    l_strands_out, bool_exists = check_file_out(
        l_strands_in,suffix,extension='strand',)
    if bool_exists == True:
        return l_strands_out

    l_strands_out = []
    for strand_in in l_strands_in:

        strand_out = strand_in+'_'+suffix

        ##
        ## 1) sort
        ##
        cmd = 'cat %s.strand' %(strand_in)
##        cmd += " | awk '{"
##        cmd += 'print $0,$%i":"$%i' %(col_chr,col_pos,)
##        cmd += "}'"
        cmd += ' | sort -k%i,%i -k%i,%i' %(col_chr,col_chr,col_pos,col_pos,)
        cmd += ' > %s.strand.sorted' %(strand_in)
        execmd(cmd)
        print cmd

        ##
        ## 2) identify duplicates (chromosome + position)
        ##
        fd_in = open('%s.strand.sorted' %(strand_in,),'r')
        if os.path.isfile('%s.strand' %(strand_out,)):
            os.remove('%s.strand' %(strand_out,))
        fd_out = open('%s.strand' %(strand_out,),'a')
        ## read first line
        line_prev = fd_in.readline()
        ## read subsequent lines
        while True:
            line = fd_in.readline()
            if line == '': break
            l = line.split()
            l_prev = line_prev.split()
            ## not a duplicate
            if not (
                l_prev[col_chr-1] == l[col_chr-1]
                and
                l_prev[col_pos-1] == l[col_pos-1]
                ):
                fd_out.write('%s' %(line_prev))
                line_prev = line
                continue
            ## both are rs (exclude current)
            elif l[col_name-1][:2] == 'rs' and l_prev[col_name-1][:2] == 'rs':
                fd_out.write('%s' %(line_prev))
            ## current line is rs (exclude previous)
            elif l[col_name-1][:2] == 'rs':
                fd_out.write('%s' %(line))
            ## prev line is rs (exclude current)
            elif l_prev[col_name-1][:2] == 'rs':
                fd_out.write('%s' %(line_prev))
            ## none are rs (exclude current)
            else:
                fd_out.write('%s' %(line))
            ## skip a line
            line = fd_in.readline()
            line_prev = line
            ## continue loop over lines
            continue
        ## append last line
        if not (
            l_prev[col_chr-1] == l[col_chr-1]
            and
            l_prev[col_pos-1] == l[col_pos-1]
            ):
            fd_out.write('%s' %(line))
        fd_in.close()
        fd_out.close()

    return l_strands_out


def update_quad_strand_to_build37(strand_in):

    strand_out = strand_in+'_build37'
    if os.path.isfile('%s.strand' %(strand_out)):
        return strand_out

    ## 1a) create liftOver input BED file from PLINK bim file
    strand2BED(strand_in,)

    ## 1b) run liftOver
    run_liftOver(strand_in,)

    ## 2) update strand file
    ## 2a) sort liftOver.BED by rsID
    cmd = 'cat %s.liftOver.BED | sort -k4,4 > %s.liftOver.BED.sorted' %(
        strand_in,strand_in,)
    execmd(cmd)
    ## 2b) sort quad strand file by rsID
    cmd = 'cat %s.strand | sort -k1,1 > %s.strand.sorted' %(
        strand_in,strand_in,)
    execmd(cmd)
    ## 2c) join by rsID / update strand
    cmd = 'join -1 1 -2 4 -o 1.1,1.2,2.2,1.4,1.5,1.6 '
    cmd += ' %s.strand.sorted %s.liftOver.BED.sorted ' %(strand_in,strand_in,)
    cmd += ' > %s.strand' %(strand_out)
    execmd(cmd)

    cmd = "sort -k1,1 HumanOmni2.5M-b36.strand > out1 ; sort -k1,1 HumanOmni2.5M-b36_build37.strand > out2 ; join -1 1 -2 1 -o 0,1.2,2.2 out1 out2 | awk '{if($2!=$3)print}"
    if os.popen(cmd).read().strip() != '':
        print 'chromsomes were updated in strand file, so code for updating chromosomes in the bim file needs to be written'
        sys.exit(0)

    return strand_out


def execmd(cmd):

    l_indexes = []
    l_cmd = cmd.split()

    if not l_cmd[0] in [
        ## unix utilities
        'cat','mv','cp','rm','paste','comm','grep','fgrep','join','sort',
        ## other
        'plink','unix2dos','wget','gnuplot','tar','convert',
        'bsub',
        'chmod',
        ##
        './liftOver',
        ]:
        print l_cmd
        stop_unknown_cmd

    if '%' in cmd:
        print cmd
        stop

    if cmd.split()[0] in ['cat','mv','cp','rm',]:
        l_indexes = [1]
    elif cmd.split()[0] == 'paste':
        l_indexes = [1,2,]
    elif cmd.split()[0] == 'comm':
        l_indexes = [2,3,]
    elif cmd.split()[0] in ['grep','fgrep',] and '-f' in cmd.split():
        l_indexes = [l_cmd.index('-f')+1]
    for index in l_indexes:
        if not os.path.isfile(cmd.split()[index]):
            print cmd
            print 'does not exist:', cmd.split()[index]
            sys.exit(0)

    print inspect.stack()[1][3]
    print cmd
    os.system(cmd)

    return


def strand2BED(strand_in,):

    fn_out = '%s.BED' %(strand_in)

    cmd = 'cat %s.strand' %(strand_in)
    cmd += ' | '
    cmd += """awk '{print "chr"$2,$3,$3+1,$1}'"""
    cmd += ' > '
    cmd += '%s' %(fn_out)
    execmd(cmd)

    return


def check_file_out(l_bfiles,suffix,extension='bed'):

    l_prefixes_out = []
    bool_exists = True
    for bfile in l_bfiles:
        prefix_out = bfile+'_'+suffix
        l_prefixes_out += [prefix_out]
        if not os.path.isfile('%s.%s' %(prefix_out,extension,)):
            bool_exists = False
            continue
        if time.time()-os.path.getmtime('%s.%s' %(prefix_out,extension,))>hours*60*60:
            bool_exists = False
            continue

    return l_prefixes_out, bool_exists


def PLINK_flip(l_bfiles,l_strands,):

    ## *before* updating from build 36 to build 37
    ## otherwise flip difference after update from build 36 to build 37
    ## 1       52396671        95.0413223140496        -       AT
    ## 1       52624083        100     +       AT

    suffix = 'flipped'

    l_bfiles_out, bool_exists = check_file_out(l_bfiles,suffix)
    if bool_exists == True:
        return l_bfiles_out

    for i in xrange(2):

        bfile = l_bfiles[i]
        strand = l_strands[i]

        bfile_out = bfile+'_'+suffix

        ##
        ## write SNPs to be flipped
        ##
        cmd = '''cat %s.strand | awk '{if($5=="-") print $1}' > %s.flip''' %(
            strand,bfile,)
        execmd(cmd)

        ##
        ## PLINK --flip
        ##
        cmd = 'plink \\\n'
        cmd += '--bfile %s \\\n' %(bfile,)
        cmd += '--make-bed --out %s \\\n' %(bfile_out,)
        cmd += '--flip %s.flip \\\n' %(bfile)
        cmd += '--noweb --nonfounders --allow-no-sex \\\n'
        execmd(cmd)
        
        ##
        ## check that .strand and .bim files alleles match
        ##
        cmd = 'join -1 1 -2 2 -o 1.6,2.5,2.6 '
        cmd += ' %s.strand.sorted %s.bim.sorted' %(strand,bfile,)
        ## init awk
        cmd += " | awk '{"
        ## init if
        cmd += 'if ('
        ## init negat
        cmd += '!( '
        cmd += '(($2==0 || substr($1,1,1)==$2) && substr($1,2,1)==$3)'
        cmd += ' || '
        cmd += '(substr($1,1,1)==$3 && ($2==0 || substr($1,2,1)==$2))'
        cmd += ' || '
        cmd += '($2==0 && $3==0)'
        ## term negat, if, awk
        cmd += ")) print}'"
        print cmd
        if os.popen(cmd).read().strip() != '':
            stop

    for bfile_out in l_bfiles_out:
        bool_out = check_that_output_was_generated(bfile_out)
        if bool_out == False:
            print 'not generated: %s.bed' %(bfile_out)
            sys.exit(0)

    return l_bfiles_out


def count_SNPs(l_bfiles,l_strands,):

    fn = 'SNPstats.txt'

    if os.path.isfile(fn):
        if time.time()-os.path.getmtime(fn)<hours*60*60:
            return

    fd = open(fn,'w')
    fd.close()

    ## quad unlifted, autosomal or X
    cmd = '''cat omni2.5-4_20120904_agv_gtu_sexupdated_exclconsent.unlifted.BED | awk 'NR%2==0{if ($1!="chr24" && $1!="chr25" && $1!="chr26") print}' | wc -l'''
    append_stats(cmd,fn,)

    ##
    ## bim/strand
    ##
    l_fn = []
    for i in xrange(2):
        bfile = l_bfiles[i]
        bfile = bfile.replace('_commonpos','')
        strand = l_strands[i]
        strand = strand.replace('_commonpos','')
        l_fn += [
            '%s.bim.positions' %(bfile,),
            '%s.bim.%s.strand.common.positions' %(bfile,strand,),
            ]
    l_fn += ['bim.strand.common.positions',]
    for fn2 in l_fn:
        for setname,condition in [
            ['all','',],
            ['autosomal','if($1>=1&&$1<=22)',],
            ['X','if($1==23)',],
            ]:
            cmd = '''cat %s | awk 'BEGIN{FS=":"}{%s print}' | wc -l ''' %(fn2,condition,)
            append_stats(cmd,fn)

    fd = open(fn,'a')
    fd.write('\n')
    fd.close()

    ##
    ## strand
    ##
    fd = open(fn,'a')
    fd.write('\n')
    fd.close()

    l_suffixes = ['','build37','exclmissmult','deduplicated','common','flipped','renamed','commonpos',]
    l_strands = [
        'HumanOmni2.5M-b36',
        'HumanOmni2.5-8v1_A-b37',
        ]
    for suffix in l_suffixes:
        for i in xrange(2):
            if i == 1 and suffix == 'build37':
                continue
            if suffix != '':
                l_strands[i] += '_'+suffix
            if suffix == 'renamed':
                continue
            strand = l_strands[i]
            append_stats("cat %s.strand | awk '{print $2}' | wc -l" %(strand),fn)
            append_stats("cat %s.strand | awk '{if($2>=1 && $2<=22) print $2}' | wc -l" %(strand),fn)
            append_stats('''cat %s.strand | awk '{if($2=="X"||$2==23) print $2}' | wc -l''' %(strand),fn)
            append_stats('''cat %s.strand | awk '{if($2=="Y"||$2==24) print $2}' | wc -l''' %(strand),fn)

    ##
    ## bim
    ##
    fd = open(fn,'a')
    fd.write('\n')
    fd.close()

    l_bfiles = [
        'omni2.5-4_20120904_agv_gtu',
        'omni2.5-8_agv_20120910_gtu',
        ]
    for suffix in [
##        '','sexupdated','build37','flipped','exclmissmult','deduplicated','common','renamed',
        '','sexupdated','exclconsent','commonID','build37','renamed','flipped','commonpos',
        ]:
        for i in xrange(2):
            if i==1 and suffix in ['sexupdated','build37','renamed',]:
                continue
            if suffix != '':
                l_bfiles[i] += '_'+suffix
            if suffix in ['exclconsent','renamed',]:
                continue
            bfile = l_bfiles[i]
            if not os.path.isfile('%s.bim' %(bfile)):
                append_stats("echo 0",fn)
                continue
            append_stats("cat %s.bim | awk '{print $2}' | wc -l" %(bfile),fn)
            append_stats("cat %s.bim | awk '{if($1>=1 && $1<=22) print $2}' | wc -l" %(bfile),fn)
            append_stats("cat %s.bim | awk '{if($1==23) print $2}' | wc -l" %(bfile),fn)
            append_stats("cat %s.bim | awk '{if($1==24) print $2}' | wc -l" %(bfile),fn)

    return


def parse_sample_and_population_from_Excel():

    if not os.path.isdir('samples'):
        os.mkdir('samples')

    l_xlsx = [
        '/lustre/scratch107/projects/agv/data/AGV_Omni-quad_720samples-to-recall_final_150812.xlsx',
        '/lustre/scratch107/projects/agv/data/AGV_Omni-8_1465samples_to-recall_final_150812.xlsx',
        '/lustre/scratch107/projects/agv/data/1000genomes_all2141_150812.xlsx',
        'GPC_8072_trb_IDlinks_TC_021012.xlsx',
        ]

    d_group2sample = {}
##        lines_keep = []
##        lines_remove = []
    for workbook in l_xlsx:
##        print 'workbook', workbook
        if not os.path.isfile(workbook):
            continue
        wb = xlrd.open_workbook(workbook)
        sh = wb.sheet_by_index(0)
        l_header = sh.row_values(0)
        ## Uganda tribes
        if workbook == 'GPC_8072_trb_IDlinks_TC_021012.xlsx':
            col_pop = l_header.index(u'trb')
            col_sample = l_header.index(u'sangerID')
        ## agv quad/octo
        elif workbook in [
            '/lustre/scratch107/projects/agv/data/AGV_Omni-quad_720samples-to-recall_final_150812.xlsx',
            '/lustre/scratch107/projects/agv/data/AGV_Omni-8_1465samples_to-recall_final_150812.xlsx',
            ]:
            col_pop = l_header.index(u'Ethnicity')
            col_sample = l_header.index(u'Sanger Sample Name')
        ## 1000g
        elif workbook == '/lustre/scratch107/projects/agv/data/1000genomes_all2141_150812.xlsx':
            col_pop = l_header.index(u'Population')
            col_sample = l_header.index(u'Broad_ID')
        else:
            stop
        ## loop over rows
        for rownum in range(1,sh.nrows):
            row = sh.row_values(rownum)
            sample = row[col_sample]
            group = row[col_pop].replace(' ','')
            ## convert float to int
            if type(sample) == float:
                sample = int(sample)
            ## convert to str
            sample = str(sample)
            if sample == 'APP5117332':
                print workbook
            ## N.B. ERROR IN EXCEL SHEET!!!
            if sample == '324662':
                group = 'GBR'
            if str(group) == '':
                if workbook == 'GPC_8072_trb_IDlinks_TC_021012.xlsx':
                    group = 'UgandaUnknown'
                else:
                    print row
                    stop
            if not group in d_group2sample.keys():
                d_group2sample[group] = []
            ## 29 Baganda duplicate samples on omni and quad
            if not sample in d_group2sample[group]:
                d_group2sample[group] += [sample]

    ##
    ## write samples from each population/tribe to file
    ##
    for group in d_group2sample.keys():
        fp_out = 'samples/%s.samples' %(group.replace(' ',''))
        fd = open(fp_out,'w')
        fd.write('\n'.join(d_group2sample[group]))
        fd.close()

    return


def PLINK_split_by_population(l_bfiles,l_populations,d_populations,):

    parse_sample_and_population_from_Excel()

    l_bfiles_out = []

    for bfile in l_bfiles:
        for population in l_populations:

            population_out = d_populations[population]

            bfile_out = '%s_%s' %(bfile,population_out,)

            if os.path.isfile('%s.bed' %(bfile_out)):
                continue
            if os.path.isfile('../pops/%s_quad/%s_quad.bed' %(population_out,population_out)):
                continue
            if os.path.isfile('../pops/%s_octo/%s_octo.bed' %(population_out,population_out)):
                continue
            if os.path.isfile('../pops/%s/%s.bed' %(population_out,population_out)):
                if time.time()-os.path.getmtime('../pops/%s/%s.bed' %(population_out,population_out))<hours*60*60:
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
                execmd(cmd)
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
            cmd += '--noweb \\\n'
            cmd += '--nonfounders \\\n'
            cmd += '--make-bed \\\n'
            cmd += '--allow-no-sex \\\n'
            cmd += '--out %s' %(bfile_out)
            print cmd

##            print count
##            print int(os.popen('cat %s.fam | wc -l' %(bfile)).read())
            bsub = "bsub -M3000000 -R'select[mem>3000] rusage[mem=3000]' "
##            bsub += '-P agv '
            bsub += '-G agv '
            bsub += '-q normal '
            bsub += "-J'%s' " %(population_out)
            bsub += '%s' %(cmd)

            execmd(bsub)

            l_bfiles_out += [bfile_out]

            ## continue loop over populations
            continue

        ## continue loop over bfiles
        continue

    ##
    ## check that output has been generated
    ##
##    while True: ## infinite loop
    for x in xrange(8*60/5):
        bool_break = True
        for bfile_out in l_bfiles_out:
            for extension in ['bed','bim','fam',]:
                fn_out = '%s.%s' %(bfile_out,extension,)
                if not os.path.isfile(fn_out):
                    print fn_out
                    bool_break = False
                    break
                if os.path.getsize(fn_out) == 0:
                    bool_break = False
                    break
                continue
            ## break loop over bfiles
            if bool_break == False:
                break
            for extension in ['keep','log','nof',]:
                fn = '%s.%s' %(bfile_out,extension,)
                if os.path.isfile(fn):
                    execmd('mv %s trash/%s' %(fn,fn))
            continue
        ## break loop over x
        if bool_break == True:
            break
        print x
        time.sleep(5*60)
        continue

    return


def append_stats(cmd,fn):

    print cmd
    i = int(os.popen(cmd).read())
    fd = open(fn,'a')
    fd.write('%s\n' %(cmd))
    fd.write('%i\n' %(i))
    fd.close()

    return


def PLINK_bmerge(l_bfiles,l_populations,d_populations,):

##    l_bfiles = [
##        'omni2.5-8_agv_20120910_gtu_flipped_deduplicated',
##        'omni2.5-4_20120904_agv_gtu_flipped_build37_sexupdated_deduplicated_renamed',
##        ]

    ##
    ## define dic
    ##
    d_bfiles = {}
    d_bfiles[l_bfiles[0]] = 'quad'
    d_bfiles[l_bfiles[1]] = 'octo'
    l_populations_out = []

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

        population_out = d_populations[population]

        bool_both_exist = True
        for bfile in l_bfiles:
            out = '%s_%s' %(bfile,population_out,)
            if not os.path.isfile('%s.bed' %(out)):
                bool_both_exist = False
                break
            continue

        ## 1) all data present on one chip
        if bool_both_exist == False:
            dn_out = '../pops/%s' %(population_out)
            if not os.path.isdir(dn_out):
                os.mkdir(dn_out)
            for bfile in l_bfiles:
                out = '%s_%s' %(bfile,population_out,)
                if not os.path.isfile('%s.bed' %(out)):
                    continue
                for extension in ['bed','fam','bim',]:
                    execmd('mv %s.%s ../pops/%s/%s.%s' %(
                        out,extension,
                        population_out,
                        population_out,extension,
                        ))
                    continue
                break
            pass
        ## 2) keep separate
##        elif population_out in ['Baganda', 'Banyarwanda',]:
        elif True:
            for bfile in l_bfiles:
                fn1 = '%s_%s' %(bfile,population_out,)
                fn2 = '%s_%s' %(population_out,d_bfiles[bfile])
                dn_out = '../pops/%s_%s' %(population_out,d_bfiles[bfile])
                if not os.path.isdir(dn_out):
                    os.mkdir(dn_out)
                for extension in ['bed','fam','bim',]:
                    execmd('mv %s.%s ../pops/%s_%s/%s.%s' %(
                        fn1,extension,
                        population_out,d_bfiles[bfile],
                        fn2,extension,
                        ))
                    continue
                continue
            pass
        ## 3) merge
        else:

            dn_out = '../pops/%s' %(population_out)
            if not os.path.isdir(dn_out):
                os.mkdir(dn_out)

            l_prefixes = ['%s_%s' %(bfile,population_out,) for bfile in l_bfiles]

            ##
            ## merge quad and octo (--bmerge)
            ##
            l_cmd = ['plink']
            l_cmd += ['--bmerge']
            l_cmd += ['%s.bed' %(l_prefixes[0])]
            l_cmd += ['%s.bim' %(l_prefixes[0])]
            l_cmd += ['%s.fam' %(l_prefixes[0])]
            l_cmd += ['--bfile %s' %(l_prefixes[1])]
            l_cmd += ['--make-bed --out %s' %(population_out)]
            l_cmd += ['--allow-no-sex']
            l_cmd += ['--nonfounders']
            l_cmd += ['--noweb']
            cmd = ' \\\n'.join(l_cmd)
            print cmd
            for extension in ['bed','fam','bim',]:
                for prefix in l_prefixes:
                    cmd += '\n\nrm %s.%s' %(prefix,extension,)
                cmd += '\n\nmv %s.%s ../pops/%s/%s.%s' %(
                    population_out,extension,
                    population_out,
                    population_out,extension,
                    )
            fd = open('%s_merge.sh' %(population_out),'w')
            fd.write(cmd)
            fd.close()
            execmd('chmod +x %s_merge.sh' %(population_out))

            bsub = "bsub -M3000000 -R'select[mem>3000] rusage[mem=3000]' "
            bsub += '-G agv '
            bsub += '-q normal '
            bsub += "-J'%s' " %(population_out)

            execmd('%s ./%s_merge.sh' %(bsub, population_out))

            pass

        ## continue loop over populations
        continue

    ##
    ## check that output has been generated
    ##
##    while True: ## infinite loop
    for x in xrange(8*60/5):
        bool_break = True
        for population_out in l_populations_out:
            for extension in ['bed','bim','fam',]:
                fn_out = '../%s/%s.%s' %(population_out,extension,)
                if not os.path.isfile(fn_out):
                    print fn_out
                    bool_break = False
                    break
                if os.path.getsize(fn_out) == 0:
                    bool_break = False
                    break
                continue
            ## break loop over bfiles
            if bool_break == False:
                break
            continue
        ## break loop over x
        if bool_break == True:
            break
        print x
        time.sleep(5*60)
        continue

    return


def PLINK_update_map_name(bfile,l_strands,):

    suffix = 'renamed'

    ## check out
    bfile_out = bfile+'_'+suffix
    if os.path.isfile('%s.bed' %(bfile_out)):
        if time.time()-os.path.getmtime('%s.bed' %(bfile_out))<hours*60*60:
            return bfile_out

    ## 1) create PLINK name dictionary
    ## i.e. write build37.name
    cmd = 'paste %s.strand %s.strand' %(l_strands[0],l_strands[1],)
    cmd += " | awk '{print $1,$7}'"
    cmd += ' > build37name.dic'
    execmd(cmd)
    
    ## 2) update quad SNP positions from build 36 to build 37
    cmd = 'plink \\\n'
    cmd += '--bfile %s \\\n' %(bfile,)
    cmd += '--make-bed --out %s \\\n' %(bfile_out,)
    cmd += '--update-map build37name.dic --update-name \\\n'
    cmd += '--noweb --nonfounders --allow-no-sex \\\n'
    execmd(cmd)

    ## check out
    bool_out = check_that_output_was_generated(bfile_out)
    if bool_out == False:
        print 'not generated: %s.bed' %(bfile_out)
        sys.exit(0)

    return bfile_out


def PLINK_update_map_pos(bfile,):

    suffix = 'build37'

    bfile_out = bfile+'_'+suffix
    if os.path.isfile('%s.bed' %(bfile_out)):
        if time.time()-os.path.getmtime('%s.bed' %(bfile_out))<hours*60*60:
            return bfile_out

    ## 1a) create liftOver input BED file from PLINK bim file
    bim2BED(bfile,)

    ## 1b) run liftOver
    run_liftOver(bfile,)

    ## 1c) create PLINK input file from liftOver output BED file
    ## i.e. write position dictionary build37pos.dic
    cmd = 'cat %s.liftOver.BED' %(bfile)
    cmd += ' | '
    ## rsID, pos
    cmd += "awk '{print $4,$2}'"
##    cmd += "awk '{print $4,$3}'"
    cmd += ' | sort -k1,1 '
    cmd += ' > '
    cmd += 'build37pos.dic'
    execmd(cmd)

    ## also create input for --extract option...
    cmd = 'cat %s.liftOver.BED' %(bfile)
    cmd += ' | '
    cmd += "awk '{print $4}'"
    cmd += ' > '
    cmd += 'build37.SNPs'
    execmd(cmd)

##    THIS IS WRONG!!! 1:206231264, 15:85452276, 6:52112717, 9:133223308
##    ## 1) create PLINK position dictionary
##    ## i.e. write build37.pos
##    cmd = 'cat %s.strand ' %(l_strands[0],)
##    cmd += " | awk '{print $1,$3}'"
##    cmd += ' > build37.pos'
##    execmd(cmd)
    
    ## 2) update quad SNP positions from build 36 to build 37
    cmd = 'plink \\\n'
    cmd += '--bfile %s \\\n' %(bfile)
    cmd += '--update-map build37pos.dic \\\n'
    cmd += '--extract build37.SNPs \\\n'
    cmd += '--make-bed --out %s \\\n' %(bfile_out)
    cmd += '--noweb --allow-no-sex --nonfounders \\\n'
    cmd += '\n'
    execmd(cmd)

    bool_out = check_that_output_was_generated(bfile_out)
    if bool_out == False:
        print 'not generated: %s.bed' %(bfile_out)
        sys.exit(0)

    return bfile_out


def bim2BED(bfile):

    fn_out = '%s.BED' %(bfile)

    cmd = 'cat %s.bim' %(bfile)
    cmd += " | awk '{"
    cmd += 'sub(/23/,"X",$1);sub(/24/,"Y",$1);sub(/26/,"MT",$1); print "chr"$1,$4,$4+1,$2'
    cmd += "}'"
    cmd += ' > '
    cmd += '%s' %(fn_out)
    execmd(cmd)

    return


def check_that_output_was_generated(bfile):

    '''maybe we ran out of memory and the output failed to be generated'''

    bool_out = True
    if not os.path.isfile('%s.bed' %(bfile)):
        bool_out = False
    
    return bool_out


def PLINK_update_sex_remove_consent(l_bfiles):

    l_suffixes = [
        'sexupdated_exclconsent',
        'exclconsent',
        ]

    l_bfiles_out = []
    for i in xrange(2):
        bfile = l_bfiles[i]
        suffix = l_suffixes[i]
        l, bool_exists = check_file_out([bfile],suffix,)
        bfile_out = bfile+'_'+suffix
        l_bfiles_out += [bfile_out]
        if bool_exists == True:
            continue

        if i == 0:
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

        ##
        ## run PLINK
        ##
        cmd = 'plink \\\n'
        cmd += '--bfile %s \\\n' %(bfile)
        cmd += '--make-bed --out %s \\\n' %(bfile_out)
        cmd += '--noweb --allow-no-sex --nonfounders \\\n'
        ## update sex
        if i == 0:
            cmd += '--update-sex %s_update-sex.txt \\\n' %(bfile)
        ## remove samples without consent
        cmd += '--remove ../cp8_remove.fam \\\n'
        cmd += '\n'
        execmd(cmd)

        bool_out = check_that_output_was_generated(bfile_out)
        if bool_out == False:
            print 'not generated: %s.bed' %(bfile_out)
            sys.exit(0)

    return l_bfiles_out


def run_liftOver(prefix,):

    cmd = './liftOver \\\n'
    cmd += '%s.BED \\\n' %(prefix)
    cmd += 'hg18ToHg19.over.chain \\\n'
    ## output
    cmd += '%s.liftOver.BED \\\n' %(prefix)
    cmd += '%s.unlifted.BED \\\n' %(prefix)
    cmd += '\n'
    execmd(cmd)

    return


def init():

    if not 'preQC' in os.getcwd():
        print 'you are located in the wrong dir'
        sys.exit(0)

    l_bfiles = [
        'omni2.5-4_20120904_agv_gtu',
        'omni2.5-8_agv_20120910_gtu',
        ]
    for bfile in l_bfiles:
        if not os.path.isfile('%s.bed' %(bfile)):
            print 'does not exist:', bfile
            sys.exit(0)

    if not os.path.isfile('liftOver'):
        execmd('wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver')
        execmd('chmod +x liftOver')
        pass

    ##
    ## .chain file
    ##
    fn = 'hg18ToHg19.over.chain'
    if not os.path.isfile(fn):
        execmd('wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/%s.gz' %(fn))
        execmd('gunzip %s.gz' %(fn))
        pass

    ##
    ## .strand files
    ##
    l_fn_zip = [
        'HumanOmni2.5M-b36-strand.zip',
        'HumanOmni2.5-8v1_A-b37-strand.zip',
        ]
    l_strands = [
        'HumanOmni2.5M-b36',
        'HumanOmni2.5-8v1_A-b37',
        ]

    l_miss = [
        'Strand-HumanOmni2.5M-b36.miss',
        'HumanOmni2.5-8v1_A-b37.miss',
        ]

    l_multiple = [
        'Strand-HumanOmni2.5M-b36.multiple',
        'HumanOmni2.5-8v1_A-b37.multiple',
        ]
    for i in xrange(len(l_fn_zip)):
        fn_zip = l_fn_zip[i]
        fn_strand = l_strands[i]+'.strand'
        if not os.path.isfile(fn_strand):
            execmd('wget http://www.well.ox.ac.uk/~wrayner/strand/%s' %(fn_zip))
            execmd('unzip %s' %(fn_zip))
            os.remove('%s' %(fn_zip))
            execmd('dos2unix %s' %(fn_strand))
            execmd('dos2unix %s' %(l_miss[i]))
            execmd('dos2unix %s' %(l_multiple[i]))
            pass
        continue

    l_populations = [
        u'Muganda', u'KIKUYU', u'Mandinka', u'ZULU', u'KALENJIN', u'Fula', u'Murundi', u'Munyarwanda', u'Sotho', u'Wolloff', u'Jola', u'GA-ADANGBE',
        u'IBO',
        'Ethiopia', ## SOMALI,AMHARA,OROMO
        ## 1000g
        'YRI','TSI','PUR','PEL','MXL','LWK','KHV','JPT','IBS','GIH','GBR','FIN','CLM','CHS','CHB','CEU','CDX','ASW','ACB',
    ##    'MKK', ## 31 samples
    ##    'CHD', ## 1 sample
        ]

    d_populations = {
        ## Uganda (East)
        u'Muganda':'Baganda', u'Murundi':'Barundi', u'Munyarwanda':'Banyarwanda',
        ## Kenya (East)
        u'KALENJIN':'Kalenjin', u'KIKUYU':'Kikuyu',
        ## South Africa (South)
        u'ZULU':'Zulu', u'Sotho':'Sotho',
        ## Nigeria (West)
        u'IBO':'Igbo',
        ## Ghana (West)
        u'GA-ADANGBE':'Ga-Adangbe',
        ## West
        u'Mandinka':'Mandinka', u'Wolloff':'Wolof', u'Fula':'Fula', u'Jola':'Jola',
        ## Ethiopia
        'Ethiopia':'Ethiopia', ## SOMALI,AMHARA,OROMO
        ## 1000g
        'YRI': 'YRI', 'CHB': 'CHB', 'ASW': 'ASW', 'TSI': 'TSI', 'CDX': 'CDX', 'CLM': 'CLM', 'CEU': 'CEU', 'KHV': 'KHV', 'PEL': 'PEL', 'LWK': 'LWK', 'MXL': 'MXL', 'CHS': 'CHS', 'GBR': 'GBR', 'ACB': 'ACB', 'IBS': 'IBS', 'FIN': 'FIN', 'JPT': 'JPT', 'PUR': 'PUR', 'GIH': 'GIH',
        }

    return (
        l_bfiles, l_strands, l_miss, l_multiple,
        l_populations, d_populations,
        )



if __name__ == '__main__':
    main()
