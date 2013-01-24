#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import os, inspect, sys, time
sys.path.append('/nfs/users/nfs_t/tc9/github/sandbox')
import gnuplot

def main():

    fn_strand = 'HumanOmni2.5-8v1_A-b37.strand'
    bfile = 'omni2.5-8_20120809_gwa_uganda_gtu'
    fn_miss = 'HumanOmni2.5-8v1_A-b37.miss'
    fn_multiple = 'HumanOmni2.5-8v1_A-b37.multiple'
    fp1000g = '/nfs/african_diversity/data/1KG/genotypes/release20110527/working_data/omni25_b37_bed_autosomes/1kg_r20110527_omni25_b37_autosomes_flipped'

    ## sort
    sort(bfile,fn_strand)
    ## 1) rsIDs in strand but not in bim file
    strand_bim_nonintersection(bfile,fn_strand,)
    ## 2) write SNPs with mismatched chromosome and/or position (and/or allele)
    ## between strand file and bim file
    strand_bim_mismatch_position(bfile,fn_strand,)
    ## clean up
    remove_sort(bfile,fn_strand,)
    ## 3) write SNPs in miss and multiple files
    strand_miss(bfile,fn_miss,)
    ## 4) write SNPs in miss and multiple files
    strand_multiple(bfile,fn_multiple,)
    ## Venn of miss, multiple, etc. from steps 1-4
    venn(bfile,)
    ## 5) write SNPs with duplicate chromosomal positions
    bim_duplicates(bfile,)
    ## 6) remove selected samples
    ## exclude SNPs
    PLINK_remove_and_exclude_and_flip(bfile,fn_strand,)
    ## summary
    flip_summary(bfile)
    ## 7) QC
    cmd = '/software/bin/python-2.7.3 ~/github/tc9/QC/QC.py --bfile %s_flipped --project uganda_gwas' %(bfile)
    if os.path.isfile('%s_flipped.bed' %(bfile)):
        time.sleep(300)
        execmd(cmd)

##    ##
##    ## 1000g prep
##    ##
##    prep1000g(bfile,fp_1000g,)

    return


def venn(bfile,):

    cmd = "cat %s.bim | awk '{if($1>=1&&$1<=22) print $2}' | sort > autosomal.SNPs" %(bfile)
    execmd(cmd)

    gnuplot.venn4(
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        fn1='bim_not_strand.SNPs',
        fn2='mismatches.SNPs',
        fn3='miss.SNPs',
        fn4='multiple.SNPs',
        fn_intersection='autosomal.SNPs',
        text1='unique to bim file',
        text2='position mismatches', ## incl. unplaced and allele mismatches (i.e. Ins/Del)
        text3='miss file',
        text4='multiple file',
        suffix='%s' %(bfile),
        )

    os.remove('autosomal.SNPs')

    return


def flip_summary(bfile):

    ## sort bim by rsID
    cmd = '''cat %s.bim | sort -k2,2 > %s.bim.sorted''' %(bfile,bfile,)
    execmd(cmd)

    cmd = 'join -1 1 -2 2 -v2 exclude.SNPs %s.bim.sorted' %(bfile)
    cmd += ' | sort'
    cmd += ' > %s_nonflipped.bim' %(bfile)
    execmd(cmd)

    s = ''

    ##
    ## count
    ##
    for alleles in [
        ['A','C',],
        ['T','G',],
        ['A','G',],
        ['T','C',],
        ['C','A',],
        ['G','T',],
        ['G','A',],
        ['C','T',],
        ['A','T',],
        ['T','A',],
        ['C','G',],
        ['G','C',],
        ['0','A',],
        ['0','T',],
        ['0','C',],
        ['0','G',],
        ['0','0',],
        ]:
        s += '%s' %(alleles)
        for bim_suffix in ['nonflipped','flipped',]:
            cmd = 'cat %s_%s.bim' %(bfile,bim_suffix,)
            cmd += " | awk '{"
            cmd += 'if($5=="%1s"&&$6=="%1s") print' %(alleles[0],alleles[1],)
            cmd += "}'"
            cmd += ' | wc -l'
            s += ' %i' %(int(os.popen(cmd).read()))
        s += '\n'

    print s
    fd = open('flip_summary.txt','w')
    fd.write(s)
    fd.close()

    ##
    ## clean up
    ##
    os.remove('%s.bim.sorted' %(bfile))
    os.remove('%s_nonflipped.bim' %(bfile))

    return


def remove_sort(bfile,strand,):

    os.remove('%s.sorted' %(strand,))
    os.remove('%s.bim.sorted' %(bfile,))

    return


def sort(bfile,strand,):

    ## sort strand by rsID
    cmd = '''cat %s | sort -k1,1 | ''' %(strand,)
    cmd += '''awk '{sub(/X/,23,$2);sub(/Y/,24,$2)'''
    cmd += ''';sub(/XY/,25,$2);sub(/MT/,26,$2)'''
    cmd += ''';print $0}' > %s.sorted''' %(strand,)
    execmd(cmd)
    ## sort bim by rsID
    cmd = '''cat %s.bim | sort -k2,2 > %s.bim.sorted''' %(
        bfile,bfile,
        )
    execmd(cmd)

    return


def prep1000g(bfile,fp_1000g,):

    cmd = '''cat %s.bim | awk '{print $1":"$4,$5,$6}' | sort -k1,1 > tmp1;  cat %s_flipped.bim | awk '{print $1":"$4,$5,$6}' | sort -k1,1 > tmp2; join -1 1 -2 1 tmp1 tmp2 > tmp3; head tmp3''' %(fp_1000g,bfile,)
    print cmd
    stop
    execmd(cmd)
    cmd = "cat tmp3 | awk '{if($3!=$5&&$4!=$5)print}'"
    print cmd
    stop

    return


def strand_bim_nonintersection(bfile,strand,):

    ## parse strand rsIDs
    cmd = 'cat %s.sorted' %(strand)
    cmd += " | awk '{"
####    cmd += 'sub(/X/,23,$2);sub(/Y/,24,$2);sub(/XY/,25,$2);sub(/MT/,26,$2);'
####    cmd += 'print $1":"$2":"$3}'"
    cmd += 'print $1'
    cmd += "}'"
####    cmd += ' | sort'
    cmd += ' > %s.SNPs' %(strand,)
    execmd(cmd)
    ## parse bim rsIDs
    cmd = 'cat %s.bim.sorted' %(bfile)
    cmd += " | awk '{"
####    cmd += 'print $2":"$1":"$4'
    cmd += 'print $2'
    cmd += "}'"
####    cmd += ' | sort'
    cmd += ' > %s.bim.SNPs' %(bfile)
    execmd(cmd)
    ## bim not strand
    cmd = 'comm -23 %s.bim.SNPs %s.SNPs > bim_not_strand.SNPs' %(bfile,strand,)
    execmd(cmd)
    ## strand not bim (None!)
    cmd = 'comm -23 %s.SNPs %s.bim.SNPs > strand_not_bim.SNPs' %(strand,bfile,)
    execmd(cmd)

    os.remove('%s.SNPs' %(strand))
    os.remove('%s.bim.SNPs' %(bfile))


    return


def strand_bim_mismatch_position(bfile,strand,):

    cmd = 'join -1 1 -2 2 -o 0,1.2,2.1,1.3,2.4,1.6,2.5,2.6,1.1,2.2'
    cmd += ' %s.sorted %s.bim.sorted' %(strand,bfile,)
    cmd += " | awk '{if("
    ## rsID mismatch
    cmd += ' ($9!=$10)'
    ## chromosome/position mismatch
    cmd += ' || ($2!=$3) || ($4!=$5)'
    ## allele mismatch (i.e. insertions/deletions)
    cmd += ' || (substr($6,1,1)!=$8 && substr($6,2,1)!=$8 && !($7==0&&$8==0))'
    cmd += " ) print $1}'"
    cmd += " | sort"
    cmd += " > mismatches.SNPs"
    execmd(cmd)

    return


def PLINK_remove_and_exclude_and_flip(bfile,strand,):

    ## concatenate SNP exclusion lists
    cmd = 'cat'
    for prefix in [
        'strand_not_bim', ## None
        'bim_not_strand','mismatches','miss','multiple', ## union.SNPs
        'duplicates',
        ]:
        cmd += ' %s.SNPs' %(prefix)
    cmd += ' | sort -u > exclude.SNPs'
    execmd(cmd)

    ## strand flipping
    cmd = '''cat %s | awk '{if($5=="-") print $1}' | sort > flip.SNPs''' %(strand,)
    execmd(cmd)

    cmd = 'plink \\\n'
    cmd += '--bfile %s \\\n' %(bfile)
    cmd += '--make-bed --out %s_flipped \\\n' %(bfile)
    cmd += '--noweb --allow-no-sex --nonfounders \\\n'
    ## remove samples without consent
    cmd += '--remove cp8_remove.fam \\\n'
    ## exclude SNPs in concatenated SNP exclusion list
    cmd += '--exclude exclude.SNPs \\\n'
    cmd += '--flip flip.SNPs \\\n'
    execmd(cmd)

    return


def execmd(cmd):

    print inspect.stack()[1][3]
    print cmd
    os.system(cmd)

    return


def strand_miss(bfile,miss,):

    cmd = "cat %s | awk 'NR>1{print $3}' | sort > miss.SNPs" %(miss,)
    execmd(cmd)

    return


def strand_multiple(bfile,multiple,):

    cmd = "cat %s | awk 'NR>1{print $1}' > multiple.SNPs" %(multiple,)
    execmd(cmd)

    return


def bim_duplicates(bfile,):

    '''this function assumes that the only multiplets are duplets'''

    ##
    ## concatenate prior exclusion lists
    ##
    cmd = 'cat'
    for prefix in ['mismatches','miss','multiple','bim_not_strand','strand_not_bim',]:
        if os.path.getsize('%s.SNPs' %(prefix)) > 0:
            cmd += ' %s.SNPs' %(prefix)
    cmd += ' | sort -u'
    cmd += ' > preduplicate_exclude.SNPs'
    execmd(cmd)

    ##
    ## sort bim by rsID before join by rsID
    ##
    cmd = 'cat %s.bim | sort -k2,2' %(bfile)
    cmd += ' > %s.bim.sorted' %(bfile,)
    execmd(cmd)

    ##
    ## exclude rsIDs in concatenated exclusion list
    ##
    cmd = 'join -1 1 -2 2 -v2'
    cmd += ' preduplicate_exclude.SNPs %s.bim.sorted' %(bfile)
    cmd += ' > %s.bim.sorted.joined' %(bfile)
    execmd(cmd)

    ##
    ## sort bim by chromosome and position
    ## before loop over chromosome and position
    ##
    cmd = 'cat %s.bim.sorted.joined | sort -k1,1 -k4,4' %(bfile)
    cmd += ' > %s.bim.sorted.joined.sorted' %(bfile,)
    execmd(cmd)

    fd_in = open('%s.bim.sorted.joined.sorted' %(bfile),'r')
    fd_out = open('duplicates.SNPs','w')
    line_prev = 'x x x x'
    while True:
        line = fd_in.readline()
        if not line:
            break
        l = line.split()
        l_prev = line_prev.split()
        if (
            ## chr
            l_prev[0] == l[0]
            and
            ## pos
            l_prev[3] == l[3]
            and
            ## not unplaced
            l[0] != '0'
            ):
            rsID = l[1]
            rsID_prev = l_prev[1]
            if rsID[:2] == 'rs':
                fd_out.write('%s\n' %(rsID_prev))
            elif rsID_prev[:2] == 'rs':
                fd_out.write('%s\n' %(rsID))
            else:
                fd_out.write('%s\n' %(rsID))
        elif l[0] == '0':
            print line
            stop
        line_prev = line
    fd_in.close()
    fd_out.close()

    os.remove('%s.bim.sorted' %(bfile))
    os.remove('%s.bim.sorted.joined' %(bfile))
    os.remove('%s.bim.sorted.joined.sorted' %(bfile))
    execmd('cp preduplicate_exclude.SNPs union.SNPs')

    execmd('sort duplicates.SNPs -o duplicates.SNPs')

    return


if __name__ == '__main__':
    main()
