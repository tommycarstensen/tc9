#!/software/bin/python

## Tommy Carstensen, Wellcome Trust Sanger Institute, November 2013, September 2014

import os
import argparse
import inspect
import sys
import time
sys.path.append('/nfs/users/nfs_t/tc9/github/sandbox')
import gnuplot
import subprocess

def main():

    d_arg = argparser()

    fn_strand = d_arg['strand']
    bfile = d_arg['bfile']
    fn_miss = d_arg['miss']
    fn_multiple = d_arg['multiple']
    plink = d_arg['plink']

    assert os.path.isfile(fn_strand)
    assert os.path.isfile('{}.bed'.format(bfile))
    assert os.path.isfile(fn_miss)
    assert os.path.isfile(fn_multiple)
    assert os.path.isfile(plink)

    ## sort
    sort(bfile, fn_strand)
    ## 1) rsIDs in strand but not in bim file
    strand_bim_nonintersection(bfile, fn_strand,)
    ## 2) write SNPs with mismatched chromosome and/or position (and/or allele)
    ## between strand file and bim file
    strand_bim_mismatch_position(bfile, fn_strand,)
    ## clean up
    remove_sort(bfile, fn_strand,)
    ## 3) write SNPs in miss and multiple files
    strand_miss(fn_miss,)
    ## 4) write SNPs in miss and multiple files
    strand_multiple(fn_multiple,)
    ## Venn of miss, multiple, etc. from steps 1-4
    venn(bfile,)
    ## 5) write SNPs with duplicate chromosomal positions
    bim_duplicates(bfile,)
    ## 6)
    X_XY_duplicates(bfile,)
    ## 7) remove selected samples
    ## exclude SNPs
    PLINK_remove_and_exclude_and_flip(bfile, fn_strand, plink)
    ## summary
    flip_summary(bfile)
    ## 8) QC

    return


def venn(bfile,):

    cmd = "cat %s.bim | awk '{if($1>=1&&$1<=22) print $2}' | sort" %(bfile)
    cmd += " > autosomal.SNPs"
    execmd(cmd)

    for fn_intersection in (None,'autosomal.SNPs'):

        if fn_intersection == 'autosomal.SNPs':
            suffix='%s_autosomal' %(os.path.basename(bfile))
        else:
            suffix='%s' %(os.path.basename(bfile))

        gnuplot.venn4(
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            fn1='bim_not_strand.SNPs',
            fn2='mismatches.SNPs',
            fn3='miss.SNPs',
            fn4='multiple.SNPs',
            fn_intersection=fn_intersection,
            text1='unique to bim file (not in strand file)',
            text2='position mismatches', ## incl. unplaced and allele mismatches (i.e. Ins/Del)
            text3='miss file',
            text4='multiple file',
            suffix=suffix,
            bool_remove=True,
            bool_sorted=False,
            )

    os.remove('autosomal.SNPs')

    return


def X_XY_duplicates(bfile):

    cmd = "join -1 1 -2 1 -o 2.2"
    cmd += " <(awk '$1==25{print $4}' %s.bim | sort)" %(bfile)
    cmd += " <(awk '$1==23{print $4,$2}' %s.bim | sort -k1,1)" %(bfile)
    cmd += " > X.XY.duplicates.SNPs"
    subprocess.call(cmd, shell = True, executable='/bin/bash')

    return


def flip_summary(bfile):

    basename = os.path.basename(bfile)

    ## sort bim by rsID
    cmd = '''cat %s.bim | sort -k2,2 > %s.bim.sorted''' %(bfile,basename,)
    execmd(cmd)

    cmd = 'join -1 1 -2 2 -v2 -o 2.1,2.2,2.3,2.4,2.5,2.6 exclude.SNPs %s.bim.sorted' %(basename)
    cmd += ' | sort'
    cmd += ' > %s_nonflipped.bim' %(basename)
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
            cmd = 'cat %s_%s.bim' %(basename,bim_suffix,)
            cmd += " | awk '{"
            cmd += 'if($5=="%1s"&&$6=="%1s") print' %(alleles[0],alleles[1],)
            cmd += "}'"
            cmd += ' | wc -l'
            s += ' %i' %(int(os.popen(cmd).read()))
        s += '\n'

    print(s)
    fd = open('flip_summary.txt','w')
    fd.write(s)
    fd.close()

    ##
    ## clean up
    ##
    os.remove('%s.bim.sorted' %(basename))
    os.remove('%s_nonflipped.bim' %(basename))

    return


def remove_sort(bfile,strand,):

    os.remove('%s.sorted' %(os.path.basename(strand,)))
    os.remove('%s.bim.sorted' %(os.path.basename(bfile,)))

    return


def sort(bfile,strand,):

    basename = os.path.basename(bfile)

    ## sort strand by rsID
    cmd = '''cat %s | sort -k1,1 | ''' %(strand,)
    cmd += '''awk '{sub(/X/,23,$2);sub(/Y/,24,$2)'''
    cmd += ''';sub(/XY/,25,$2);sub(/MT/,26,$2)'''
    cmd += ''';print $0}' > %s.sorted''' %(os.path.basename(strand),)
    execmd(cmd)
    ## sort bim by rsID
    cmd = '''cat %s.bim | sort -k2,2 > %s.bim.sorted''' %(
        bfile, basename,
        )
    execmd(cmd)

    return


def strand_bim_nonintersection(bfile,strand,):

    basename = os.path.basename(bfile)

    ## parse strand rsIDs
    cmd = 'cat %s.sorted' %(os.path.basename(strand))
    cmd += " | awk '{"
####    cmd += 'sub(/X/,23,$2);sub(/Y/,24,$2);sub(/XY/,25,$2);sub(/MT/,26,$2);'
####    cmd += 'print $1":"$2":"$3}'"
    cmd += 'print $1'
    cmd += "}'"
####    cmd += ' | sort'
    cmd += ' > %s.SNPs' %(os.path.basename(strand))
    execmd(cmd)
    ## parse bim rsIDs
    cmd = 'cat %s.bim.sorted' %(basename)
    cmd += " | awk '{"
####    cmd += 'print $2":"$1":"$4'
    cmd += 'print $2'
    cmd += "}'"
####    cmd += ' | sort'
    cmd += ' > %s.bim.SNPs' %(basename)
    execmd(cmd)
    ## bim not strand
    cmd = 'comm -23 %s.bim.SNPs %s.SNPs > bim_not_strand.SNPs' %(
        basename,os.path.basename(strand),)
    execmd(cmd)
    ## strand not bim (None!)
    cmd = 'comm -23 %s.SNPs %s.bim.SNPs > strand_not_bim.SNPs' %(
        os.path.basename(strand),basename,)
    execmd(cmd)

    os.remove('%s.SNPs' %(os.path.basename(strand)))
    os.remove('%s.bim.SNPs' %(basename))


    return


def strand_bim_mismatch_position(bfile,strand,):

    basename = os.path.basename(bfile)

    cmd = 'join -1 1 -2 2 -o 0,1.2,2.1,1.3,2.4,1.6,2.5,2.6,1.1,2.2'
    cmd += ' %s.sorted %s.bim.sorted' %(os.path.basename(strand), basename)
    # change chromosome 25 (XY, PAR) to chromosome 23 (X)
    cmd += " | awk '{sub(25,23,$3);"
    cmd += " chrom_strand=$2; chrom_bim=$3; if("
    ## rsID mismatch
    cmd += ' ($9!=$10)'
    ## chromosome/position mismatch
    cmd += ' || (chrom_strand!=chrom_bim) || ($4!=$5)'
    ## allele mismatch (i.e. insertions/deletions)
    cmd += ' || (substr($6,1,1)!=$8 && substr($6,2,1)!=$8 && !($7==0&&$8==0))'
    cmd += " ) print $1}'"
    cmd += " | sort"
    cmd += " > mismatches.SNPs"
    execmd(cmd)

    return


def PLINK_remove_and_exclude_and_flip(bfile,strand,plink):

    basename = os.path.basename(bfile)

    ## concatenate SNP exclusion lists
    cmd = 'cat'
    for prefix in [
        'strand_not_bim', ## None
        'bim_not_strand','mismatches','miss','multiple',
        'duplicates', 'X.XY.duplicates',
        ]:
        cmd += ' %s.SNPs' %(prefix)
    cmd += ' | sort -u > exclude.SNPs'
    execmd(cmd)

    ## strand flipping
    cmd = '''cat %s | awk '{if($5=="-") print $1}' | sort > flip.SNPs''' %(
        strand,)
    execmd(cmd)

    cmd = '%s \\\n' %(plink)
    cmd += '--bfile %s \\\n' %(bfile)
    cmd += '--make-bed --out %s_flipped \\\n' %(basename)
    cmd += '--noweb --allow-no-sex --nonfounders \\\n'
    ## exclude SNPs in concatenated SNP exclusion list
    cmd += '--exclude exclude.SNPs \\\n'
    cmd += '--flip flip.SNPs \\\n'
    execmd(cmd)

    ## Change Paternal ID and Maternal ID from -9 to 0.
    cmd = 'cat %s.fam' %(bfile)
    cmd += " | awk '{if($3==-9&&$4==-9) {$3=0; $4=0}; print $1,$2,$3,$4,$5,$6}"
    cmd += ' > %s.fam.tmp' %(bfile)

    os.rename('%s.fam.tmp %s.fam' %(bfile, bfile))

    return


def execmd(cmd):

    print(inspect.stack()[1][3])
    print(cmd)
    os.system(cmd)

    return


def strand_miss(miss,):

    cmd = "cat %s | awk 'NR>1{print $3}' | sort > miss.SNPs" %(miss,)
    execmd(cmd)

    return


def strand_multiple(multiple,):

    cmd = "cat %s | awk 'NR>1{print $1}' | sort > multiple.SNPs" %(multiple,)
    execmd(cmd)

    return


def bim_duplicates(bfile,):

    '''this function assumes that the only multiplets are duplets'''

    basename = os.path.basename(bfile)

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
    cmd += ' > %s.bim.sorted' %(basename,)
    execmd(cmd)

    ##
    ## exclude rsIDs in concatenated exclusion list
    ##
    cmd = 'join -1 1 -2 2 -v2 -o 2.1,2.2,2.3,2.4,2.5,2.6'
    cmd += ' preduplicate_exclude.SNPs %s.bim.sorted' %(basename)
    cmd += ' > %s.bim.sorted.joined' %(basename)
    execmd(cmd)

    ##
    ## sort bim by chromosome and position
    ## before loop over chromosome and position
    ##
    cmd = 'cat %s.bim.sorted.joined | sort -k1,1 -k4,4' %(basename)
    cmd += ' > %s.bim.sorted.joined.sorted' %(basename,)
    execmd(cmd)

    fd_in = open('%s.bim.sorted.joined.sorted' %(basename),'r')
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
            print(line)
            stop
        line_prev = line
    fd_in.close()
    fd_out.close()

    os.remove('%s.bim.sorted' %(basename))
    os.remove('%s.bim.sorted.joined' %(basename))
    os.remove('%s.bim.sorted.joined.sorted' %(basename))

    execmd('sort duplicates.SNPs -o duplicates.SNPs')

    return


def argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument('--bfile', required=True)
    parser.add_argument('--strand', required=True, help='')
    parser.add_argument('--multiple', required=True,)
    parser.add_argument('--miss', required=True,)
    parser.add_argument('--plink', required=True,)
    ## http://docs.python.org/3/library/functions.html#vars
    d_arg = vars(parser.parse_args())

    return d_arg


if __name__ == '__main__':
    main()
