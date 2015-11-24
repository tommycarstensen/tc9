#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, July 2013

import argparse
import os
import urllib.request
import sys
import fileinput
##import gzip

#cat joined.txt | awk 'BEGIN{above=0;max=0;maxchrom=0;maxpos=0} {print $0; if($3>0.45) {above=1; if($3>max) {max=$3;maxchrom=$1;maxpos=$2}} else {if(above==1) print maxchrom,maxpos,max; above=0;max=0}}'

fp_fai='/lustre/scratch111/resources/ref/Homo_sapiens/GRCh37_53/Homo_sapiens.GRCh37.dna.all.fa.fai'
url_ensembl='ftp://ftp.ensembl.org/pub/release-72/gtf/homo_sapiens/Homo_sapiens.GRCh37.72.gtf.gz'

def main():

    l_fp_gemma, title = argparser()

    download_ensembl()

    sort()

    gemma_ensembl_annotation(l_fp_gemma,title,)

    plot(title)

    return


def download_ensembl():

    basename = os.path.basename(url_ensembl)
    if not os.path.isfile(basename):
        urllib.request.urlretrieve(url_ensembl, basename)

    return


def sort():

    basename = os.path.basename(url_ensembl)
    if not os.path.isfile('%s.sorted' %(basename[:-3])):
        print('sort', basename)
        cmd = 'zcat %s' %(basename)
        ## autosomes
        cmd += " | awk '{if($1>=1&&$1<=22) {"
        ## print start
        cmd += ' print $1,$4,substr($16,2,length($16)-3);'
        ## print end
        cmd += ' print $1,$5,substr($16,2,length($16)-3)}'
        cmd += " }'"
        cmd += ' | sort -u | sort -k1n,1 -k2n,2 > %s.sorted' %(basename[:-3])
        os.system(cmd)

##    if bool_sorted = False:
####    basename = os.path.basename(fp_gemma)
####    print('sort', basename)
####    cmd = 'cat %s' %(fp_gemma)
####    cmd += " | awk '{print $1,$3,$10}' | sort -k1n,1 -k2n,2 "
####    cmd += ' > %s.sorted' %(basename)
####    os.system(cmd)
##    
####    print(basename)
####    l_unsorted = []
####    with gzip.open(basename) as file:
####        for line in file:
####            l = line.decode().rstrip().split()
####            chrom = l[0]
####            pos1 = l[3]
####            pos2 = l[4]
####            gene_name = l[15]
####            print(chrom,pos1,pos2,gene_name)
####            l_unsorted += [(chrom,pos1,gene_name,)]
####            l_unsorted += [(chrom,pos2,gene_name,)]
####            print(line)
####            stop

    return


def argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--title',
        dest='title',
        type=str,default='',
        required = True,
        )

    parser.add_argument(
        '--gemma',
        dest='gemma', nargs='+', type=str, required=True,
        help = '/lustre/scratch113/projects/uganda_gwas/users/sf14/imputedAnalysisWithSnpAnnotation/imputedAnalysisALLsnp/Cholesterol/output/GemmaImputedDataAssoAnalyisChrALL_Cholesterol.assoc.txt',
        )

    namespace_args = parser.parse_args()

    d = vars(namespace_args)
    l_fp_gemma = d['gemma']
    title = d['title']

    bool_found = False

    return l_fp_gemma, title


def plot(affix,):

    ## set paths
    fp_plt='%s.plt' %(affix)
##    assoc=/lustre/scratch113/projects/uganda_gwas/users/sf14/imputedAnalysisWithSnpAnnotation/imputedAnalysisALLsnp/$trait/output/GemmaImputedDataAssoAnalyisChrALL_$trait.assoc.txt
    fp_gnuplot='/nfs/team149/Software/bin/gnuplot'
    ## touch labels file in case it is not created
    os.system('touch %s.labels' %(affix))

    plt = ''
    ## terminal
    plt += "set term postscript eps enhanced color font 'Helvetica,36'\n"
    plt += "set output '%s.eps'\n" %(affix)
    plt += 'set encoding iso_8859_1\n'
    plt += 'set size 5,5\n'
    ## xtics
    plt += 'unset xtics\n'
    plt += "set xtics format ' '\n"
    plt += 'set xtics scale 0\n'
    ## xtics
    plt += 'set xtics('
    with open(fp_fai) as fd_fai:
        lines = fd_fai.readlines()
    plt += ','.join('"%i" %i.0' %(
        i+1,int(int(lines[i].split()[2])+int(lines[i].split()[1])/2)) for i in range(22))
    plt += ')\n'
    ## xlabel
##    plt += 'set xlabel "position/chromosome"\n'
    plt += 'set ylabel "-log_1_0({/Helvetica-Italic p})"\n'
    plt += 'set title "%s"\n' %(affix)
    plt += "set palette model RGB defined ( 0 '#A9A9A9', 1 '#696969', 2 '#7F7F7F', 3 '#FF0000' )\n"
    plt += 'unset colorbox\n'
    plt += 'set pm3d map\n'
    ## arrow (http://gnuplot.info/docs_4.2/gnuplot.html#x1-15800043.2)
    plt += 'set arrow from graph(0),first(8-log(5)/log(10)),4'
    plt += ' to graph(1),first(8-log(5)/log(10)),4'
    plt += ' nohead lt 0 lc rgb "black"\n'
    ## splot
    s = ','.join(lines[i].split()[2] for i in range(22))
##    plt += """splot "<cat %s | awk -v str=%s '""" %(gemma,s)
    plt += """splot "<cat %s.annotated | awk -v str=%s '""" %(affix,s)
    plt += ' BEGIN{split(str, a, \\",\\"); ymax=(8-log(5)/log(10))}'
##    plt += " {chrom=$1; pos=$3; x=pos+a[chrom]; y=-log($10)/log(10);"
    plt += " {chrom=$1; pos=$2; x=pos+a[chrom]; y=-log($3)/log(10); gene=$4;"
    plt += " if(y>ymax) {z=2+chrom%2} else {z=chrom%2};"
    plt += ' print x,y,z;'
    plt += ' if(y>ymax) {print x,y+0.50,z,gene>\\"%s.labels\\"}' %(affix)
    plt += """}'" u 1:2:3 w p pt 7 palette t ''"""
    plt += ''',"%s.labels" u 1:2:3:4 w labels font "Helvetica,18" t ''\n''' %(affix)

    with open(fp_plt,'w') as fd_plt:
        fd_plt.write(plt)

    ## gnuplot
    print('gnuplot')
    sys.stdout.flush()
    os.system('%s %s' %(fp_gnuplot,fp_plt))
##    os.remove('%s.labels' %(affix))
##    os.remove('%s.plt' %(affix))
    ## convert
    print('convert')
    sys.stdout.flush()
    os.system('convert %s.eps %s.png' %(affix,affix))
##    os.remove('%s.eps' %(affix))

    return


def gemma_ensembl_annotation(l_fp_gemma,title,):

    l_chroms = [str(chrom) for chrom in range(1,23)]

    bool_EOF_ensembl = False

    fp_ensembl = '%s.sorted' %(os.path.basename(url_ensembl)[:-3])
##    basename_gemma = os.path.basename(fp_gemma)
##    with open(fp_ensembl) as f_ensembl, open('%s.sorted' %(basename_gemma)) as f_gemma, open('%s.annotated' %(basename_gemma),'w') as f_out:
    with open(fp_ensembl) as f_ensembl, fileinput.input(l_fp_gemma) as f_gemma, open('%s.annotated' %(title),'w') as f_out:

        chrom_ensembl_prev, pos_ensembl_prev, gene_prev = next(parse_ensembl(f_ensembl))
        chrom_gemma, pos_gemma, prob = next(parse_gemma(f_gemma))
        chrom_ensembl, pos_ensembl, gene = next(parse_ensembl(f_ensembl))
        
        while True:
            ## different chromosomes
            if chrom_ensembl != chrom_gemma or chrom_ensembl_prev != chrom_gemma:
                if (
                    l_chroms.index(chrom_ensembl) < l_chroms.index(chrom_gemma)
                    or
                    l_chroms.index(chrom_ensembl_prev) < l_chroms.index(chrom_gemma)
                    ):
                    while chrom_ensembl != chrom_gemma or chrom_ensembl_prev != chrom_gemma:
                        chrom_ensembl_prev = chrom_ensembl
                        pos_ensembl_prev = pos_ensembl
                        gene_prev = gene
                        chrom_ensembl, pos_ensembl, gene = next(parse_ensembl(f_ensembl))
#                        print('B',chrom_ensembl,chrom_ensembl_prev)
                else:
##                    print('a',chrom_gemma,pos_gemma,pos_ensembl_prev,gene_prev)
                    f_out.write('%s %s %s %s\n' %(chrom_gemma,pos_gemma,prob,gene))
                    try:
                        chrom_gemma, pos_gemma, prob = next(parse_gemma(f_gemma))
                    except StopIteration:
                        break
                continue
            ## GEMMA position greater than ENSEMBL position
            if bool_EOF_ensembl == False and pos_gemma > pos_ensembl:
                chrom_ensembl_prev = chrom_ensembl
                pos_ensembl_prev = pos_ensembl
                gene_prev = gene
                while (
                    (
                        chrom_ensembl == chrom_ensembl_prev
                        and
                        pos_ensembl == pos_ensembl_prev)
#                    or
#                    chrom_ensembl != chrom_ensembl_prev
                    ):
#                    chrom_ensembl_prev = chrom_ensembl
#                    pos_ensembl_prev = pos_ensembl
#                    gene_prev = gene
                    try:
                        chrom_ensembl, pos_ensembl, gene = next(parse_ensembl(f_ensembl))
                        if pos_gemma < 190237 and pos_gemma > 190237-20000:
                            print(pos_gemma,pos_ensembl)
                    except StopIteration:
                        bool_EOF_ensembl = True
                        break
#                    print('A', chrom_ensembl, pos_ensembl, chrom_ensembl_prev, pos_ensembl_prev)
#                print(chrom_ensembl,pos_ensembl,pos_ensembl_prev)
                continue
            if chrom_ensembl != chrom_gemma:
##                print('d',chrom_gemma,pos_gemma,pos_ensembl_prev,gene_prev,chrom_ensembl,chrom_ensembl_prev)
                f_out.write('%s %s %s %s\n' %(chrom_gemma,pos_gemma,prob,gene))
            elif pos_ensembl-pos_gemma < pos_gemma-pos_ensembl_prev:
##                print('b',chrom_gemma,pos_gemma,pos_ensembl,gene,chrom_ensembl,chrom_ensembl_prev)
                f_out.write('%s %s %s %s\n' %(chrom_gemma,pos_gemma,prob,gene))
            else:
##                print('c',chrom_gemma,pos_gemma,pos_ensembl_prev,gene_prev,chrom_ensembl,chrom_ensembl_prev)
                f_out.write('%s %s %s %s\n' %(chrom_gemma,pos_gemma,prob,gene))
            try:
                chrom_gemma_prev = chrom_gemma
                pos_gemma_prev = pos_gemma
                chrom_gemma, pos_gemma, prob = next(parse_gemma(f_gemma))
                if chrom_gemma == chrom_gemma_prev and pos_gemma_prev > pos_gemma:
                    print('your input file',l_fp_gemma,'is not sorted')
                    print('prev pos',pos_gemma_prev,'curr pos',pos_gemma)
                    sys.exit()
                continue
            except StopIteration:
                break

    return


def parse_gemma(file):

    line = file.readline().rstrip()

    ## EOF
    if line == '': return

    ## skip header
    if line[:3] == 'chr':
        line = file.readline().rstrip()

    ## parse line
    l = line.split()
    prob = l[9]
    chrom = l[0]
    pos = int(l[2])

    yield chrom, pos, prob


def parse_ensembl(file):

    s = file.readline().rstrip()
    if s == '': return
    l = s.split()
    chrom = l[0]
    pos = int(l[1])
##    gene = l[15]
    gene = l[2]

    yield chrom, pos, gene


if __name__ == '__main__':
    main()
