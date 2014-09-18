#!/software/bin/python

## Tommy Carstensen, Wellcome Trust Sanger Institute, 2012-2013

## built-ins
import os
import sys
import time
import math
import inspect
import argparse
import subprocess
## add-ons
import numpy
import pwd
## http://docs.scipy.org/doc/scipy/reference/stats.html
from scipy import stats

## todo 2013-01-23: get rid of any references to uganda_gwas and agv and tc9 and lustre in final version --- the importance of this was highlighted on 2013-02-08, when the script failed for Deepti, when she used the project uganda

## todo 2013-02-08: add --max-mem option to never allow script to ask for more than 1GB of memory, to make use of plentiful hosts with 1GB if few samples in pop

## todo 2013-02-08: optimize memory management - a scatter of memory vs samples did not reveal a linear relationship...

## todo 2013-02-08: add a --no-pca option (or if not EIGENSOFT executable specified perhaps)

## todo 2013-02-08: find out how much memory EIGENSOFT uses (3628MB for uganda_gwas 4788 samples)

## todo 2013-02-08: add a --no-pca-with-1000g flag

## todo 2013-02-08: add a check that the provided 1000G bim file has SNPs that match the bim file (do a join and print any diffs and do wc -l) --- if different then inform about the --no-pca-with-1000g flag

## todo 2013-11-11: get rid of python path in this script (/software/python-3.3.2/bin/python3) - instead parse from os or sys

## todo 2013-11-11: calculate heterozygosity straight from bed file instead of creating intermediate raw file (easy to fix!!!) - saves disk space and time

## todo 2013-11-11: make plink and EIGENSOFT executable ('/nfs/users/nfs_t/tc9/opt/plink-1.07-x86_64/plink') an argument

## todo 2013-11-11: check assumption that FID and IID are identical

## todo 2013-11-11: integrate alkes price long range ldregions into script: self.cluster_exclude = '/nfs/t149_influenza_exomes/working/analysis2/ldregions.txt'

## todo 2013-11-17: replace recodeA method for calculating with het with a script that reads directly from bed

## todo 2013-11-21: get rid of *all* magic numbers (and paths) - /nfs/team149/Software/bin/gnuplot

## todo 2013-11-25: make a 3d PCA plot

class main:

    ## Essential document describing the order of operations in PLINK:
    ## http://pngu.mgh.harvard.edu/~purcell/plink/flow.shtml

    '''This script assumes that FID and IID are always identical''' ## well, check it then!!!

    def main(self,):

        self.init()

        self.check_logs()

        if self.verbose == True: print('############ execute ############')
        self.plink_execution()

        if self.verbose == True: print('############ figures ############')
        self.plink_figures()

        if self.verbose == True: print('############ tables ############')
        self.plink_tables()

        if self.verbose == True: print('############ pdf ############')
        self.pdflatex()

        return


    def pdflatex(self,):

        if not os.path.isfile('%s.postQC.autosomes.bed' %(self.o)):
            return

        if not os.path.isfile('QCtemplate.tex'):
            return

        if os.path.isfile('%s.pdf' %(self.o)):
            return

        ## write tex document to string
        s_tex = self.pdflatex_document()

        ## write tex document to file
        fd = open('%s.tex' %(self.o),'w')
        fd.write(s_tex)
        fd.close()

        ## execute pdflatex
        self.execmd('pdflatex %s.tex' %(self.o))

        print('exit')
        sys.exit(0)

        ## send QC report by mail
        cmd = '(cat "Here is your %s %s QC report"; uuencode %s.pdf %s.pdf)' %(
            self.args.project, self.i, self.o, self.o,)
        cmd += ' | mail -s "%s %s"' %(self.args.project,self.i,)
        address = '%s@sanger.ac.uk' %(pwd.getpwuid(os.getuid())[0],)
        cmd += ' %s' %(address)
        self.execmd(cmd)

        os.remove('%s.log' %(self.o))

        return


    def pdflatex_document(self,):

##        self.execmd('wget http://www.tommycarstensen.com/QCtemplate.tex')
        fd = open('QCtemplate.tex','r')
        s = fd.read()
        fd.close()
##        os.remove('QCtemplate.tex')

        l_keys_tex = []
        l_vals_tex = []

        ## project description
        uid = pwd.getpwuid(os.getuid())[0]
        s = s.replace('${uid}',uid)

##        s = s.replace('${project}',self.args.project)
##        s = s.replace('${bfile_in}',self.i)
        s = s.replace('${project}',self.args.project.replace('_','\_')) ## tmp!!! underscore problem
        s = s.replace('${bfile_in}','') ## tmp!!! underscore problem

        ##
        ## data description, samples
        ##
        n_samples = str(int(os.popen('cat %s.fam | wc -l' %(self.i)).read()))
        s = s.replace('${samples_preQC}',n_samples)

        cmd = 'cat %s.fam' %(self.i)
        cmd += " | awk '"
        cmd += ' BEGIN{cntmale=0;cntfemale=0}'
        cmd += ' {'
        cmd += ' cnt++;'
        cmd += ' if($5==1) {cntmale++}'
        cmd += ' else if($5==2) {cntfemale++}'
        cmd += ' }'
        cmd += ' END{print cnt,cntmale,cntfemale}'
        cmd += "'"
        l = os.popen(cmd).read().strip().split()
        l_keys_tex += [
            'samples_preQC','samples_males_preQC','samples_females_preQC',]
        l_vals_tex += l

        ##
        ## data description, SNPs
        ##
        cmd = 'cat %s.bim' %(self.i)
        cmd += " | awk '"
        cmd += ' BEGIN{cntAUTO=0;cntX=0;cntY;cntXY=0;cntMT=0;cnt0=0}'
        cmd += ' {'
        cmd += ' cnt++;'
        cmd += ' if($1>=1&&$1<=22) {cntAUTO++}' ## magic number
        cmd += ' else if($1==23) {cntX++}'
        cmd += ' else if($1==24) {cntY++}'
        cmd += ' else if($1==25) {cntXY++}'
        cmd += ' else if($1==26) {cntMT++}'
        cmd += ' else if($1==0) {cnt0++}'
        cmd += ' }'
        cmd += ' END{print cnt,cntAUTO,cntX,cntY,cntXY,cntMT,cnt0}'
        cmd += "'"
        l = os.popen(cmd).read().strip().split()
        l_vals_tex += l
        l_keys_tex += [
            'SNPs_preQC_all','SNPs_preQC_autosomal',
            'SNPs_preQC_X','SNPs_preQC_Y','SNPs_preQC_XY','SNPs_preQC_MT',
            'SNPs_preQC_unplaced',
            ]

        ## sample call rate
        s = s.replace('${threshold_imiss}',str(100*self.args.threshold_imiss))

        cmd = 'cat %s.imiss' %(self.o)
        cmd += " | awk '"
        cmd += ' BEGIN{sum1=0;sum2=0}'
        cmd += ' NR>1{sum1+=$6;NR1++;if($6<1-%s) {NR2++;sum2+=$6}}' %(self.args.threshold_imiss)
        cmd += ' END{print 1-(sum1/NR1), 1-(sum2/NR2)}'
        cmd += "'"
        l = os.popen(cmd).read().strip().split()
        l_vals_tex += l
        l_keys_tex += ['imiss_presampleQC','imiss_postsampleQC',]

        ## heterozygosity
        s = s.replace(
            '${threshold_het_stddev}',str(self.args.threshold_het_stddev))
        l_cmds,average,stddev,het_min,het_max = self.het2stddev(
            float(self.args.threshold_imiss),bool_execute=True)
        s = s.replace('${het_num_avg}',str(average))
        s = s.replace('${het_num_stddev}',str(stddev))

        cmd = 'sort %s.imiss.samples > %s.imiss.samples.sorted' %(
            self.o,self.o,)
        self.execmd(cmd)
        cmd = 'sort %s.het.samples > %s.het.samples.sorted' %(self.o,self.o,)
        self.execmd(cmd)
        cmd = 'comm -13 %s.imiss.samples.sorted %s.het.samples.sorted | wc -l' %(
            self.o,self.o,)
        s = s.replace('${samples_remove_het}',str(int(os.popen(cmd).read())))

        ## sex check
        cmd = 'sort %s.sexcheck.samples > %s.sexcheck.samples.sorted' %(
            self.o,self.o,)
        self.execmd(cmd)
        cmd = 'cat %s.imiss.samples %s.het.samples | sort -u' %(self.o,self.o)
        cmd += ' > %s.imiss.het.samples.sorted' %(self.o)
        self.execmd(cmd)
        cmd = 'comm -13'
        cmd += ' %s.imiss.het.samples.sorted' %(self.o)
        cmd += ' %s.sexcheck.samples.sorted' %(self.o)
        cmd += ' | wc -l'
        s = s.replace('${samples_remove_sex}',str(int(os.popen(cmd).read())))

        ## IBD
        s = s.replace(
            '${threshold_IBD_postHWE}',
            str(self.args.threshold_pi_hat_max_postHWE),)
        cmd = 'cat %s.prehardy.prune.in | wc -l' %(self.o)
        s = s.replace('${SNPs_pruned}',str(int(os.popen(cmd).read())))

        ## SNP call rate
        s = s.replace('${threshold_lmiss}',str(100*self.args.threshold_lmiss),)

        ## HWE
        s = s.replace('${threshold_hwe}',str(self.args.threshold_hwe_min),)
        s = s.replace(
            '${threshold_IBD_preHWE}',
            str(self.args.threshold_pi_hat_max_preHWE),)

        ## SNP call rate, X chromosome
        cmd = "cat %s.postQC.X.fam | awk '{if($5==1) print}' | wc -l" %(self.o)
        s = s.replace('${samples_lmiss_X_males}',str(int(os.popen(cmd).read())))
        cmd = "cat %s.postQC.X.fam | awk '{if($5==2) print}' | wc -l" %(self.o)
        s = s.replace('${samples_lmiss_X_females}',str(int(os.popen(cmd).read())))

        ##
        ## various simple counts
        ##
        for suffix, var_tex in [
            ## sample QC
            ['sampleQC.samples','samples_remove_QC',],
            ## IBD (pre HWE)
            [
                'genome.%.2f.samples' %(self.args.threshold_pi_hat_max_preHWE),
                'samples_remove_IBD_preHWE',],
            ## IBD (post HWE)
            [
                'genome.%.2f.samples' %(self.args.threshold_pi_hat_max_postHWE),
                'samples_remove_IBD_postHWE',],
            ## SNP QC
            ['lmiss.hwe.SNPs','SNPs_exclude_SNPQC',],
            ## SNP call rate, autosomal
            ['lmiss.SNPs','SNPs_exclude_lmiss',],
            ## HWE, autosomal
            ['hwe.SNPs','SNPs_exclude_hwe',],
            ## samples post QC
            ['postQC.autosomes.fam','samples_postQC',],
            ## SNPs post QC, autosomal
            ['postQC.autosomes.bim','SNPs_postQC_autosomal',],
            ## sample call rate
            ['imiss.samples','samples_remove_imiss',],
            ## SNP call rate, X
            ['X.males.lmiss.SNPs','SNPs_lmissXMale',],
            ['X.females.lmiss.SNPs','SNPs_lmissXFemale',],
            ['X.lmiss.union.SNPs','SNPs_lmissXunion',],
            ## HWE, X
            ['X.females.hwe.SNPs','SNPs_exclude_X_HWE',],
            ## SNPs post QC, X
            ['postQC.X.bim','SNPs_postQC_X',],
            ]:
            cmd = 'cat %s.%s | wc -l' %(self.o,suffix)
            n = int(os.popen(cmd).read())
            s = s.replace('${%s}' %(var_tex),str(n))

        ##
        ## tables
        ##
        d_tables = {
            'imiss':'samples_imiss_tr',
            'het':'samples_het_tr',
            'sexcheck':'samples_sexcheck_tr',
            'lmiss':'SNPs_lmiss_tr',
            'genome_outliers':'samples_genome_postHWE_tr',
            'genome':'samples_genome_tr',
            'hwe':'SNPs_hwe_tr',
            'lmiss.X.males':'SNPs_lmissXMale_tr',
            'lmiss.X.females':'SNPs_lmissXFemale_tr',
            'hwe.X.females':'SNPs_hwe_X_tr',
            }
        for table_prefix,tex_key in list(d_tables.items()):
            cmd = 'cat %s.table' %(table_prefix)
            cmd += ' | cut -f2-'
            l = os.popen(cmd).read().strip().split('\n')[1:]
            for i in range(len(l)):
                l[i] = ' & '.join(l[i].replace('_','\\_').split())
            s = s.replace('${%s}' %(tex_key),' \\\\\n'.join(l))

        ##
        ## figures
        ##
        s = s.replace('${fig_imiss}','{%s.imiss}.png' %(self.o))
        s = s.replace('${fig_het}','{%s.het}.png' %(self.o))
        s = s.replace('${fig_lmiss}','{%s.lmiss}.png' %(self.o))
        s = s.replace('${fig_imiss_het}','{%s.het.call}.png' %(self.o))
        s = s.replace('${fig_genome}','{%s.genome.max}.png' %(self.o))
        s = s.replace('${fig_hwe}','{%s.hwe}.png' %(self.o))
        s = s.replace('${fig_hwehet}','{%s.hwe.scatter}.png' %(self.o))
        if os.path.isfile('%s.EIGENSOFT.evec' %(self.o)):
            s = s.replace('${fig_PCAwithin12}','{%s.EIGENSOFT.evec.pc1.pc2}.png' %(self.o))
            s = s.replace('${fig_PCAwithin34}','{%s.EIGENSOFT.evec.pc1.pc2}.png' %(self.o))
        else:
            s = s.replace('${fig_PCAwithin12}','{%s.imiss}.png' %(self.o))
            s = s.replace('${fig_PCAwithin34}','{%s.imiss}.png' %(self.o))

        for i in range(len(l_keys_tex)):
            key_tex = l_keys_tex[i]
            print(key_tex)
            print(l_vals_tex)
            print(l_keys_tex)
            s = s.replace('${%s}' %(key_tex),l_vals_tex[i])

##        if '${' in s:
##            print s
##            stop

        return s


    def gnuplot_histogram(
        self,
        ## in
        prefix_in,
        column = '$1',
        s_plot = None,
        ## out
        prefix_out = None,
        bool_remove = True,
        ## style
        color = 'blue',
        xlabel=None, ylabel=None, title=None,
        x_min=None,x_max=None,y_min=None,y_max='',
        x_step = None, ## width of boxes
        tic_step = None, ## tics on axis
        bool_timestamp = False,
        lines_extra = None,
        ):

        if not prefix_out:
            prefix_out = prefix_in

        if bool_timestamp == True:
            title += '\\n%s' %(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime()))

        sett = []
        sett += [
            'set terminal postscript eps enhanced color "Helvetica" 24\n',
            'set output "%s.ps"\n' %(prefix_out),
            'set size 3,3\n', ## scale 300%
            'set encoding iso_8859_1\n',
            'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
            ]
        if x_min and x_max:
            sett += [
                'min=%f\n' %(float(x_min)),
                'max=%f\n' %(float(x_max)),
                'set xrange [min:max]\n',
                ]
            if tic_step:
                sett += ['set xtics min,%f,max\n' %(tic_step),]
        elif x_max:
            sett += ['set xrange [:%f]\n' %(float(x_max))]

        sett += [
            'width=%s	#interval width\n' %(x_step),
            '#function used to map a value to the intervals\n',
            'hist(x,width)=width*floor(x/width)+width/2.0\n',
            'set yrange [0:]\n',
    ##        'set xtics min,(max-min)/%i,max\n' %(n_tics), ## min and max must be floats...
            'set boxwidth width*0.9\n',
            'set style fill solid 0.5	#fillstyle\n',
            'set tics out nomirror\n',
            'set title "%s" noenhanced\n' %(title),
            '#count and plot\n',
    ##        'unset ytics\n',
            ]
        if xlabel:
            sett += ['set xlabel "{/=36 %s}"\n' %(xlabel),]
        if ylabel:
            sett += ['set ylabel "{/=36 %s}"\n' %(ylabel),]
        if lines_extra:
            sett += lines_extra

        ## The frequency option makes the data monotonic in x; points with the same
        ## x-value are replaced by a single point having the summed y-values.  The
        ## resulting points are then connected by straight line segments.
        ## See also
        ## smooth.dem
        if s_plot:
            sett += [s_plot]
        else:
            s_plot = 'plot [:][:%s]"%s"' %(y_max, prefix_in)
            s_plot += ' u (hist(%s,width)):(1.0)' %(column)
            s_plot += ' smooth freq w boxes lc rgb"%s" notitle\n' %(color)
            sett += [s_plot]

        self.gnuplot(sett,prefix_out,bool_remove,)

        return


    def gnuplot(self,sett,prefix_out,bool_remove,):

        ## remove png if it exists
        if os.path.isfile('%s.png' %(prefix_out)):
            os.remove('%s.png' %(prefix_out))

        ## write gnuplot settings
        fd = open('%s.plt' %(prefix_out),'w')
        fd.writelines(sett)
        fd.close()

        ##
        ## execute gnuplot settings
        ##
##        os.system('/usr/bin/gnuplot %s.plt' %(prefix_out))
        ## set object circle only works with gnuplot 4.3 and higher
        os.system('/nfs/team149/Software/bin/gnuplot %s.plt' %(prefix_out))
        ## convert postscript to portable network graphics
        ## convert postscript to png if it was created
        if os.path.isfile('%s.ps' %(prefix_out)):
            os.system('convert %s.ps %s.png' %(prefix_out,prefix_out,))
            os.remove('%s.ps' %(prefix_out))

        ## clean up
        if bool_remove == True:
            os.remove('%s.plt' %(prefix_out))

        return


    def gnuplot_scatter(
        self,
        prefix,
        xlabel='',ylabel='',
        title=None,
        x_min = '', x_max = '',
        y_min = '', y_max = '',
        terminal = 'postscript',
        column1 = '1', column2 = '2',
        line_plot = None, ## manual plot line
        s_plot = None, ## manual def of stuff to be plotted
        bool_remove = True,
        lines_extra = None,
        prefix_out=None,
        bool_title_enhanced=True,
        bool_timestamp = False,
        ):

        if not prefix_out:
            prefix_out = prefix

        d_output = {
            'postscript':'ps',
            'png':'png',
            }
        d_terminal = {
            'postscript':'%s eps enhanced color "Helvetica" 48' %(terminal),
            'png':'png',
            }
        
        sett = []
               
        sett += [
            'set terminal %s\n' %(d_terminal[terminal]),
            'set output "%s.%s"\n' %(prefix_out, d_output[terminal]),
            'set size 4,4\n', ## scale 400%
            'set encoding iso_8859_1\n', ## postscript encoding *necessary* for special characters (e.g. Angstrom)
            'set xlabel "%s"\n' %(xlabel),
            'set ylabel "%s"\n' %(ylabel),
        ]

        if lines_extra:
            sett += lines_extra

        if title:
            if bool_timestamp == True:
                title += '\\n%s' %(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime()))
            if bool_title_enhanced == True:
                sett += ['set title "%s"\n' %(title)]
            else:
                sett += ['set title "%s" noenhanced\n' %(title)]

        if not line_plot:

            line_plot = 'plot '
            if not s_plot:
                line_plot += '[%s:%s]' %(x_min,x_max)
                line_plot += '[%s:%s]' %(y_min,y_max)
                line_plot += '"%s" u %s:%s ' %(prefix,column1,column2,)
            else:
                line_plot += s_plot

            if regression_data:
                ps = 3
            else:
                ps = 2
            line_plot += ' lc 0 lt 1 ps %s lw 3 pt 7 t ""' %(
                ps,
                )

            if errorbars == True:
                line_plot += ' w errorb'
            if regression == True:
                line_plot += ', f(x) lt 1 lc 1 lw 10'
                if regression_title:
                    line_plot += ' t "%s"' %(regression_title)
                else:
                    line_plot += ' t ""'
                    
            line_plot += '\n'

        sett += [line_plot]

        self.gnuplot(sett,prefix_out,bool_remove)

        return


    def histogram_snpweight(self,):

        '''This function is not called from anywhere at the moment...'''

        for pc in range(1,11):
            self.gnuplot_histogram(
                '%s.posthardy.EIGENSOFT.snpweight' %(self.o),
                prefix_out = '%s.posthardy.EIGENSOFT.%i.snpweight' %(self.o,pc,),
                x_step=0.05,
                column='$%i' %(3+pc),
                xlabel='snpweight',
                title='%s\\nPC %i' %(
                    self.o.replace('_','\\\\_'),
                    pc,
                    ),
                color = 'cyan',
                )

        return


    def check_logs(self,):

        print('checking logs - temporary function')

        ## this is just a temporary function...

        l = os.listdir(os.getcwd())
        for fn in l:
            if fn[-4:] != '.log':
                continue
##            if os.path.isdir(fn):
##                continue
            fd = open(fn)
            s = fd.read()
            fd.close()
##            if 'Warning' in s:
##                print s
##                stop_warning
            if 'ERROR' in s:
                print(s)
                print(fn)
####                os.remove(fn)
                stop_error_log
        l = os.listdir('%s/stderr' %(os.getcwd()))
        for fn in l:
            if fn[-4:] != '.err':
                continue
##            if os.path.isdir(fn):
##                continue
            fd = open('stderr/%s' %(fn))
            s = fd.read()
            fd.close()
##            if 'Warning' in s:
##                print s
##                stop_warning
            if 'ERROR' in s.upper():
                print(s)
                print('stderr/%s' %(fn))
####                os.remove('stderr/%s' %(fn))
                stop_error_stderr
##            if os.path.getsize('stderr/%s' %(fn)) > 0:
##                print 'stderr/%s' %(fn)
##                stop
            continue
        l = os.listdir('%s/stdout' %(os.getcwd()))
        for fn in l:
            if '.genome.out' in fn: continue
            if '.indep-pairwise.out' in fn: continue
            if fn[-4:] != '.out':
                continue
##            if os.path.isdir(fn):
##                continue
            fd = open('stdout/%s' %(fn))
            s = fd.read()
            fd.close()
##            if 'Warning' in s:
##                print s
##                stop_warning
            if 'TERM_MEMLIMIT' in s.upper():
                print(s)
                print(fn)
####                os.remove('stdout/%s' %(fn))
                stop_error_MEMLIMIT

        print('checked logs')

        return


    def scatter_hwe(self,):

        if os.path.isfile('%s.qq.hwe.png' %(self.o)):
            return
        if os.path.isfile('%s.hwe.scatter.png' %(self.o)):
            return

        n_samples = int(os.popen('cat %s.fam | wc -l' %(self.i)).read())
        cmd = 'cat %s.sampleQC.nonfounders.samples | wc -l' %(self.o)
        n_samples -= int(os.popen(cmd).read())

        n_SNPs = (int(os.popen('cat %s.hwe | wc -l' %(self.o)).read())-1)/3

##        s = '#!/software/bin/R\n'
##        s += 'df<-read.table("%s.hwe", header=T)\n' %(self.o)
##        s += 'png("%s.hwe.qqplot.png")\n' %(self.o)
##        s += 'o=-log10(sort(df$p_lrt, decreasing=F))\n'
##        s += 'e=-log10(ppoints(length(o)))\n'
##        s += 'qqplot(e, o, pch=20, cex=1, col="black", xlim=c(0, max(e)), ylim=c(0, max(o)))\n'
##        s += 'lines(e,e, col="grey")\n'
##        s += 'dev.off()\n'
##        fd = open('%s.qq.r' %(self.o),'w')
##        fd.write(s)
##        fd.close()
##        self.execmd('chmod +x %s.qq.r' %(self.o))
##        self.execmd('./%s.qq.r' %(self.o))

##        cmd = 'cat %s.hwe' %(self.o)
##        cmd += ' | '
##        cmd += "awk 'NR>1{"
##        cmd += 'if ($9!="NA") {logp=-log($9)/log(10); print logp} '
##        cmd += "}'"
##        cmd += ' | sort -r -n '
##
##        ## http://gettinggeneticsdone.blogspot.co.uk/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html
##        ## Thanks to Daniel Shriner at NHGRI for providing this code for creating expected and observed values
##        cmd += " | awk '{print $1,-log(NR/%i)/log(10)}' " %(n_SNPs,)
##
##        cmd += '> %s.qq.hwe.dat' %(self.o)
##        self.execmd(cmd)

##        from scipy.stats import norm
##        n = norm.ppdf((0.5+arange(len(p))/len(p), loc=m, scale=s)

##        if self.bool_verbose == True:
##            bool_timestamp = True
##
##        gnuplot_scatter(
##            '%s.qq.hwe' %(self.o),
####            s_plot = '"< paste %s.het %s.imiss" u (($5-$3)/$5):(1-$12)' %(self.o,self.o,),
####            line_plot = line_plot,
##            column1 = 1, column2 = 2,
##            xlabel = '-log(p_{HWE,observed})',
##            ylabel = '-log(p_{HWE,expected})',
####            xmin=0,xmax=8,
##            title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
##                self.i.replace('_','\\\\_'),n_samples,n_SNPs,
##                ),
##            bool_timestamp=bool_timestamp,
##            )

        ##
        ## 2nd plot...
        ##
        line_plot = 'plot '
        for s_ineq, suffix,color in [
            ['<','fail',1,],
            ['>','pass',0,],
            ]:
            cmd = 'cat %s.hwe' %(self.o)
            cmd += " | awk '"
            cmd += 'NR>1{if ($9!="NA" && $9 %s %.1e) print}' %(
                s_ineq,float(self.args.threshold_hwe_min))
            cmd += "'"
            cmd += ' > %s.hwe.%s' %(self.o,suffix)
            self.execmd(cmd)
            line_plot += '"%s.hwe.%s" u 7:8 ps 1 pt 7 lc %i t "",' %(
                self.o,suffix,color,)
        line_plot = line_plot[:-1]+'\n'

        if self.bool_verbose == True:
            bool_timestamp = True

        self.gnuplot_scatter(
            '%s.hwe.scatter' %(self.o),
            column1 = '$7', column2 = '$8',
            line_plot = line_plot,
            xlabel = 'O(HET)',
            ylabel = 'E(HET)',
##            xmin=0,xmax=8,
            title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                self.o.replace('_','\\\\_'),n_samples,n_SNPs,
                ),
            prefix_out = '%s.hwe.scatter' %(self.o),
            bool_timestamp=bool_timestamp,
            )

        os.remove('%s.hwe.fail' %(self.o))
        os.remove('%s.hwe.pass' %(self.o))

        return


##    def scatter_lmiss_frq(self):
##
##        if os.path.isfile('%s.lmiss.frq.png' %(self.o)):
##            return
##
##        ##
##        ## join
##        ##
##        cmd = "sort -k2,2 %s.SNPQC.lmiss > %s.SNPQC.lmiss.sorted" %(self.o,self.o,)
##        self.execmd(cmd)
##        cmd = "sort -k2,2 %s.frq > %s.frq.sorted" %(self.o,self.o,)
##        self.execmd(cmd)
##        cmd = 'join -1 2 -2 2 %s.SNPQC.lmiss.sorted %s.frq.sorted > %s.lmiss.frq.joined' %(
##            self.o,self.o,self.o,
##            )
##        self.execmd(cmd)
##
##        os.remove('%s.SNPQC.lmiss.sorted' %(self.o))
##        os.remove('%s.frq.sorted' %(self.o))
##
##        line_plot = 'plot [.9:1][0:.5]'
##        line_plot += '"%s.lmiss.frq.joined" u (1-$5):10 lc 0 ps 1 pt 7 t ""' %(self.o,)
##
##        l_arrows = [
##            'set arrow from %f,graph(0) to %f,graph(1) nohead lc 2 lt 1 lw 3 front\n' %(
##                .99,.99,
##                ),
##            'set arrow from %f,graph(0) to %f,graph(1) nohead lc 3 lt 1 lw 3 front\n' %(
##                .95,.95,
##                ),
##            'set arrow from 0,%f to 1,%f nohead lc 4 lt 2 lw 2 front\n' %(
##                .05,.05,
##                ),
##            ]
##
##        n_samples = int(os.popen('cat %s.lmiss.frq.joined | wc -l' %(self.o)).read())
##                            
##        cmd = '''cat %s.bim | awk '{if ($1<=22 && $1>=1) print}' | wc -l''' %(self.o)
##        cmd = '''cat %s.SNPQC.lmiss | wc -l''' %(self.o)
##        n_SNPs = int(os.popen(cmd).read())
##        ## subtract header
##        n_SNPs -= 1
##
##        gnuplot_scatter(
##            '%s.lmiss.frq.joined' %(self.o),
####            s_plot = '"< paste %s.lmiss %s.frq" u (1-$5):$10' %(self.o,self.o,),
####            s_plot = s_plot,
##            line_plot = line_plot,
##            column1 = '(1-$5)', column2 = 9,
##            xlabel = 'F MISS',
##            ylabel = 'MAF',
##            title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
##                self.i.replace('_','\\\\_'),n_samples,n_SNPs,
##                ),
####            ymax=.5,ymin=0,
######            xmin=0,xmax=1,
####            xmin=0.9,xmax=1,
##            lines_extra = l_arrows,
##            prefix_out='%s.lmiss.frq' %(self.o),
##            )
##
##        return


    def plink_tables(self):

        for l_suffixes_in,function in [
            [['imiss'],self.table_imiss,],
            [['het'],self.table_het,],
            [['sexcheck'],self.table_sexcheck,],
            [['postIBD.frq'],self.table_frq,],
            [['SNPQC.lmiss'],self.table_lmiss,],
            [['prehardy.genome'],self.table_genome,],
            [['hwe','X.females.hwe',],self.table_hwe,],
            ]:
            bool_continue = False
            for in_suffix in l_suffixes_in:
                fp_in = '%s.%s' %(self.o,in_suffix,)
                out_suffix = in_suffix
                fp_out = '%s.table' %(out_suffix)
                if not os.path.isfile(fp_in):
                    bool_continue = True
                    break
                if not os.path.getsize(fp_in) > 0:
                    bool_continue = True
                    break
                if os.path.isfile(fp_out):
                    bool_continue = True
                    break
                bool_input = self.check_file_in(fp_in)
                if bool_input == False:
                    bool_continue = True
                    break
                continue
            if bool_continue == True:
                continue

            ## lock
            if os.path.isfile('%s.%s.lock' %(self.o,in_suffix,)):
                sys.exit(0)
            else:
                self.execmd('touch %s.%s.lock' %(self.o,in_suffix,))
            ## tabulate
            function()
            ## unlock
            os.remove('%s.%s.lock' %(self.o,in_suffix,))


        ##
        ## combined
        ##
        print('\ncounts')
        for suffix in ['imiss','het','sexcheck',]:
            if not os.path.isfile('%s.%s.samples' %(self.o,suffix,)):
                continue
            cmd = 'cat %s.%s.samples | wc -l' %(self.o,suffix,)
            print(suffix, os.popen(cmd).read(), self.o)
        if (
            os.path.isfile('%s.imiss' %(self.o))
            and
            os.path.isfile('%s.het' %(self.o))
            and
            os.path.isfile('%s.sexcheck' %(self.o))
            ):
            set_combined = set()
            for suffix in ['imiss','het','sexcheck',]:
                l = os.popen('cat %s.%s.samples' %(self.o,suffix,)).readlines()
##                print len(l)
                set_combined |= set(l)
            print('combined', len(set_combined))

        return


    def table_frq(self):

        affix = 'postIBD'
        fn_table = fn = 'frq.%s.table' %(affix)
        
        if os.path.isfile(fn_table): return

        cmd = "awk 'NR>1{print $5}' %s.%s.frq" %(self.o,affix,)
        if self.verbose == True:
            print(cmd)
        l_frq = [float(s) for s in os.popen(cmd).readlines()]

##        l_scores = [.0,.1,.2,.3,.4,.5,1.,2.,3.,4.,5.,10.,20.,30.,40.,50.,]
        l_scores = [.0,.1,.2,.5,1.,2.,5.,10.,20.,50.,]

        s = '%s' %(self.i)
        for score in l_scores:
            score /= 100
            ## N.B. "kind" works with Python 2.7.3 but not Python 2.5.2
            n = int(len(l_frq)*stats.percentileofscore(l_frq, score, kind='weak',))/100
            if self.verbose == True:
                print(score, n)
            s += '\t%i' %(n)
        s += '\n'

        if not os.path.isfile(fn):
            l = [str(score) for score in l_scores]
            s_header = '#\t'+'\t'.join(l)+'\n'
            fd = open(fn,'w')
            fd.write(s_header)
            fd.close()

        fd = open(fn,'a')
        fd.write(s)
        fd.close()

        self.execmd('sort %s -k1,1 -u -o %s' %(fn,fn))

        return


    def table_sexcheck(self):

        fn_table = fn_out = 'sexcheck.table'

        if os.path.isfile(fn_table):
            return

        cmd = '''awk '{if($5=="PROBLEM") print $3,$4,$6}' %s.sexcheck''' %(
            self.o,
            )
        l_sexcheck = [s.split() for s in os.popen(cmd).readlines()]
        if self.verbose == True:
            print()
            print('sexcheck')
            print('PEDSEX SNPSEX F')
            print(l_sexcheck)

        cmd = 'cat %s.sexcheck' %(self.o,)
        cmd += ''' | awk '{if($5=="PROBLEM") print $2,$3,$4,$6}' '''
        lines = os.popen(cmd).readlines()
        lines = ['%s\t%s\n' %(self.i,'\t'.join(line.split())) for line in lines]

        if not os.path.isfile(fn_out):
            fd = open(fn_out,'w')
            fd.write('# IID PEDSEX SNPSEX F\n')
            fd.close()

        fd = open(fn_out,'a')
        fd.writelines(lines)
        fd.close()

        cmd = 'sort -u %s | sort -k1,1 -o %s' %(fn_out,fn_out,)
        self.execmd(cmd)

        return l_sexcheck


    def table_hwe(self):

        l_range = list(range(4,8+1))

        for suffix in ['','.X.females',]:

            fn_table = 'hwe%s.table' %(suffix)
            fn_hwe = '%s%s.hwe' %(self.o,suffix,)

            if os.path.isfile(fn_table): continue
            if not os.path.isfile(fn_hwe): continue

            ## do loop over thresholds even though it's not optimal...
            s = '%s' %(self.i)
            for i in l_range:
                f = 10**-i
                cmd = 'grep ALL %s' %(fn_hwe)
                cmd += " | awk '{if ($9 < %.1e) print $2}' " %(f,)
                cmd += ' | wc -l'
                n = int(os.popen(cmd).read())
                print(i, n)
                s += '\t%i' %(n)
            s += '\n'

            if not os.path.isfile(fn_table):
                l = ['']
                for i in l_range:
                    l += ['10^-%i' %(i)]
                s_header = '\t'.join(l)+'\n'
                fd = open(fn_table,'w')
                fd.write(s_header)
                fd.close()

            fd = open(fn_table,'a')
            fd.write(s)
            fd.close()

            self.execmd('sort %s -k1,1 -u -o %s' %(fn_table,fn_table,))

        return


    def table_lmiss(self):

        for suffix in ['SNPQC','X.males','X.females',]:

            if suffix == 'X.males' and self.bool_filter_females == True: continue

            if not os.path.isfile('%s.%s.lmiss' %(self.o,suffix,)):
                continue

            fn_table = 'lmiss.%s.table' %(suffix)

            if os.path.isfile(fn_table):
                continue

            s = ''
            cmd = '''awk 'NR>1{print 1-$5}' %s.%s.lmiss''' %(
                self.o,suffix,
                )
            l_lmiss = [float(val) for val in os.popen(cmd).readlines()]
    ##        s += '%s\n' %(str(stats.relfreq(l_lmiss,numbins=6,defaultreallimits=(.97,1,))))
    ##        s += '%s\n' %(str(stats.relfreq(l_lmiss,numbins=10,defaultreallimits=(.98,1,))))
    ##        s += '%s\n' %(str(stats.relfreq(l_lmiss,numbins=10,defaultreallimits=(.99,1,))))
            lowerreallimit = .950
            binsize = 0.005
            numbins = int((1-lowerreallimit)/binsize)
            l_cumfreq = stats.cumfreq(
                a=l_lmiss, numbins=numbins, defaultreallimits=(lowerreallimit,1,),
                )

            if self.bool_verbose == True:
                print('cumfreq+extrapoints', l_cumfreq[0])
                print('lowerreallimit', l_cumfreq[1])
                print('binsize', l_cumfreq[2])
                print('extrapoints', l_cumfreq[3])

            self.table_header_cumfreq(fn_table,lowerreallimit,binsize,)

            s = self.cumfreq2string(l_cumfreq,)
            line = '%s\t%s\n' %(self.i,s,)
            fd = open(fn_table,'a')
            fd.write(line)
            fd.close()

            self.execmd('sort %s -k1,1 -u -o %s' %(fn_table,fn_table,))

        return


    def table_header_cumfreq(self,fn_table,lowerreallimit,binsize,):

        if not os.path.isfile(fn_table):
            l = ['#']
            for i in range(int((1-lowerreallimit)/binsize)+1):
                l += ['%.3f' %(lowerreallimit+i*binsize)]
            s = '\t'.join(l)+'\n'
            fd = open(fn_table,'w')
            fd.write(s)
            fd.close()

        return
    

    def table_genome(self):

        if (
            os.path.isfile('genome.table')
            and
            os.path.isfile('genome_outliers.table')
            ):
            return

        ## make sure other processes don't write this table
        fd = open('genome.table','w')
        fd.close()

        ##
        ## define pi hat range
        ##
        numbins = 4
        l_pi_hat_max = []
        for i in range(numbins):
            pi_hat_max = (i+1)*((.20-0.)/numbins)
            l_pi_hat_max += [pi_hat_max]
        l_pi_hat_max += [0.5,0.9,1.0,]

        ##
        ## related samples at different thresholds
        ##

        n_samples = int(os.popen('cat %s.fam | wc -l' %(self.i)).read())
        n_samples -= int(os.popen('cat %s.sampleQC.samples | wc -l' %(self.o)).read())

        l = [self.i]
        l_header = ['#']
        for pi_hat_max in l_pi_hat_max:
            if pi_hat_max == 1:
                n = n_samples
            else:
                cmd = 'python %s/QC_IBD_prune.py ' %(os.path.dirname(sys.argv[0]))
                cmd += '--pi_hat_max %.2f --genome %s.prehardy --imiss %s --out %s' %(
                    pi_hat_max, self.o, self.o, self.o,
                    )
                cmd += '\n\n'
                self.execmd(cmd)
                cmd = 'cat %s.genome.%.2f.samples | wc -l' %(self.o,pi_hat_max,)
                n = n_samples-int(os.popen(cmd).read())
            l += ['%s' %(n)]
            l_header += ['%.2f' %(pi_hat_max)]

        pi_hat_max = 1
        n = n_samples
        l += ['%s' %(n)]
        l_header += ['%.2f' %(pi_hat_max)]

        s = '\t'.join(l)+'\n'
        s_header = '\t'.join(l_header)+'\n'

        fn_table = 'genome.table'

        if not os.path.isfile(fn_table):
            fd = open(fn_table,'w')
            fd.write(s_header)
            fd.close()

        fd = open(fn_table,'a')
        fd.write(s)
        fd.close()

        self.execmd('sort %s -k1,1 -u -o %s' %(fn_table,fn_table,))

        ##
        ## highly related samples
        ##
        cmd = '''awk 'NR>1{if($10>0.9) print $1,$3,$10}' %s.prehardy.genome''' %(self.o)
        lines = os.popen(cmd).readlines()
        s = '## someheader\n'
        for line in lines:
            s += '%s\t%s\n' %(self.i,'\t'.join(line.split()))
        fn_out = 'genome_outliers.table'
        fd = open(fn_out,'a')
        fd.write(s)
        fd.close()
        self.execmd('sort -u %s | sort -k1,1 -o %s' %(fn_out,fn_out,))

        ##
        ## highly related samples
        ##
        fn = 'genome.high.table'
        cmd = "awk 'NR>1{if($10>0.90) print $1,$3,$10}'"
        cmd += ' %s.prehardy.genome > %s.genome.high' %(self.o,self.o,)
        self.execmd(cmd)
        ## sort
        cmd = 'sort -k2,2 %s.imiss > %s.imiss.sorted' %(self.o,self.o,)
        self.execmd(cmd)
        ## sort and join 1
        cmd = 'sort -k1,1 %s.genome.high > %s.genome.high.sorted' %(self.o,self.o)
        self.execmd(cmd)
        cmd = 'join -1 2 -2 1 -o 0,2.3,1.6 '
        cmd += '%s.imiss.sorted %s.genome.high.sorted >> %s' %(
            self.o,self.o,fn,)
        self.execmd(cmd)
        ## sort and join 2
        cmd = 'sort -k2,2 %s.genome.high > %s.genome.high.sorted' %(self.o,self.o)
        self.execmd(cmd)
        cmd = 'join -1 2 -2 2 -o 0,2.3,1.6 '
        cmd += '%s.imiss.sorted %s.genome.high.sorted >> %s' %(
            self.o,self.o,fn,)
        self.execmd(cmd)
        ## get rid of duplicates
        cmd = 'sort -u %s -o %s' %(fn,fn)
        self.execmd(cmd)
        ## clean up
        os.remove('%s.genome.high' %(self.o))
        os.remove('%s.genome.high.sorted' %(self.o))
        os.remove('%s.imiss.sorted' %(self.o))

        return


    def table_imiss(self):

        fn_table = 'imiss.table'
        if os.path.isfile(fn_table): return

        cmd = '''awk 'NR>1{print 1-$6}' %s.imiss''' %(self.o)
        l_imiss = [float(imiss) for imiss in os.popen(cmd).readlines()]
##        ## how many in each bin?
##        print stats.relfreq(l_miss,numbins=5,defaultreallimits=(.975,1,))
        lowerreallimit = .970
        binsize = 0.005
        numbins = int((1-lowerreallimit)/binsize)
        l_cumfreq = stats.cumfreq(
            a=l_imiss, numbins=numbins, defaultreallimits=(lowerreallimit,1,),
            )

        if self.bool_verbose == True:
            print('cumfreq', l_cumfreq[0])
            print('lowerreallimit', l_cumfreq[1])
            print('binsize', l_cumfreq[2])
            print('extrapoints', l_cumfreq[3])

        self.table_header_cumfreq(fn_table,lowerreallimit,binsize,)

        s = self.cumfreq2string(l_cumfreq,)
        line = '%s\t%s\n' %(self.i,s,)
        fd = open(fn_table,'a')
        fd.write(line)
        fd.close()

        self.execmd('sort %s -k1,1 -u -o %s' %(fn_table,fn_table))

        return l_imiss


    def cumfreq2string(self,l_cumfreq,):

        s = '\t'.join(
            [str(l_cumfreq[3])]+[str(int(v)) for v in l_cumfreq[0]+l_cumfreq[3]]
            )

        return s


    def table_het(self):

        fn_table = 'het.table'

        if os.path.isfile(fn_table):
            return

        l_cmds,average,stddev,het_min,het_max = self.het2stddev(
            float(self.args.threshold_imiss),
            bool_execute=True,bool_remove=False,)

        ## parse heterozygosities for stats.histogram
        cmd = 'cat %s.het.joined ' %(self.o,)
##        cmd += " | awk '{het=($5-$3)/$5; print het}'"
        cmd += " | awk '{print $5}'"
        l_het = [float(line) for line in os.popen(cmd).readlines()]
        os.remove('%s.het.joined' %(self.o))

        ## write stats to file
        fd = open('het.stats','w')
        fd.write('%s %f %f\n' %(
            self.i,average,self.args.threshold_het_stddev*stddev,))
        fd.close()

        ## histogram
        l_histogram = stats.histogram(
            l_het, numbins=2*int(self.args.threshold_het_stddev),
            defaultlimits=(
                average-int(self.args.threshold_het_stddev)*stddev,
                average+int(self.args.threshold_het_stddev)*stddev,
                ),
            )

        ## cumulated frequency...
        l_cumfreq = stats.cumfreq(
            l_het,numbins=1,defaultreallimits=(
                0,
                average+int(self.args.threshold_het_stddev)*stddev,
                ),
            )

        histogram_extrapoints = l_histogram[3]
        cumfreq_extrapoints = l_cumfreq[3]

        if self.bool_verbose == True:
            print('histogram', l_histogram[0])
            print('less than lower', int(
##                sum(l_cumfreq[0])-sum(l_histogram[0][:])
                histogram_extrapoints-cumfreq_extrapoints
                ))
            print('greater than upper', int(
##                len(l_het)-sum(l_cumfreq[0])
                cumfreq_extrapoints
                ))

        ## less than lower
        l = [str(int(histogram_extrapoints-cumfreq_extrapoints))]
        l += [str(int(v)) for v in l_histogram[0]]
        ## greater than upper
        l += [str((cumfreq_extrapoints))]
        l += [str(len(l_het))]

        s = '\t'.join(l)
        tr = '%s\t%s\n' %(self.i,s,)

        l = ['-$\\infty$']
        l += [str(i) for i in range(
            -int(self.args.threshold_het_stddev),
            +int(self.args.threshold_het_stddev)+1,
            )]
        l += ['$\\infty$']
        l = ['[%s:%s)' %(l[i],l[i+1]) for i in range(len(l)-1)]
        th = '#\t'+'\t'.join(l)+'\n'

        fd = open(fn_table,'a')
        fd.write(th+tr)
        fd.close()

        self.execmd('sort %s -k1,1 -u -o %s' %(fn_table,fn_table))

        return l_het


    def xfrange(self, start, stop, step):
        while start < stop:
            yield start
            start += step


    def plink_figures(self):

        for l_suffixes_in,function in [
            [['imiss'],self.histogram_imiss,],
            [['SNPQC.lmiss'],self.histogram_lmiss,],
            [['het'],self.histogram_het,],
            [['postIBD.frq'],self.histogram_frq,],
            [['prehardy.genome'],self.histogram_genome,],
            [['imiss','het','sexcheck',],self.scatter_het_call,],
            [['postIBD.frq','hwe',],self.scatter_frq_hwe,],
            [['hwe',],self.histogram_hwe,],
            [['hwe',],self.scatter_hwe,],
##            [['SNPQC.lmiss','frq',],self.scatter_lmiss_frq,],
##            [['mds',],self.scatter_mds,],
##            [['prehardy.mds',],self.scatter_mds_excl_1000g,],
            [['posthardy.mds',],self.scatter_mds_excl_1000g,],
            [['%s.mds' %(self.fn1000g),],self.scatter_mds_incl_1000g,],
            [['posthardy.EIGENSOFT.evec',],self.scatter_PCA,],
            [['%s.EIGENSOFT.evec' %(self.fn1000g),],self.scatter_PCA,],
            ]:
            bool_continue = False
            for in_suffix in l_suffixes_in:
                fp_in = '%s.%s' %(self.o,in_suffix,)
                if not os.path.isfile(fp_in):
                    bool_continue = True
                    break
                if not os.path.getsize(fp_in) > 0:
                    bool_continue = True
                    break
                bool_input = self.check_file_in(fp_in)
                if bool_input == False:
                    bool_continue = True
                    break
            if bool_continue == True:
                continue
            ## lock
            if os.path.isfile('%s.%s.lock' %(self.o,in_suffix,)):
                print('other process plotting. exiting.')
                sys.exit(0)
            else:
                self.execmd('touch %s.%s.lock' %(self.o,in_suffix,))
            ## plot
            if l_suffixes_in[0][-4:] == 'evec':
                function(l_suffixes_in[0][:-len('.EIGENSOFT.evec')],)
            else:
                function()
            ## unlock
            os.remove('%s.%s.lock' %(self.o,in_suffix,))

##        if (
##            not os.path.isfile('venn3_%s.png' %(self.o))
##            and
##            os.path.isfile('%s.%s.samples' %(self.o,'imiss',))
##            and
##            os.path.isfile('%s.%s.samples' %(self.o,'het',))
##            and
##            os.path.isfile('%s.%s.samples' %(self.o,'sexcheck',))
##            ):
##            cmd = 'cat %s.het.samples' %(self.o)
##            l_het = os.popen(cmd).read().split()
##            cmd = 'cat %s.imiss.samples' %(self.o)
##            l_imiss = os.popen(cmd).read().split()
##            cmd = 'cat %s.sexcheck.samples' %(self.o)
##            l_sexcheck = os.popen(cmd).read().split()
####            print len(l_het)
####            print len(l_imiss)
####            print len(l_sexcheck)
####            print 'imiss only', len(set(l_imiss)-set(l_sexcheck+l_het))
####            print 'het only', len(set(l_het)-set(l_imiss+l_sexcheck))
####            print 'sexcheck only', len(set(l_sexcheck)-set(l_imiss+l_het))
####            print 'imiss and het', len(set(set(l_imiss)&set(l_het))-set(l_sexcheck))
####            print 'het and sexcheck', len(set(set(l_sexcheck)&set(l_het))-set(l_imiss))
####            print 'imiss and sexcheck', 0, len(set(set(l_imiss)&set(l_sexcheck))-set(l_het))
####            print 'combined', len(set(l_sexcheck)|set(l_imiss)|set(l_het))
##            if l_het != [] and l_imiss != [] and l_sexcheck != []:
##                print('venn3', len(l_het), len(l_imiss), len(l_sexcheck))

        return


    def histogram_hwe(self):

        if os.path.isfile('%s.hwe.png' %(self.o)):
            return
        
        cmd = 'cat %s.hwe' %(self.o)
        cmd += ' | '
        cmd += "awk 'NR>1{"
        cmd += 'if ($9!="NA") {logp=-log($9)/log(10); if(logp<=10) print logp}'
        cmd += "}'"
        cmd += '> %s.hwe.dat' %(self.o)
        self.execmd(cmd)

        n_samples = int(os.popen('cat %s.fam | wc -l' %(self.i)).read())
        cmd = 'cat %s.sampleQC.nonfounders.samples | wc -l' %(self.o)
        n_samples -= int(os.popen(cmd).read())

        n_SNPs = (int(os.popen('cat %s.hwe | wc -l' %(self.o)).read())-1)/3

        if self.bool_verbose == True:
            bool_timestamp = True

        self.gnuplot_histogram(
##                '%s.hwe' %(bfile),
            '%s.hwe.dat' %(self.o),
            prefix_out = '%s.hwe' %(self.o),
##            x_step=0.01,
            x_step=0.1,
##            x_max=0.1,
##            column='$9',
            column='$1',
            xlabel='-log_{10}({/Helvetica-Italic p}_{HWE})',
            title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                self.o.replace('_','\\\\_'),n_samples,n_SNPs,
                ),
            color = 'yellow',
            bool_timestamp = True,
            )

        os.remove('%s.hwe.dat' %(self.o))

        return


    def wc(self,fn):

        with open(fn) as fd:
            n = len(fd.readlines())

        return n


    def fam2samples(self,bfile):

        fam = '%s.fam' %(bfile)
        with open(fam) as fd:
            l_samples = [line.split()[1] for line in fd.readlines()]

        return l_samples


    def scatter_PCA(self,affix,):

        prefix = '%s.%s' %(self.o,affix,)

        fn_evec = '%s.EIGENSOFT.evec' %(prefix)
        if os.path.isfile('%s.pc1.pc2.png' %(fn_evec)):
            return

        ##
        ## count samples and SNPs
        ##
        n_samples = self.wc(fn_evec)-1
        n_SNPs = self.wc('%s.prune.in' %(prefix))

        ## is the pop dic file already sorted by population?
        cmd = 'sort -k2,2 %s -o %s.sorted' %(
            self.args.fn_pops,os.path.basename(self.args.fn_pops))
        self.execmd(cmd)

        ## assign number to each population
        ## convert PLINK sample ID to EIGENSOFT sample ID
        cmd = 'cat %s.sorted' %(os.path.basename(self.args.fn_pops))
        cmd += " | awk '"
        cmd += ' BEGIN{i=0;pop=""}'
        cmd += ' {'
        cmd += ' ID1=$1;sub(/-RECLUSTER/,"",ID1);'
        cmd += ' sub(/-RESCAN/,"",ID1);'
        cmd += ' ID2=substr(ID1,length(ID1)-18,19);'
        cmd += ' if($2!=pop) {pop=$2; i++};'
        cmd += ' print ID2":"ID2,i,$2'
        cmd += ' }'
        cmd += " '"
        ## sort by sample ID
        cmd += ' | sort -k1,1 -u'
        cmd += ' > %s.EIGENSOFT.sorted' %(os.path.basename(self.args.fn_pops))
        self.execmd(cmd)
        os.remove('%s.sorted' %(os.path.basename(self.args.fn_pops)))

        ##
        ## sort eigenvectors by sample ID
        ##
        cmd = 'sort -k1,1 %s > %s.sorted' %(fn_evec,fn_evec)
        self.execmd(cmd)

        ##
        ## join (pop str and pop integer with eigenvectors via sample ID)
        ##
        cmd = 'join -1 1 -2 1 %s.EIGENSOFT.sorted %s.sorted' %(
            os.path.basename(self.args.fn_pops),fn_evec)
        ## sort by population in previous sorted order
        cmd += ' | sort -k2n,2'
        cmd += ' > %s.joined' %(fn_evec)
        self.execmd(cmd)
        ## remove sorted files after join
        os.remove('%s.EIGENSOFT.sorted' %(os.path.basename(self.args.fn_pops)))
        os.remove('%s.sorted' %(fn_evec))

        ##
        ## check that all samples/pops accounted for
        ##
        if self.wc('%s.joined' %(fn_evec))+1 != self.wc('%s' %(fn_evec)):
            print('evec joined with poplist',self.wc('%s.joined' %(fn_evec)),'evec prejoin',self.wc('%s' %(fn_evec)))
            print('%s.joined' %(fn_evec),fn_evec)
            stop_add_to_pop_dic

        ##
        ## parse populations
        ##
        cmd = "cat %s.joined | awk '{print $3}' | uniq" %(fn_evec)
        l_pops = os.popen(cmd).read().strip().split()

        ##
        ## specify colors
        ##
        l_colors1 = [[15*i,15*i,15*i] for i in range(17)]
        l_colors2 = [
            [255,0,0,], ## red
            [255,51,0,],
            [255,102,0,],
            [255,153,0,],
            [255,204,0,],
            [255,255,0,],
            [204,255,0,],
            [153,255,0,],
            [102,255,0,],
            [51,255,0,],
            [0,255,0,], ## green
            [0,255,51,],
            [0,255,102,],
            [0,255,153,],
            [0,255,204,],
            [0,255,255,],
            [0,204,255,],
            [0,153,255,],
            [0,102,255,],
            [0,51,255,],
            [0,0,255,],
            [51,0,255,],
            [102,0,255,],
            [153,0,255,],
            [204,0,255,],
            [255,0,255,],
            [255,0,204,],
            [255,0,153,],
            [255,0,102,],
            [255,0,51,],
            ]
        l_colors = l_colors2+l_colors1

        ## point types
##        l_pt = [5,7,9,11,]
        l_pt = [17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]

        ## what pops part of this project?
        l_samples_project = self.fam2samples(self.i)
        with open(self.args.fn_pops) as fd:
            d_samples2pops = {line.split()[0]:line.split()[1] for line in fd.readlines()}
        set_pops = set()
        for sample in l_samples_project:
            try:
                set_pops |= set([d_samples2pops[sample]])
            except KeyError:
                continue

        for pc1,pc2 in [[1,2],[3,4]]:
            line_plot = ''
##            ## http://www.gnuplotting.org/using-a-palette-as-line-color/
##            line_plot += 'h1 = 0/360.\nh2 = 320/360.\nset palette model HSV functions (1-gray)*(h2-h1)+h1,1,1\n'
##    ##        ## http://gnuplot.sourceforge.net/demo/pm3dcolors.html
##    ##        line_plot += 'set palette model HSV defined ( 0 0 1 1, max 1 1 1 )\n'
##    ####        line_plot += 'set palette defined ('
##    ####        line_plot += ' 1 "#ff000", 2 "#ff5500", 3 "#ffaa00", 4 "#ffff00",'
##    ####        line_plot += ' 5 "#aaff0", 6 "#55ff00", 7 "#00ff00", 8 "#00ff55",'
##    ####        line_plot += ' 9 "#00ffaa", 10 "#00ffff", 11 "#00aaff", 12 "#0055ff",'
##    ####        line_plot += ' 13 "#0000ff", 14 "#5500ff", 15 "#aa00ff", 16 "#ff00ff",'
##    ####        line_plot += ' 17 "#ff00aa", 18 "#ff0055", 19 "#ff5555", 20 "#ffaaaa",'
##    ####        line_plot += ' 21 "#55ff55", 22 "#aaffaa", 23 "#5555ff", 24 "#aaaaff",'
##    ####        line_plot += ' 25 "#000000", 26 "#111111", 27 "#222222", 28 "#333333",'
##    ####        line_plot += ' 29 "#444444", 30 "#555555", 31 "#666666", 32 "#777777",'
##    ####        line_plot += ' 33 "#888888", 34 "#999999", 35 "#aaaaaa", 36 "#bbbbbb",'
##    ####        line_plot += ' 37 "#cccccc", 38 "#dddddd", 39 "#eeeeee",'
##    ####        line_plot = line_plot[:-1]
##    ####        line_plot += ' )\n'
##            line_plot += 'unset colorbox\n'
            line_plot += 'plot'
            i_style = 0
            for i_pop in range(len(l_pops)):
                pop = l_pops[i_pop]
                cmd = '''awk '{if($3=="%s")''' %(pop)
                cmd += " print $0}'"
                cmd += ' %s.joined' %(fn_evec)
                print(cmd)
                if int(os.popen('%s | wc -l' %(cmd)).read()) == 0: continue
                color = l_colors[i_style]
                if pop == 'Luo': color = [255,0,0]
                elif pop == 'LWK': color = [0,0,255]
                elif pop == 'Banyarwanda': color = [0,255,0]
                elif pop == 'Baganda': color = [255,255,0]
                elif pop == 'Kalenjin': color = [0,255,255]
                elif pop == 'Munyankole': color = [255,0,255]
                elif pop == 'Other': color = [127,127,127]
                else: color = [0,0,0]
                i_style += 1
                line_plot += ''' "< awk '{if($3==\\"%s\\")''' %(pop)
                line_plot += " print $0}'"
                line_plot += ' %s.joined"' %(fn_evec)
                ## http://gnuplot.sourceforge.net/demo_4.2/histograms.html
    ##            line_plot += ' u %i:%i:%i:key(%i)' %(pc1+4,pc2+4,4,3,)
                line_plot += ' u %i:%i' %(pc1+3,pc2+3,)
##                line_plot += ' u %i:%i:%i' %(pc1+3,pc2+3,2,)
                if (pop == self.project or pop in set_pops) and set_pops != set(l_pops):
                    line_plot += ' ps 3'
                else:
                    line_plot += ' ps 2'
                line_plot += ' pt %i' %(l_pt[i_style%len(l_pt)])
##                ## http://gnuplot.sourceforge.net/demo_cvs/varcolor.html
##                line_plot += ' lc pal z'
                line_plot += ' lc rgb %i' %(65536*color[0]+256*color[1]+color[2])
                line_plot += ' t "%s", ' %(pop)
                continue
            ##
            ## add labels
            ##
            line_plot += ' "< cat %s.joined' %(fn_evec)
            line_plot += ''' | awk -F ':' '{print $1,$0}' | awk '{if($6>0.05||$6<-0.05)'''
            line_plot += ''' print $4$0}'"'''
##            line_plot += ' u %i:%i:($%i+0.02) w labels' %(1+pc1+3,1+pc2+3,1,)
            line_plot += ' u %i:%i:($%i+0.002) w labels' %(1+pc1+3,1+pc2+3,1,)
            line_plot += ' font "Helvetica,28" noenhanced'
            line_plot += ' notitle'
            line_plot += ' lc rgb "#000000", '
            ## EOL
            line_plot = line_plot[:-2]+'\n'

            if self.bool_verbose == True: bool_timestamp = True
            else: bool_timestamp = False

            self.gnuplot_scatter(
                '%s.evec' %(self.o),
                line_plot = line_plot,
                xlabel = 'PC%i' %(pc1),
                ylabel = 'PC%i' %(pc2),
                title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                    self.o.replace('_','\\\\_'),n_samples,n_SNPs,
                    ),
                prefix_out='%s.pc%i.pc%i' %(fn_evec,pc1,pc2,),
                lines_extra=['set key out\n'],
                bool_timestamp=bool_timestamp,
                )

        ## clean up
        os.remove('%s.joined' %(fn_evec))

        return


    def scatter_mds_excl_1000g(self):

        mds_prefix = '%s.posthardy' %(self.o)
        fn = '%s.mds' %(mds_prefix)
        if os.path.isfile('%s.1.2.mds.png' %(fn)):
            return

        '''this function needs to be rewritten.
it's ugly and I will not understand it 1 year form now.'''

        l_colors = [
            [255,0,0,],
            [255,85,0,],
            [255,170,0,],
            [255,255,0,],
            [170,255,0,],
            [85,255,0,],
            [0,255,0,],
            [0,255,85,],
            [0,255,170,],
            [0,255,255,],
            [0,170,255,],
            [0,85,255,],
            [0,0,255,],
            [85,0,255,],
            [170,0,255,],
            [255,0,255,],
            [255,0,170,],
            [255,0,85,],
            ]

        ## give option to add --samplefile flag instead!!!
        if self.project == 'uganda_gwas':
            l_tribes = [
                'Muganda',
                'Munyarwanda',
                'Rwandese Ugandan',
                'Murundi',
                'Munyankole',
                'Mukiga',
                'Unknown',
                'Other',
                'Mutanzania',
                'Musoga',
                'Mutooro',
                'Mufumbira',
##                'Nyanjiro (Tanzania)',
                ]
        else:
            l_tribes = [self.o,]

        if not os.path.isdir('mds'):
            os.mkdir('mds')
            pass

        ##
        ## map samples to their tribe/ethnicity
        ##
        if self.project == 'uganda_gwas':
            d_IID2tribe = {}
            for tribe in l_tribes:
                fd = open('samples/%s.samples' %(tribe.replace(' ','')),'r')
                s = fd.read()
                fd.close()
                l_IIDs = s.split('\n')
                for IID in l_IIDs:
                    if IID == 'NA':
                        continue
                    d_IID2tribe[IID] = tribe
                    continue
                ## empty files
                fd = open('mds/%s.mds' %(tribe),'w')
                fd.close()
                continue
            ##
            ## write individual mds files for each tribe/ethnicity to be plotted
            ##
            fd = open(fn,'r')
            lines = fd.readlines()
            fd.close()
            for line in lines[1:]:
                FID = line.split()[0][-10:]
                if FID in list(d_IID2tribe.keys()):
                    pop_prefix = d_IID2tribe[FID]
                else:
                    pop_prefix = 'other'
                    print('no tribe assigned', FID, os.popen('grep %s samples/*' %(FID)).read())
                if FID == 'APP5339441':
                    print(tribe)
                    stop
                fd = open('mds/%s.mds' %(pop_prefix),'a')
                fd.write(line)
                fd.close()
                continue
        else:
            self.execmd('cp %s mds/%s.mds' %(fn,self.o))

        ##
        ## count samples and SNPs
        ##
        n_samples = self.wc(fn)-1
        n_SNPs = self.wc('%s.prune.in' %(mds_prefix))

        ##
        ## loop over combinations of principal components
        ##
        d_columns = {1:4,2:5,3:6,4:7,}

        ##
        ## parse .mds file for subsequent calculation of MD outliers
        ##
        cmd = "cat %s | awk 'NR>1{print $2}'" %(fn)
        l_IIDs = os.popen(cmd).read().strip().split('\n')
        cmd = "cat %s | awk 'NR>1{print $4,$5,$6,$7}'" %(fn)
        s = os.popen(cmd).read()
        ## http://docs.scipy.org/doc/numpy/reference/generated/numpy.fromstring.html
        array_components = numpy.fromstring(s,sep=' ').reshape(n_samples,4,)

        ## point types
        l_pt = [5,7,9,11,]

        for pc1 in range(1,5):
            for pc2 in range(1,5):
                if (
                    (pc1 != 1 or pc2 != 2)
                    and
                    (pc1 != 3 or pc2 != 4)
                    ):
                    continue
                if self.verbose == True:
                    print('mds', pc1, pc2)
                if pc2 <= pc1: continue
        
                line_plot = 'plot '
        ##        line_plot += '"%s.mds" u 4:5 lc 0 ps 2 pt 7 t ""\n' %(bfile)
                for i_tribe in range(len(l_tribes)):
                    tribe = l_tribes[i_tribe]
                    color = l_colors[i_tribe]
##                    l_hex = [hex(col) for col in color]
                    line_plot += '"mds/%s.mds" ' %(tribe)
                    line_plot += 'u %i:%i ' %(d_columns[pc1],d_columns[pc2],)
                    line_plot += 'lc rgb %i ps 2 pt %i t "%s",' %(
                        65536*color[0]+256*color[1]+color[2],
                        l_pt[i_tribe%len(l_pt)],
                        tribe,
                        )
##                line_plot = line_plot[:-1]+'\n'

                ##
                ## add labels for MD outliers
                ##
                line_plot = self.add_MDS_outlier_labels(
                    array_components,l_IIDs,pc1,pc2,line_plot,)

                ## finalize plot line
                line_plot = line_plot[:-1]+'\n'

                if self.bool_verbose == True:
                    bool_timestamp = True

                self.gnuplot_scatter(
                    '%s.mds' %(self.o),
        ##            s_plot = '"< paste %s.het %s.imiss" u (($5-$3)/$5):(1-$12)' %(self.o,self.o),
        ##            s_plot = s_plot,
                    line_plot = line_plot,
##                    column1 = 4, column2 = 5,
                    xlabel = 'D%i' %(pc1),
                    ylabel = 'D%i' %(pc2),
                    title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                        self.o.replace('_','\\\\_'),n_samples,n_SNPs,
                        ),
                    prefix_out='%s.%i.%i.mds' %(fn,pc1,pc2,),
                    lines_extra=['set key out\n'],
##                    bool_remove=False,
                    bool_timestamp=bool_timestamp,
                    )

        return


    def add_MDS_outlier_labels(self,array_components,l_IIDs,pc1,pc2,line_plot,):

        u = array_components[:,pc1-1]
        v = array_components[:,pc2-1]
        ## http://docs.scipy.org/doc/numpy/reference/generated/numpy.cov.html
##                cov_pooled = np.cov(u,v,ddof=-1)
        cov_pooled = numpy.cov(u,v)
        cov_pooled_inverted = numpy.linalg.inv(cov_pooled)
        lines = ['#PC1\tPC2\tIID\tMD F p\n']
        for i in range(len(u)):
            diff = numpy.array([numpy.mean(u)-u[i],numpy.mean(v)-v[i]])
            MD_squared = mahalanobis_squared = numpy.dot(
                numpy.dot(
                    numpy.transpose(diff),cov_pooled_inverted,
                    ),
                diff,
                )
            p = 2
            n = len(u)
##                    F = MD_squared*((n-p-1)/(p*(n-2)))*((n-1)/(n))
            F = MD_squared*(n-1)*(n-p-1)/(p*n*(n-2))
##                    prob = 1-instance_tests.fdist(F,p,n-p-1)
            prob = 'NA'
            IID = l_IIDs[i]
##                    print n,MD_squared,F,prob
            if MD_squared > 5**2:
                lines += ['%f\t%f\t%s\t%s %s %s\n' %(
                    u[i],v[i],IID,math.sqrt(MD_squared),F,prob,)]
            else:
                lines += ['%f\t%f\t \n' %(
                    u[i],v[i],)]
        fd = open('%s.mds.labels' %(self.o),'w')
        fd.writelines(lines)
        fd.close()
        
        line_plot += '"%s.mds.labels" u 1:2:3 w labels ' %(self.o)
        line_plot += 'font "Helvetica,20" noenhanced t ""'

        return line_plot


    def scatter_mds_incl_1000g(self):

        '''this function needs to be rewritten.
it's ugly and I will not understand it 1 year form now.'''

        if os.path.isfile('%s.%s.mds.1.2.png' %(self.o,self.fn1000g,)):
            return

        if not os.path.isdir('mds'):
            os.mkdir('mds')

        if self.project != 'uganda_gwas':
            return

        d_sample2pop = {}
        d_pops = {}
        for fn in os.listdir('samples'):
            ethnicity = fn[:-8]
            fd = open('samples/%s' %(fn),'r')
            s = fd.read().strip()
            fd.close()
            l = s.split('\n')
            d_pops[ethnicity] = []
            for s in l:
                if s == 'NA':
                    continue
                if ethnicity == '':
                    continue
                if s in list(d_sample2pop.keys()):
                    if ethnicity == 'Unknown':
                        continue
                    print(s, d_sample2pop[s], ethnicity)
                    stop
                d_sample2pop[s] = ethnicity
                d_pops[ethnicity] += [s]
            fd = open('mds/%s_%s.mds' %(self.o,ethnicity),'w')
            fd.close()
            continue

        fn_mds = fn = '%s.%s.mds' %(self.o,self.fn1000g,)

        fd = open(fn,'r')
        lines = fd.readlines()
        fd.close()
        x = 0
        for line in lines[1:]:
            x += 1
            FID = line.split()[1]
            if 'APP' in FID:
                FID = FID[-10:]
            if '_' in FID:
                FID = FID[FID.rindex('_')+1:]
            prefix = d_sample2pop[FID]
##            print prefix
            fd = open('mds/%s_%s.mds' %(self.o,prefix),'a')
            fd.write(line)
            fd.close()
##            print FID, os.popen('grep %s samples/*' %(FID)).read()
##        stop

        i = 0
        d_colors = {}
        ## this is specific to the uganda_gwas project and should be a command line option...
        ## get rid of magic numbers/lists!!!
        l_tribes = [
            'Muganda',
            'Munyarwanda',
            'RwandeseUgandan',
            'Murundi',
            'Munyankole',
            'Mukiga',
            'Unknown',
            'Other',
            'Mutanzania',
            'Musoga',
            'Mutooro',
            'Mufumbira',
##                'Nyanjiro (Tanzania)',
            ]
        for tribe in l_tribes:
            d_colors[tribe] = [i*20,i*20,i*20,]
            i += 1

        l_colors = [
            [255,0,0,],
            [255,85,0,],
            [255,170,0,],
            [255,255,0,],
            [170,255,0,],
            [85,255,0,],
            [0,255,0,],
            [0,255,85,],
            [0,255,170,],
            [0,255,255,],
            [0,170,255,],
            [0,85,255,],
            [0,0,255,],
            [85,0,255,],
            [170,0,255,],
            [255,0,255,],
            [255,0,170,],
            [255,0,85,],
            ]

        l_pts = [5,7,9,11,]

        l_pops = list(set(d_pops.keys())-set(l_tribes))
        l_pops.sort()
        l_pops += l_tribes

        for c1,c2 in [
            [1,2,],
            [3,4,],
            ]:
            line_plot = 'set key out\nplot '
    ##        line_plot = 'plot [:0.25]'
    ##        line_plot = 'plot [:0][-0.01:0.02]'
            i_color = 0
            d_plot = {'uganda':'','1000g':'',}
            for i_pop in range(len(l_pops)):
                pop = l_pops[i_pop]
                fn_mds_pop = 'mds/%s_%s.mds' %(self.o,pop)
                if not os.path.isfile(fn_mds_pop):
                    print('mds skipping', pop)
                    continue
                if os.path.getsize(fn_mds_pop) == 0:
                    print('mds skipping', pop)
                    continue
                if pop in list(d_colors.keys()):
                    color = d_colors[pop]
                    ps = 2
                    pt = 7
                else:
                    print(pop)
                    color = l_colors[i_color%len(l_colors)]
    ##                color = l_colors[i_color]
                    i_color += 1
                    ps = 2
                    pt = l_pts[i_color%len(l_pts)]
                s = '"%s" u %i:%i lc rgb "#%s" ps %i pt %s t "%s",' %(
                    fn_mds_pop,
                    4-1+c1,4-1+c2,
                    "".join(map(chr, color)).encode('hex'),
                    ps,pt,
                    pop,
                    )
                if pop in list(d_colors.keys()):
                    d_plot['uganda'] += s
                else:
                    d_plot['1000g'] += s
            line_plot += d_plot['uganda']
            line_plot += d_plot['1000g']
            line_plot = line_plot[:-1]+'\n'
            n_samples = self.wc(fn_mds)
            n_samples -= 1
            cmd = 'cat %s.%s.prune.in | wc -l' %(self.o,self.fn1000g,)
            n_SNPs = int(os.popen(cmd).read())

            if self.bool_verbose == True:
                bool_timestamp = True

            self.gnuplot_scatter(
                '%s.%s.mds' %(self.o,self.fn1000g,),
                line_plot = line_plot,
                column1 = 4-1+c1, column2 = 4-1+c2,
                xlabel = 'D%i' %(c1),
                ylabel = 'D%i' %(c2),
                title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                    self.o.replace('_','\\\\_'),n_samples,n_SNPs,
                    ),
                prefix_out = '%s.%s.mds.%i.%i' %(self.o,self.fn1000g,c1,c2,),
                bool_timestamp=bool_timestamp,
                )

        return


    def calc_samples(self,fn_subtract,):

        n_samples = int(os.popen('cat %s.fam | wc -l' %(self.i)).read())-1
        cmd = 'cat %s | wc -l' %(fn_subtract)
        n_samples -= int(os.popen(cmd).read())

        return n_samples


    def scatter_frq_hwe(self): ## perhaps do contour instead...

        if os.path.isfile('%s.frq.hwe.png' %(self.o)):
            return

        n_samples = self.calc_samples('%s.sampleQC.IBD.samples' %(self.o))
        n_SNPs = int(os.popen('cat %s.postIBD.frq | wc -l' %(self.o)).read())-1

        ## sort and join SNPs
        cmd = "sed '1d' %s.postIBD.frq | sort -k2,2 > %s.postIBD.frq.sorted" %(self.o,self.o,)
        self.execmd(cmd)
        cmd = "cat %s.hwe | awk '{if(NR%%3==2) print}' | sort -k2,2" %(self.o,)
        cmd += " | awk '{logp=-log($9)/log(10); print $2,logp}'"
        cmd += " > %s.hwe.sorted" %(self.o,)
        self.execmd(cmd)
        cmd = 'join -1 2 -2 1 -o 0,1.5,2.2'
        cmd += ' %s.postIBD.frq.sorted %s.hwe.sorted' %(self.o,self.o)
        cmd += ' > %s.frq.hwe.joined' %(self.o)
        self.execmd(cmd)

        os.remove('%s.postIBD.frq.sorted' %(self.o))
        os.remove('%s.hwe.sorted' %(self.o))

        line_plot = 'plot [:][4:8]'
        line_plot += '"%s.frq.hwe.joined" ' %(self.o)
        line_plot += 'u 2:3 lc 0 ps 2 pt 7 t ""'
        line_plot += '\n'

        if self.bool_verbose == True:
            bool_timestamp = True

        self.gnuplot_scatter(
            '%s.frq.hwe' %(self.o),
            line_plot = line_plot,
            column1 = '(1-$6)', column2 = 12,
            xlabel = 'MAF',
            ylabel = '-log(p_{HWE})',
            title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                self.o.replace('_','\\\\_'),n_samples,n_SNPs,
                ),
            bool_timestamp=bool_timestamp,
            )

        os.remove('%s.frq.hwe.joined' %(self.o))

        return


    def scatter_het_call(self,bool_with_stddev=True,):

        '''this function is a mess. i was in a hurry, when I wrote it'''

        if os.path.isfile('%s.het.call.png' %(self.o)):
            return

        ## check that number of samples are equal
        ## which they should be if checks are run in *parallel*
        if (
            os.popen('cat %s.het | wc -l' %(self.o)).read()
            !=
            os.popen('cat %s.imiss | wc -l' %(self.o)).read()
            ):
            stop
        ## sort and join lists of samples in case they are not of equal length
        cmd = "sed '1d' %s.het | sort > %s.het.sorted" %(self.o,self.o,)
        self.execmd(cmd)
        cmd = "sed '1d' %s.imiss | sort > %s.imiss.sorted" %(self.o,self.o,)
        self.execmd(cmd)
        cmd = 'join %s.imiss.sorted %s.het.sorted > %s.imiss.het.joined' %(
            self.o,self.o,self.o,
            )
        self.execmd(cmd)

        os.remove('%s.het.sorted' %(self.o))
        os.remove('%s.imiss.sorted' %(self.o))

        if (
            os.popen('cat %s.het | wc -l' %(self.o)).read()
            !=
            os.popen('cat %s.sexcheck | wc -l' %(self.o)).read()
            ):
            stop

        ##
        ## join het and imiss (sexcheck samples only)
        ##
        cmd = "cat %s.sexcheck.samples | awk '{print $1}' " %(self.o)
        cmd += '> %s.sexcheck.samples.1column' %(self.o)
        self.execmd(cmd)
        ## use comm and join instead of grep... safer if IIDs are not unique...
        ## join1 (imiss)
        cmd = 'grep -f %s.sexcheck.samples.1column %s.imiss | ' %(self.o,self.o,)
        cmd += "awk '{print $1, 1-$6}' > %s_join1.txt" %(self.o,)
        self.execmd(cmd)
        ## join2 (het)
        cmd = 'grep -f %s.sexcheck.samples.1column %s.het' %(self.o,self.o,)
        cmd += ' | '
        cmd += "awk '{print $1,$5}' > %s_join2.txt" %(self.o,)
        self.execmd(cmd)
        ## join het and imiss
        cmd = "join %s_join1.txt %s_join2.txt > %s.sex.dat" %(
            self.o,self.o,self.o,
            )
        self.execmd(cmd)

        ## clean up
        os.remove('%s_join1.txt' %(self.o))
        os.remove('%s_join2.txt' %(self.o))
        os.remove('%s.sexcheck.samples.1column' %(self.o))

        ## opposite sex
        cmd = 'cat %s.sexcheck' %(self.o,)
        cmd += ''' | awk '{if($5=="PROBLEM" && $4!="0") print $1}' '''
        cmd += '> %s.opposite.txt' %(self.o)
        self.execmd(cmd)

        ## unknown sex
        cmd = 'cat %s.sexcheck' %(self.o,)
        cmd += ''' | awk '{if($5=="PROBLEM" && $4=="0") print $1}' '''
        cmd += '> %s.unknown.txt' %(self.o)
        self.execmd(cmd)

        ## use join instead of grep... safer in case of non-unique IDs...
        for suffix in ['opposite','unknown',]:
            ## join1 (imiss)
            cmd = "grep -f %s.%s.txt %s.imiss" %(self.o,suffix,self.o,)
            cmd += "| awk '{print $1, 1-$6}'"
            cmd += "> %s_join1.txt" %(self.o,)
            self.execmd(cmd)
            ## join2 (het)
            cmd = "grep -f %s.%s.txt %s.het" %(self.o,suffix,self.o,)
            cmd += "| awk '{print $1,$5}'"
            cmd += "> %s_join2.txt" %(self.o,)
            self.execmd(cmd)
            ## rm
            os.remove('%s.%s.txt' %(self.o,suffix,))
            ## join imiss het
            cmd = 'join %s_join1.txt %s_join2.txt > %s.sex.%s.dat' %(
                self.o,self.o,self.o,suffix,
                )
            self.execmd(cmd)
            os.remove('%s_join1.txt' %(self.o))
            os.remove('%s_join2.txt' %(self.o))

        ##
        ##
        ##
        n_samples = int(os.popen('cat %s.het | wc -l' %(self.o)).read())-1
        if bool_with_stddev == True:

            l_cmds,average,stddev,het_min,het_max = self.het2stddev(
                float(self.args.threshold_imiss),bool_execute=True)

        cmd = 'cat %s.autosomes.SNPs | wc -l' %(self.o)
        n_SNPs = int(os.popen(cmd).read())

        ## horizontal lines
        l_arrows = [
##            'set arrow from graph(0),%f to graph(1),%f nohead lc 0\n' %(
            'set arrow from 0,%f to 1,%f nohead lc 0 lt 2 lw 2\n' %(
                float(self.args.threshold_imiss),
                float(self.args.threshold_imiss),
                ),
            'set arrow from 0,%f to 1,%f nohead lc 0 lt 2 lw 2\n' %(
                0.99,0.99,
                ),
            'set arrow from 0,%f to 1,%f nohead lc 0 lt 2 lw 2\n' %(
                0.98,0.98,
                ),
            'set arrow from 0,%f to 1,%f nohead lc 0 lt 2 lw 2\n' %(
                0.97,0.97,
                ),
            ]
        ## vertical lines
        if bool_with_stddev == True:
            l_arrows += ['set arrow from %f,graph(0) to %f,graph(1) ' %(
                average-self.args.threshold_het_stddev*stddev,
                average-self.args.threshold_het_stddev*stddev,
                )]
            l_arrows += ['nohead lc 0 lt 1 lw 3\n']
            l_arrows += ['set arrow from %f,graph(0) to %f,graph(1) ' %(
                average+self.args.threshold_het_stddev*stddev,
                average+self.args.threshold_het_stddev*stddev,
                )]
            l_arrows += ['nohead lc 0 lt 1 lw 3\n']

        if bool_with_stddev == True:
            d_colors = {
                .97:'green',.99:'blue',
                .00:'black',
                #.98:'orange',
                }
            l_thresholds = list(d_colors.keys())
            l_thresholds.sort()
            l_thresholds.reverse()
            for threshold in l_thresholds:

                try:
                    (
                        l_cmds,average_w_exclusion,stddev_w_exclusion,het_min,het_max,
                        ) = self.het2stddev(
                            threshold,bool_execute=True,)
                ## perhaps nothing above threshold
                except IndexError:
                    continue
                
                color = d_colors[threshold]
                l_arrows += [
                    'set arrow from %f,graph(0) to %f,graph(1) ' %(
                        average_w_exclusion-self.args.threshold_het_stddev*stddev_w_exclusion,
                        average_w_exclusion-self.args.threshold_het_stddev*stddev_w_exclusion,
                        ),
                    'nohead lc rgb"%s" lt 2 lw 3\n' %(color),
                    'set arrow from %f,graph(0) to %f,graph(1) ' %(
                        average_w_exclusion+self.args.threshold_het_stddev*stddev_w_exclusion,
                        average_w_exclusion+self.args.threshold_het_stddev*stddev_w_exclusion,
                        ),
                    'nohead lc rgb"%s" lt 2 lw 3\n' %(color),
                    ]
    ##            print threshold, len(l_het), average_w_exclusion, 3*stddev_w_exclusion

        if bool_with_stddev == True:
            line_plot = 'plot [%f:%f][:1]' %(
                min(het_min,average-self.args.threshold_het_stddev*stddev,)-0.001,
                max(het_max,average+self.args.threshold_het_stddev*stddev,)+0.001,
                )
        else:
            line_plot = 'plot [:][:1]'
        line_plot += '"%s.imiss.het.joined" ' %(self.o)
        line_plot += 'u ($10):(1-$6) lc 0 ps 2 pt 7 t ""'

        ## sexcheck dots and labels
        l_colors = ['red','orange',]
        for i in range(2):
            suffix = ['opposite','unknown',][i]
            color = l_colors[i]
            if not os.path.getsize('%s.sex.%s.dat' %(self.o,suffix,)) > 0:
                continue
            line_plot += ','
            line_plot += '"%s.sex.%s.dat" u 3:2 t "" lc rgb"%s" ps 2 pt 7' %(
                self.o,suffix,color,
                ) ## het and sex points
            line_plot += ','
            line_plot += '"%s.sex.%s.dat" u 3:($2+0.002):1 t "" w labels ' %(
                self.o,suffix,
                ) ## het and sex labels
            line_plot += 'font "Helvetica,20" noenhanced lc 2'
            continue
        line_plot += '\n'

        if self.bool_verbose == True:
            bool_timestamp = True

        self.gnuplot_scatter(
            '%s.het.call' %(self.o),
##            s_plot = '"< paste %s.het %s.imiss" u (($5-$3)/$5):(1-$12)' %(bfile,bfile,),
            line_plot = line_plot,
            column1 = '(1-$6)', column2 = 12,
            xlabel = 'heterozygosity',
            ylabel = 'call rate',
            title = '%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                self.o.replace('_','\\\\_'),n_samples,n_SNPs,
                ),
            lines_extra = l_arrows,
            bool_timestamp = bool_timestamp,
            )

        os.remove('%s.sex.opposite.dat' %(self.o))
        os.remove('%s.sex.unknown.dat' %(self.o))
        os.remove('%s.sex.dat' %(self.o))
        os.remove('%s.imiss.het.joined' %(self.o))

        return


    def het2stddev(
        self,threshold_imiss,bool_execute=False,bool_remove=True,):

        l_cmds = []

        cmd = 'cat %s.imiss' %(self.o)
        cmd += " | awk 'NR>1 {if($6>%f) print $1,$2}'" %(1-threshold_imiss)
        cmd += ' | sort -k2,2'
        cmd += ' > %s.imiss.remove.sorted' %(self.o)
        l_cmds += [cmd]

        cmd = "cat %s.het | awk 'NR>1' | sort -k2,2 > %s.het.sorted" %(self.o,self.o,)
        l_cmds += [cmd]

        cmd = 'join -1 2 -2 2 -v2 -o 2.1,2.2,2.3,2.4,2.5,2.6'
        cmd += ' %s.imiss.remove.sorted %s.het.sorted ' %(self.o,self.o,)
        cmd += ' > %s.het.joined' %(self.o)
        l_cmds += [cmd]

        if bool_execute == True:
            for cmd in l_cmds:
                self.execmd(cmd)

        cmd = ''
        if bool_execute == False:
            ## use var=$(cmd) if bash...
            ## use set var=`cmd` if tcsh...
            ## set variable s equal to stdout of command in parenthesis
            cmd += 's=$( '
        cmd += ' cat %s.het.joined' %(self.o)
        cmd += " | awk '"
        ## BEGIN
        cmd += ' BEGIN {minhet=1; maxhet=0; sum=0; sumsq=0}'
        ## loop NR
        cmd += ' {'
        cmd += ' het=$5;'
        cmd += ' sum+=het; sumsq+=het**2;'
        cmd += ' if(het<minhet) minhet=het; if(het>maxhet) maxhet=het}'
        ## END
        cmd += ' END {'
        cmd += ' n=NR;' ## without header
        cmd += ' mean=sum/n;'
        cmd += ' stddev=sqrt(sumsq/n - (sum/n)**2);'
        cmd += ' print mean, stddev, minhet, maxhet, n}'
        cmd += "'"
        if bool_execute == False:
            cmd += ' )'
        l_cmds += [cmd]

        if bool_execute == True:
            l = os.popen(cmd).read().split()
            average = float(l[0])
            stddev = float(l[1])
            het_min = float(l[2])
            het_max = float(l[3])
        else:
            average = None
            stddev = None
            het_min = None
            het_max = None

        cmd = 'rm %s.het.sorted %s.imiss.remove.sorted' %(self.o,self.o,)
        l_cmds += [cmd]
        if bool_execute == True:
            self.execmd(cmd)

        if bool_remove == True:
            cmd = 'rm %s.het.joined' %(self.o,)
            l_cmds += [cmd]
            if bool_execute == True:
                self.execmd(cmd)

        return l_cmds, average, stddev, het_min, het_max

    

    def plink_execution(self):

        l_plink_cmds = self.l_plink_cmds()

        ##
        ## execute plink commands
        ##
        for plink_cmd_full in l_plink_cmds:

            plink_cmd = plink_cmd_full.split()[0].replace('--','')

            ##
            ## check that file input exists (LSF file dependency does not work...)
            ##
            bool_continue_in, in_prefix, fn_in = self.check_input_existence(
                plink_cmd,plink_cmd_full,
                )
            if bool_continue_in == True:
                print('in does not exist', plink_cmd, fn_in)
                continue

            ##
            ## check that file output doesn't already exist
            ## and immediately create it if it doesn't
            ##
            bool_continue_out, out_prefix, fn_out = self.check_output_existence(
                plink_cmd,plink_cmd_full,
                )
            if bool_continue_out == True:
                print('out exists', plink_cmd, fn_out)
                continue

            ##
            ## initiate command
            ##
            cmd = ''
            if self.bool_verbose == True:
                cmd += self.mail(plink_cmd,out_prefix,'initiated',)

            ##
            ## before plink
            ##
            cmds_extra_before = self.extra_commands_before(
                plink_cmd,plink_cmd_full,
                in_prefix, out_prefix,
                )
            if cmds_extra_before != '':
                cmd += '\n%s\n' %(cmds_extra_before,)
            else:
                cmd += '\n'

            ##
            ## plink cmd
            ##
            cmd_plink = self.append_plink(plink_cmd_full,out_prefix,)
            cmd += cmd_plink

            ##
            ## after plink
            ##
            cmds_extra_after = self.extra_commands_after(
                plink_cmd,in_prefix,out_prefix,
                )
            if cmds_extra_after != '':
                cmd += '\n%s\n' %(cmds_extra_after)
            else:
                cmd += '\n'

            ##
            ## LSF settings
            ##
            cmd_LSF = self.append_LSF(
                plink_cmd,plink_cmd_full=plink_cmd_full,)
            cmd_LSF += 'bash %s_%s.sh' %(out_prefix,plink_cmd,)

            if (
                self.bool_run_all == True
                and
                plink_cmd not in ['genome','indep-pairwise',]
                ):
                cmd += self.cmd_rerun(plink_cmd,)

            ##
            ## terminate command
            ##
            if self.bool_verbose == True:
                cmd += self.mail(plink_cmd,out_prefix,'finished',)

            ##
            ## write plink command and associated commands
            ## to individual shell scripts
            ##
            fd = open('%s_%s.sh' %(out_prefix,plink_cmd,),'w')
            fd.write(cmd)
            fd.close()
            ## make shell script executable
            self.execmd('chmod +x %s_%s.sh' %(out_prefix,plink_cmd,))

            if self.bool_verbose == True:
                print()
                print(plink_cmd)
##                print cmd
##                print cmd_LSF

            if not '--remove' in cmd and plink_cmd not in [
                'missing','het','recodeA','check-sex',
                ]:
                print(cmd)
                stop1
            if '%' in cmd.replace('%02g',''):
                print(cmd)
                stop2_parenthesis_in_cmd
##            print cmd
##            print cmd_LSF
##            continue
##            stop3

            self.execmd(cmd_LSF)

        return


    def mail(self,plink_cmd,out_prefix,status,):

        msg = 'job $LSB_JOBID %s' %(status,)
        subject = 'job $LSB_JOBID %s, %s %s %s' %(
            status,self.i,plink_cmd,out_prefix,)

        cmd = ''
        address = '%s@sanger.ac.uk' %(pwd.getpwuid(os.getuid())[0],)
##        address += '%s@sanger.ac.uk\n' %(os.getlogin(),)
        cmd += '\n'
        cmd += 'echo "%s" ' %(msg)
        cmd += '| mail -s "%s" ' %(subject)
        cmd += '%s\n' %(address)

        return cmd


    def cmd_rerun(self,plink_cmd,):

        ## N.B. REMEMBER NOT TO INSERT LINE BREAKS BETWEEN FUNCTIONS HERE!!!
        ## USE SEMICOLONS INSTEAD IF MULTIPLE COMMANDS!!!

        cmd = ''
        cmd += '/software/python-3.3.2/bin/python3 %s' %(
            os.path.join(
                os.path.dirname(sys.argv[0]),
                os.path.basename(sys.argv[0]),))

        for k,v in list(vars(self.args).items()):
            if k[:4] == 'bool': ## e.g. bool_filter_females, bool_verbose, bool_no_sexcheck
                if v == True:
                    cmd += ' --%s' %(k)
                else:
                    pass
            else:
                if v == None or v == 'None':
                    pass
                else:
                    cmd += ' --%s %s' %(k,v)

        return cmd


    def check_output_existence(self,plink_cmd,plink_cmd_full,):

        bool_continue = False

        if '--out' in plink_cmd_full:
            l = plink_cmd_full.split()
            if plink_cmd in ['genome','indep-pairwise',]:
                s = l[l.index('--out')+1]
                if plink_cmd == 'genome': i2 = -len('.$i.$j')
                elif plink_cmd == 'indep-pairwise': i2 = -len('.$chrom')
                out_prefix = s[s.index('/')+1:i2]
            else:
                out_prefix = l[l.index('--out')+1]
        elif '--bfile' in plink_cmd_full:
            l = plink_cmd_full.split()
            out_prefix = os.path.basename(l[l.index('--bfile')+1])
            print(l)
            print(out_prefix)
            stop_tmp
        else:
            out_prefix = self.o

        bool_file_out_exists = False
        for out_suffix in self.d_out_suffix[plink_cmd]:
            fn_out = '%s.%s' %(out_prefix,out_suffix,)
            if os.path.isfile(fn_out):
                bool_file_out_exists = True
##                if self.bool_verbose == True:
##                    print fn_out, 'exists'
                break
            ##
            ## this conflicts with file output check in shell scripts...
            ##
##            ##
##            ## create the output so another process doesn't think it's not created
##            ##
##            else:
##                fd = open(fn_out,'w')
##                fd.close()
            continue

        ## check that no other process has started writing to the log file recently
        ## to avoid the same step being carried out twice
        ## and worst case scenario output being corrupted
        if bool_file_out_exists == False:

            if plink_cmd == 'indep-pairwise':
                chromosome = 1
                fn_log = 'prune/%s.%i.log' %(out_prefix,chromosome,)
            elif plink_cmd == 'genome':
                fn_log = 'genome/%s.00.00.log' %(out_prefix)
            else:
                fn_log = '%s.log' %(out_prefix)

            if os.path.isfile(fn_log):
##                    ## Python 2.5
##                    with open(fn_log,'r') as fd:
##                        lines = fd.readlines()
                ## All Python versions
                fd = open(fn_log,'r')
                lines = fd.readlines()
                fd.close()

                ## file is about to be generated
                if len(lines) == 0:
                    bool_file_out_exists = True
                    fn_out = fn_log
                ## file is about to be generated
                elif 'Analysis finished: ' not in lines[-2]:
                    bool_file_out_exists = True
                    fn_out = fn_log
                ## some other file was generated
                ## and wrote to a log file with the same name as the current one
                else:
                    bool_file_out_exists = False

            ##
            ## create the output so another process doesn't think it's not created
            ##
            else:

                fd = open(fn_log,'w')
                fd.close()
                pass

            pass
                
        if bool_file_out_exists == True:
            bool_continue = True
            pass

        return bool_continue, out_prefix, fn_out


    def check_input_existence(self,plink_cmd,plink_cmd_full,):

        ## start out assuming all files are present
        bool_continue = False

        ## initiate list of file paths to check
        l_fp_in = []

        ## append bfile
        if '--bfile' in plink_cmd_full:
            l = plink_cmd_full.split()
            bfile_prefix = l[l.index('--bfile')+1]
            pass
        else:
            bfile_prefix = self.i
            pass
        l_fp_in += ['%s.bed' %(bfile_prefix)]

        ## append extract, remove, exclude, keep
        for option in [
            'extract','remove','exclude','keep',
            'read-freq', 'read-genome',
            'bmerge',
            ]:
            if not '--%s' %(option) in plink_cmd_full:
                continue
            l = plink_cmd_full.split()
            if option == 'bmerge':
                fp_in = l[l.index('--%s' %(option))+2]
            else:
                fp_in = l[l.index('--%s' %(option))+1]
            l_fp_in += [fp_in]
            continue

        ## check                
        bool_input = True
        for fp_in in l_fp_in:
            bool_input = self.check_file_in(fp_in)
            if bool_input == True: continue
            if self.bool_verbose == True:
                print(plink_cmd, fp_in, 'does not exist or was recently modified or was not written to touch file')
            break
        if bool_input == False:
            bool_continue = True
            pass

        return bool_continue, bfile_prefix, fp_in


    def append_plink(self,plink_cmd_full,out_prefix,):

        ## http://pngu.mgh.harvard.edu/~purcell/plink/flow.shtml
        
        plink_cmd = plink_cmd_full.split()[0].replace('--','')

        ## initiate plink
##        cmd = 'plink \\\n'
        cmd = '%s \\\n' %(self.path_plink)

        ## --bfile
        if '--bfile' not in plink_cmd_full:
            cmd += '--bfile %s \\\n' %(self.i,)

        ##
        cmd += '%s \\\n' %(plink_cmd_full,)

        ## --extract SNPs
        if '--extract' in plink_cmd_full:
            pass
        elif plink_cmd == 'cluster':
            pass
        else:
            cmd += '--extract %s.autosomes.SNPs \\\n' %(self.o,)

        ## --exclude SNPs

        ## --remove samples

        ## --keep samples

        ## --out
        if '--out' not in cmd:
            cmd += '--out %s \\\n' %(out_prefix,)
        else:
            pass

        ## --noweb
        cmd += '--noweb \\\n'
        ## http://pngu.mgh.harvard.edu/~purcell/plink/faq.shtml#faq3
        ## --allow-no-sex
        cmd += '--allow-no-sex \\\n'
        ## http://pngu.mgh.harvard.edu/~purcell/plink/faq.shtml#faq3
        ## --nonfounders
        if plink_cmd not in ['hardy']:
            cmd += '--nonfounders \\\n'

        if plink_cmd in ['genome','indep-pairwise',]:
            ## backslash causes trouble when using bash -c 'command \ --option'
            cmd = cmd.replace('\\\n','')
            ## genome
            cmd = cmd.replace('$i','$0')
            cmd = cmd.replace('$j','$1')
            ## indep-pairwise
            cmd = cmd.replace('$chrom','$0')
        elif '$' in cmd:
            print(cmd)
            stop

        return cmd


    def append_LSF(
        self, plink_cmd, plink_cmd_full=None, JOBID=None, verbose=True,
        memMB = None,
        ):

        if plink_cmd_full:
            memMB = self.assign_memory(plink_cmd,plink_cmd_full,)
        ## --genome, EIGENSOFT, --indep-pairwise
        else:
            if memMB == None:
                stop
        
        cmd = ''
        cmd += 'bsub \\\n'
        cmd += "-M%i -R'select[mem>%i] rusage[mem=%i]' \\\n" %(
            memMB,memMB,memMB,
            )
        cmd += '-G %s \\\n' %(self.project)
        cmd += '-q normal \\\n'
        ## redirect stdout and stderr
        if verbose == True and self.bool_verbose == True:
            cmd += '-o stdout/%s.%s.out \\\n' %(self.o,plink_cmd,) ## tmp
            cmd += '-e stderr/%s.%s.err \\\n' %(self.o,plink_cmd,) ## tmp
        else:
            cmd += '-o /dev/null \\\n'
        ## specify JOBID
        if not JOBID:
            JOBID = '%s.%s' %(self.o,plink_cmd,)
##        cmd += '-J"%s" \\\n' %(JOBID,)
        cmd += "-J'%s' \\\n" %(JOBID,)

        return cmd


    def assign_memory(
        self,plink_cmd,plink_cmd_full=None,memMB_per_1000_samples=None,):

        ## approximate memory (MB) required per 1000 samples
        ## depends on number of SNPs as well...
        ## numbers below worked out for 2.5M SNPs
        d_memMB = {
            'check-sex':800,
            'cluster':900,
##            'het':1100,
            'recodeA':800,
            'bmerge':1200,
            'freq':800,
            'indep-pairwise':800,
            'make-bed':800,
            'missing':800,
            'genome':800,
            'hardy':800,
            }

        cmd = 'cat %s.fam | wc -l' %(self.i)
        n_samples = int(os.popen(cmd).read())
        ## additional samples if bmerge
        if plink_cmd == 'bmerge':
            cmd = 'cat %s.fam | wc -l' %(self.fp1000g)
            if self.bool_verbose == True:
                print(cmd)
            n_samples += int(os.popen(cmd).read())

        ## assign appropriate amount of memory (GB)
        if not memMB_per_1000_samples:
            if not plink_cmd in list(d_memMB.keys()):
                memMB_per_1000_samples = 800
            else:
                memMB_per_1000_samples = d_memMB[plink_cmd]
            
        memGB_calc = math.ceil(memMB_per_1000_samples*n_samples/1000000.)
        memMB = int(max(2000,1000*memGB_calc))

        return memMB


    def l_plink_cmds(self):

        l_plink_cmds = []

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune
        s_indep_pairwise = '--indep-pairwise %i %i %f --maf %f \\\n' %(
            int(self.args.threshold_indepWindow),
            int(self.args.threshold_indepShift),
            float(self.args.threshold_indepRsquared),
            float(self.args.threshold_indepMAF),
            )

        ## http://pngu.mgh.harvard.edu/~purcell/plink/strat.shtml#cluster
        s_cluster = '--cluster \\\n'
        s_cluster += '--mds-plot 10 \\\n' ## 10 components

##        ##
##        ## sequential for Manj
##        ##
##        l_plink_cmds += ['--het --remove %s.imiss.samples --out %s.sequential' %(self.o,self.o,)]
##        l_plink_cmds += ['--check-sex --extract %s.X.SNPs --remove %s.imiss.samples --out %s.sequential' %(self.o,self.o,self.o,)]

        ##
        ## 1) sample QC
        ##

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#missing
        l_plink_cmds += ['--missing']

        ## http://pngu.mgh.harvard.edu/~purcell/plink/ibdibs.shtml#inbreeding
##        cmd = '--het'
        cmd = '--recodeA'
        ## keep just to make sure it doesn't run until --missing finishes
        ## otherwise correct SD can't be calculated
        cmd += ' --keep %s.imiss' %(self.o)
        l_plink_cmds += [cmd]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#sexcheck
        if self.bool_no_sexcheck == False:
            l_plink_cmds += ['--check-sex --extract %s.X_not_XY.SNPs' %(self.o)]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed
        s = '--make-bed \\\n'
        s += '--remove %s.sampleQC.samples \\\n' %(self.o)
        s += '--out %s.sampleQC \\\n' %(self.o)
        l_plink_cmds += [s]

        ##
        ## 2) SNP QC, pre Hardy, SNP exclusion (missingness)
        ##

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#missing
        s = '--missing \\\n'
        s +='--out %s.SNPQC \\\n' %(self.o)
        s += '--remove %s.sampleQC.samples \\\n' %(self.o,)
        l_plink_cmds += [s]

        ##
        ## 3) SNP QC, pre Hardy, sample removal (IBD/relatedness)
        ##

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#freq
        s = '--freq \\\n'
        s += '--remove %s.sampleQC.samples \\\n' %(self.o)
        s += '--out %s.preIBD \\\n' %(self.o)
        s += '--exclude %s.lmiss.SNPs \\\n' %(self.o)
##        s += '--exclude %s.lmiss.exclLRLD.SNPs \\\n' %(self.o)
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune
        s = s_indep_pairwise
        s += '--read-freq %s.preIBD.frq \\\n' %(self.o)
        s += '--remove %s.sampleQC.samples \\\n' %(self.o)
        s += '--out prune/%s.prehardy.$chrom \\\n' %(self.o)
        s += '--chr $chrom \\\n' ## parallel
        s += '--exclude %s.lmiss.SNPs \\\n' %(self.o)
##        s += '--exclude %s.lmiss.exclLRLD.SNPs \\\n' %(self.o)
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/ibdibs.shtml#genome
        s = self.write_genome_cmd(
            suffix='prehardy', suffix_remove='sampleQC', bfile_in=self.i,)
        l_plink_cmds += [s]

##        ## http://pngu.mgh.harvard.edu/~purcell/plink/strat.shtml#cluster
##        s = s_cluster+'--read-genome %s.prehardy.genome \\\n' %(self.o,)
##        s += '--bfile %s \\\n' %(self.i,)
##        s += '--remove %s.sampleQC.samples \\\n' %(self.o,)
##        s += '--extract %s.prehardy.prune.in \\\n' %(self.o,)
##        s += '--out %s.prehardy \\\n' %(self.o,)
##        l_plink_cmds += [s]

        ##
        ## 4) SNP QC, Hardy, SNP exclusion (HWE)
        ##

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#hardy
        s = '--hardy \\\n'
        s += '--remove %s.sampleQC.nonfounders.samples \\\n' %(self.o,)
        s += '--exclude %s.lmiss.SNPs \\\n' %(self.o,)
##        s += '--exclude %s.lmiss.exclLRLD.SNPs \\\n' %(self.o)
        l_plink_cmds += [s]

        ##
        ## 5a) SNP QC, post Hardy, PCA without 1000G
        ## 5b) SNP QC, post Hardy, MDS without 1000G
        ##

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#freq
        s = '--freq \\\n'
        s += '--out %s.postIBD \\\n' %(self.o)
        s += '--remove %s.sampleQC.IBD.samples \\\n' %(self.o,)
        s += '--exclude %s.lmiss.hwe.LRLD.SNPs \\\n' %(self.o,)
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune
        s = s_indep_pairwise
        s += '--read-freq %s.postIBD.frq \\\n' %(self.o,)
        s += '--remove %s.sampleQC.IBD.samples \\\n' %(self.o,)
        s += '--exclude %s.lmiss.hwe.LRLD.SNPs \\\n' %(self.o,)
        s += '--out prune/%s.posthardy.$chrom \\\n' %(self.o)
        s += '--chr $chrom \\\n' ## parallel
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/ibdibs.shtml#genome
        s = self.write_genome_cmd(
            suffix='posthardy', suffix_remove='sampleQC.IBD', bfile_in=self.i,
            )
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/strat.shtml#cluster
        s = s_cluster+'--read-genome %s.posthardy.genome \\\n' %(self.o,)
        s += '--bfile %s \\\n' %(self.i,)
        s += '--remove %s.sampleQC.nonfounders.samples \\\n' %(self.o,)
        s += '--extract %s.posthardy.prune.in \\\n' %(self.o,)
        s += '--out %s.posthardy \\\n' %(self.o)
        l_plink_cmds += [s]

        ##
        ## 5c) SNP QC, post Hardy, MDS with 1000G
        ## 5d) SNP QC, post Hardy, PCA with 1000G
        ##

        suffix = self.fn1000g

        ## 5) http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#merge
        s = '--bmerge \\\n'
        s += '%s.bed \\\n' %(self.fp1000g,)
        s += '%s.bim \\\n' %(self.fp1000g,)
        s += '%s.fam \\\n' %(self.fp1000g,)
        s += '--make-bed \\\n'
        s += '--out %s.%s \\\n' %(self.o,suffix,)
        s += '--extract %s.%s.autosomes.comm.SNPs \\\n' %(self.o,self.fn1000g,)
        s += '--bfile %s \\\n' %(self.i,)
        s += '--remove %s.sampleQC.IBD.samples \\\n' %(self.o)
        s += '--exclude %s.lmiss.hwe.SNPs \\\n' %(self.o,)
##        s += '--exclude %s.lmiss.hwe.monomorphic.SNPs \\\n' %(self.o,)
        l_plink_cmds += [s]

        ## 5) http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#freq
        s = '--freq \\\n'
        s += '--bfile %s.%s \\\n' %(self.o,suffix,)
        s += '--out %s.%s \\\n' %(self.o,suffix,)
        s += '--extract %s.%s.autosomes.comm.SNPs \\\n' %(self.o,suffix,)
        s += '--remove %s.sampleQC.IBD.samples \\\n' %(self.o,)
        s += '--exclude %s.ldregions.SNPs \\\n' %(self.o,)
        l_plink_cmds += [s]

        ## 5) http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune
        ## followed by PCA; see function indep_pairwise_after
        s = s_indep_pairwise
        s += '--bfile %s.%s \\\n' %(self.o,suffix,)
        s += '--read-freq %s.%s.frq \\\n' %(self.o,suffix,)
        s += '--extract %s.%s.autosomes.comm.SNPs \\\n' %(self.o,suffix,)
        s += '--remove %s.sampleQC.IBD.samples \\\n' %(self.o,)
        s += '--exclude %s.ldregions.SNPs \\\n' %(self.o,)
        s += '--out prune/%s.%s.$chrom \\\n' %(self.o,suffix,)
        s += '--chr $chrom \\\n'
        l_plink_cmds += [s]

##        ## 5) http://pngu.mgh.harvard.edu/~purcell/plink/ibdibs.shtml#genome
##        s = self.write_genome_cmd(
##            ## posthardy
##            suffix=suffix, suffix_remove='sampleQC.IBD',
##            bfile_in='%s.%s' %(self.o,suffix,),
####            ## prehardy
####            bfile, suffix, 'genome', '%s.%s' %(bfile,suffix,),
##            )
##        l_plink_cmds += [s]
##
##        ## 5) http://pngu.mgh.harvard.edu/~purcell/plink/strat.shtml#cluster
##        s = s_cluster+'--read-genome %s.%s.genome \\\n' %(self.o,suffix,)
##        s += '--bfile %s.%s \\\n' %(self.o,suffix,)
##        s += '--extract %s.%s.prune.in \\\n' %(self.o,suffix,)
##        s += '--out %s.%s \\\n' %(self.o,suffix,)
##        ## posthardy
##        s += '--remove %s.sampleQC.nonfounders.samples \\\n' %(self.o,)
####        ## prehardy
####        s += '--remove %s.sampleQC.IBD.samples \\\n' %(self.o,)
##        l_plink_cmds += [s]

        ##
        ## 6) post SNP QC, make-bed
        ##

        ## http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed
        s = '--make-bed \\\n'
        s += '--remove %s.sampleQC.IBD.samples \\\n' %(self.o)
        s += '--out %s.postQC.autosomes \\\n' %(self.o)
        s += '--exclude %s.lmiss.hwe.SNPs \\\n' %(self.o,)
        l_plink_cmds += [s]

        ##
        ## 7) SNP QC, chromosome X
        ##

        if self.bool_no_sexcheck == False:

            for sex in ['males','females']:

                if sex == 'males' and self.bool_filter_females == True: continue

                ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#missing
                s = '--missing \\\n'
                s += '--out %s.X.%s \\\n' %(self.o,sex,)
                s += '--remove %s.sampleQC.IBD.samples \\\n' %(self.o,)
                s += '--extract %s.X_and_XY.SNPs \\\n' %(self.o)
                s += '--filter-%s \\\n' %(sex)
                l_plink_cmds += [s]

            ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#hardy
            s = '--hardy \\\n'
            s += '--out %s.X.females \\\n' %(self.o,)
            s += '--remove %s.sampleQC.nonfounders.samples \\\n' %(self.o,)
    ##        s += '--remove %s.sampleQC.IBD.samples \\\n' %(self.o,)
            s += '--extract %s.X_and_XY.SNPs \\\n' %(self.o)
            s += '--exclude %s.X.lmiss.union.SNPs \\\n' %(self.o,)
            s += '--filter-females \\\n'
            l_plink_cmds += [s]

            ## http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed
            s = '--make-bed \\\n'
            s += '--remove %s.sampleQC.IBD.samples \\\n' %(self.o)
            s += '--extract %s.X_and_XY.SNPs \\\n' %(self.o)
            s += '--out %s.postQC.X \\\n' %(self.o,)
            s += '--exclude %s.X.lmiss.hwe.SNPs \\\n' %(self.o,)
            l_plink_cmds += [s]

        return l_plink_cmds


    def write_genome_cmd(
        self, suffix, suffix_remove, bfile_in,
        ):

        out = 'genome/%s.%s.$i.$j' %(self.o,suffix,),
        remove = '%s.%s.samples' %(self.o,suffix_remove,),
        extract = '%s.%s.prune.in' %(self.o,suffix,),
        genome_lists = 'fam/%s.%s.fam.$i fam/%s.%s.fam.$j' %(
            self.o,suffix,self.o,suffix,
            )

        cmd_genome = '--genome \\\n'
        cmd_genome += '--genome-lists %s \\\n' %(genome_lists)
        cmd_genome += '--bfile %s \\\n' %(bfile_in)
        cmd_genome += '--out %s \\\n' %(out)
        cmd_genome += '--remove %s \\\n' %(remove)
        cmd_genome += '--extract %s' %(extract)

        return cmd_genome


    def check_file_in(self,fp_in,):

        bool_check = True
        ## file exists?
        if not os.path.isfile(fp_in):
            bool_check = False
        ## file genereation is finished; i.e. written to touch?
        else:
            fd = open('%s.touch' %(self.o),'r')
            s = fd.read()
            fd.close()
            l_fp_out = s.strip().split('\n')
            if not fp_in in l_fp_out:
                bool_check = False

        return bool_check


    def extra_commands_before(
        self,plink_cmd,plink_cmd_full,in_prefix,out_prefix,
        ):

        l_cmds = []

        cmd = ''
        for out_suffix in self.d_out_suffix[plink_cmd]:
            fn_out = '%s.%s' %(out_prefix,out_suffix,)
            cmd += 'if [ -f %s ]; then\nexit\nfi\n' %(fn_out)
        l_cmds += [cmd]

        l_cmds += ['touch %s\n' %(fn_out)]

        if plink_cmd == None:
            pass

        elif plink_cmd == 'genome':
            l_cmds = self.genome_before(in_prefix,out_prefix,)

        elif plink_cmd == 'indep-pairwise': ## parallel
            l_cmds = self.indep_pairwise_before(out_prefix,)

##        elif plink_cmd == 'hardy':
##            l_cmds = self.hardy_before(out_prefix,)

        cmds = '\n\n'.join(l_cmds)
        
        return cmds


    def extra_commands_after(self,plink_cmd,in_prefix,out_prefix,):

        l_cmds = []

        if plink_cmd == 'missing':
            l_cmds = self.missing_after()

        elif plink_cmd == 'recodeA':
            l_cmds = self.recodeA_after()
            cmd = self.concatenate_sampleQC_remove_lists()
            l_cmds += [cmd]

        elif plink_cmd == 'check-sex':
            cmd = 'echo %s.sexcheck >> %s.touch\n' %(self.o,self.o,)
            l_cmds += [cmd]
            cmd = '''awk '{if($5=="PROBLEM") print $1,$2}' '''
            cmd += ' %s.sexcheck > %s.sexcheck.samples' %(self.o, self.o,)
            l_cmds += [cmd]
            cmd = 'echo %s.sexcheck.samples >> %s.touch\n' %(self.o,self.o,)
            l_cmds += [cmd]
            cmd = self.concatenate_sampleQC_remove_lists()
            l_cmds += [cmd]

        elif plink_cmd == 'genome':
            l_cmds = self.genome_after(in_prefix,out_prefix,)

        elif plink_cmd == 'indep-pairwise':
            l_cmds = self.indep_pairwise_after(out_prefix,in_prefix,)

        elif plink_cmd == 'hardy':
            l_cmds = self.hardy_after(out_prefix,)

        elif plink_cmd == 'freq':
            cmd = 'echo %s.frq >> %s.touch' %(out_prefix,self.o,)
            l_cmds += [cmd]

        elif plink_cmd == 'bmerge':
            cmd = 'echo %s.%s.bed >> %s.touch' %(self.o,self.fn1000g,self.o,)
            l_cmds += [cmd]
##            l_cmds = self.bmerge_after(bfile,out_prefix,)

##        elif plink_cmd == 'cluster':
##            l_cmds += ['fi\n']

        elif plink_cmd == 'cluster':
            cmd = 'echo %s.mds >> %s.touch' %(out_prefix,self.o,)
            l_cmds += [cmd]

        cmds = '\n\n'.join(l_cmds)

        return cmds


    def hardy_after(self,out_prefix,):

        l_cmds = []

        ##
        ## autosome
        ##
        if out_prefix == self.o:

            ##
            ## if output exists
            ##
            cmd = 'if [ -s %s.hwe ]; then\n' %(self.o,)

            cmd += 'echo %s.hwe >> %s.touch\n' %(self.o,self.o,)

            cmd += "awk '{if ($9 < %.1e) print $2}' %s.hwe > %s.hwe.SNPs\n" %(
                float(self.args.threshold_hwe_min), self.o,self.o,)
            cmd += 'echo %s.hwe.SNPs >> %s.touch\n' %(self.o,self.o,)

            cmd += 'cat %s.lmiss.SNPs %s.hwe.SNPs > %s.lmiss.hwe.SNPs\n' %(
                self.o,self.o,self.o,
                )
            cmd += 'echo %s.lmiss.hwe.SNPs >> %s.touch\n' %(self.o,self.o,)

            cmd += 'cat %s.lmiss.SNPs %s.hwe.SNPs %s.ldregions.SNPs' %(
                self.o,self.o,self.o,)
            cmd += ' > %s.lmiss.hwe.LRLD.SNPs\n' %(self.o)
            cmd += 'echo %s.lmiss.hwe.LRLD.SNPs >> %s.touch\n' %(self.o,self.o,)

            cmd += 'cat %s.lmiss.SNPs %s.hwe.SNPs' %(self.o,self.o)
            cmd += " <(awk 'NR>1&&$5<=0{print $2}' %s.postIBD.frq)" %(self.o)
            cmd += ' > %s.lmiss.hwe.monomorphic.SNPs\n' %(self.o)
            cmd += 'echo %s.lmiss.hwe.monomorphic.SNPs >> %s.touch\n' %(self.o,self.o,)

            ##
            ## only merge SNPs found in 1000g and current dataset
            ##
            for fn_in,fn_out in [
                [self.fp1000g,self.fn1000g,],
                [self.i,self.o,],
                ]:
                ## 1) sort
                cmd += "\ncat %s.bim | awk '{if($1>=1&&$1<=22) print $2}' " %(fn_in)
                cmd += '| sort > %s.autosomes.SNPs.sorted' %(fn_out,)

            fn_out = '%s.%s.autosomes.comm.SNPs' %(self.o,self.fn1000g,)
            cmd += '\ncomm -12 %s.autosomes.SNPs.sorted %s.autosomes.SNPs.sorted > %s' %(
                self.o,self.fn1000g,fn_out,)
            cmd += '\necho %s >> %s.touch' %(fn_out,self.o,)

            cmd += '\nrm %s.autosomes.SNPs.sorted %s.autosomes.SNPs.sorted' %(self.o,self.fn1000g,)

            cmd += '\necho %s.hwe >> %s.touch\n' %(out_prefix,self.o,)
            
            ##
            ## fi
            ##
            cmd += '\nfi\n'
            l_cmds += [cmd]

        ##
        ## X
        ##
        if out_prefix == '%s.X.females' %(self.o):

            cmd = 'if [ -s %s.X.females.hwe' %(self.o)
            if self.bool_filter_females == False:
                cmd += ' -a -s %s.X.males.lmiss.SNPs' %(self.o)
            cmd += ' ]; then\n'

            cmd += 'cat %s.X.females.hwe | ' %(self.o)
            cmd += "awk '{if ($9 < %.1e) print $2}'" %(
                float(self.args.threshold_hwe_min))
            cmd += ' > %s.X.females.hwe.SNPs\n' %(self.o)
            cmd += 'cat '
            if self.bool_filter_females == False:
                cmd += '%s.X.males.lmiss.SNPs ' %(self.o)
            cmd += '%s.X.females.lmiss.SNPs ' %(self.o)
            cmd += '%s.X.females.hwe.SNPs ' %(self.o)
            cmd += '| sort -u '
            cmd += '> %s.X.lmiss.hwe.SNPs\n' %(self.o)

            cmd += 'echo %s.X.lmiss.hwe.SNPs >> %s.touch\n' %(self.o,self.o,)

            cmd += 'fi\n'
            l_cmds += [cmd]

        return l_cmds


    def EIGENSOFT(self,prefix_in,prefix_out,bool_removal=True,):

        prefix = prefix_out

        ## http://computing.bio.cam.ac.uk/local/doc/eigenstrat.txt
        ## http://helix.nih.gov/Applications/README.eigenstrat

        ##Estimated running time of the smartpca program is 
        ##  2.5e-12 * nSNP * NSAMPLES^2 hours            if not removing outliers.
        ##  2.5e-12 * nSNP * NSAMPLES^2 hours * (1+m)    if m outlier removal iterations.
        ##Thus, under the default of up to 5 outlier removal iterations, running time is 
        ##  up to 1.5e-11 * nSNP * NSAMPLES^2 hours.

        memMB = 4000

        if self.project == 'uganda_gwas':
            bool_removal = False ## ms23 24jan2013
        elif self.project == 'agv':
            bool_removal = False ## cp8 17jan2013
        else:
            bool_removal = False

        self.write_EIGENSOFT_parameter_file(prefix,bool_removal,)

        cmd = ''
        cmd += 'echo EIGENSOFT\n\n'

        ## EIGENSOFT has a 39 character limit for idnames; i.e. 2x19+1
        cmd += 'cat %s.fam' %(prefix_in)
        ## init awk
        cmd += " | awk '{"
        cmd += 's1=$1;s2=$2;'
        cmd += 'sub(/-RECLUSTER/,"",s1);sub(/-RECLUSTER/,"",s2);'
        cmd += 'sub(/-RESCAN/,"",s1);sub(/-RESCAN/,"",s2);'
        cmd += 'print $1,$2,'
        cmd += 'substr(s1,length(s1)-18,19),substr(s2,length(s2)-18,19)'
        ## term awk
        cmd += "}'"
        cmd += ' > %s.EIGENSOFT.recoded.txt' %(prefix)
        cmd += '\n\n'

        ##
        ## exit if duplicate present (or created)
        ##
        cmd += 'dups=$('
        cmd += 'cat %s.EIGENSOFT.recoded.txt' %(prefix)
        cmd += " | awk '{print $3}'"
        cmd += ' | sort'
        cmd += ' | uniq -d'
        cmd += ' | wc -l'
        cmd += ')\n'
        cmd += 'if [ $dups -gt 0 ]; then\n'
        cmd += 'echo duplicate sample IDs or bug\n'
        cmd += 'exit\n'
        cmd += 'fi\n'
        cmd += '\n'

        ##
        ## --pheno
        ##
        for s1,s2,s3 in [
            [' -v2','1','>',], ## fou
            ['','2','>>',], ## nof
            ]:
            cmd += 'join -1 2 -2 2 %s -o 2.1,2.2,2.3,2.4,2.5,2.6' %(s1)
            ## sort non-founders
            cmd += ' <(sort -k2,2 %s.sampleQC.nonfounders.samples)' %(self.o,)
            ## sort all samples
            cmd += ' <(sort -k2,2 %s.fam)' %(prefix_in,)
            cmd += ''' | awk '{print $1,$2,"%s"}' ''' %(s2)
            cmd += ' %s %s.EIGENSOFT.phe\n\n' %(s3,prefix,)

##        cmd += 'if [ ! -f %s.EIGENSOFT.bed ]; then\n\n' %(prefix)
        cmd += '\n%s \\\n' %(self.path_plink)
        cmd += '--bfile %s \\\n' %(prefix_in,)
        cmd += '--make-bed \\\n'
        cmd += '--noweb \\\n'
        cmd += '--extract %s.prune.in \\\n' %(prefix_out,)
        cmd += '--remove %s.sampleQC.IBD.samples \\\n' %(self.o)
##        cmd += '--out %s.tmp \\\n' %(prefix)
        cmd += '--out %s.EIGENSOFT \\\n' %(prefix)
        cmd += '--no-pheno --pheno %s.EIGENSOFT.phe \\\n' %(prefix) ## does not work for strings (e.g. fou and nof)
##        cmd += '\nfi'
        cmd += '\n\n'

        cmd += 'rm %s.EIGENSOFT.phe\n\n' %(prefix)

        ##
        ## --update-ids
        ##
        cmd += 'join -1 3 -2 2 -o 1.1,2.3,2.4,1.4,1.5,1.6,1.7'
        cmd += " <(awk '{print NR,$0}' %s.EIGENSOFT.fam | sort -k3,3)" %(prefix)
        cmd += ' <(sort -k2,2 %s.EIGENSOFT.recoded.txt)' %(prefix)
        cmd += ' | sort -k1n,1 '
        cmd += """ | awk '{if($7==1) {pheno="fou"} else {pheno="nof"}; print $2,$3,$4,$5,$6,pheno}'"""
        cmd += ' > %s.EIGENSOFT.updateids.fam \n\n' %(prefix)

        cmd += 'mv %s.EIGENSOFT.updateids.fam %s.EIGENSOFT.fam\n\n' %(prefix,prefix)

##        ## --update-ids cannot be combined with --remove
##        cmd += '%s \\\n' %(self.path_plink)
##        cmd += '--bfile %s.tmp \\\n' %(prefix,)
##        cmd += '--make-bed \\\n'
##        cmd += '--noweb \\\n'
##        cmd += '--update-ids %s.EIGENSOFT.recoded.txt \\\n' %(prefix,)
##        cmd += '--out %s.EIGENSOFT \\\n' %(prefix)
##        cmd += '\n\n'

##        if False: ## tmp!!!
##            cmd += 'rm %s.tmp.fam %s.tmp.log %s.tmp.bim %s.tmp.bed\n\n' %(
##                prefix,prefix,prefix,prefix,)

        ##
        ## create population list file
        ##
        cmd += 'echo "fou" > %s.EIGENSOFT.poplist\n\n' %(prefix)

        ##
        ## run EIGENSOFT with parameter file options (with or without sample removal)
        ##
        cmd += '\n%s \\\n' %(self.eigensoft)
        cmd += '-p %s.par\n\n' %(prefix)

        ##
        ## parse outliers (PROBABLY BUGGED AFTER MAKING CHANGES ABOVE!!!)
        ##
        if bool_removal == True:
            ## print EIGENSOFT sample names to file
            cmd += 'cat %s.EIGENSOFT.outlier' %(prefix)
            cmd += " | awk '{print $3}'"
            cmd += ''' | awk 'BEGIN{FS=":"}{print $2}' '''
            cmd += ' > %s.outlier.EIGENSOFTsamples' %(prefix)
            ## convert EIGENSOFT sample names to original sample names
            cmd += ' ; fgrep -w -f %s.outlier.EIGENSOFTsamples %s.EIGENSOFT.recoded.txt' %(prefix,prefix,)
            cmd += " | awk '{print $1,$2}'"
            cmd += ' > %s.EIGENSOFT.samples' %(prefix)

            cmd += '\n\n'

            ## clean-up
            cmd += 'rm %s.outlier.EIGENSOFTsamples\n\n' %(prefix)

        else:

            fd = open('%s.EIGENSOFT.samples' %(prefix),'w')
            fd.close()

        ## signal that we are done
        cmd += 'echo %s.EIGENSOFT.evec >> %s.touch\n\n' %(prefix,self.o)

##        ##
##        ## append outliers
##        ##
##        ## after hardy
##        cmd += ' cat %s.sampleQC.IBD.samples %s.EIGENSOFT.samples' %(bfile,bfile)
##        cmd += ' > %s.sampleQC.IBD.EIGENSOFT.samples' %(bfile)

        fn_sh = '%s_EIGENSOFT.sh' %(prefix)
        fd = open(fn_sh,'w')
        fd.write(cmd)
        fd.close()
        self.execmd('chmod +x %s' %(fn_sh))
        cmd_LSF = self.append_LSF(prefix,'EIGENSOFT',memMB=memMB,)
        cmd_LSF += 'bash %s' %(fn_sh)

        return cmd_LSF


##    def cluster_after(self,bfile,in_prefix,out_prefix,):
##
##        l_cmds = []
##        cmd = '\n;\n'
##
##        cmd += '"'
##        l_cmds += [cmd]
##
##        return l_cmds


    def write_EIGENSOFT_parameter_file(self,prefix,bool_removal,):

        ## http://helix.nih.gov/Applications/README.popgen
        ## http://computing.bio.cam.ac.uk/local/doc/popgen.txt
        ##
        ## input
        par = 'genotypename: %s.EIGENSOFT.bed\n' %(prefix,)
        par += 'snpname: %s.EIGENSOFT.bim\n' %(prefix,)
        par += 'indivname: %s.EIGENSOFT.fam\n' %(prefix,)
        ##
        ## output
        par += 'evecoutname: %s.EIGENSOFT.evec\n' %(prefix,)
        par += 'evaloutname: %s.EIGENSOFT.eval\n' %(prefix,)
        ##
        ## optional parameters
        ##
        ## numoutlieriter: maximum number of outlier removal iterations.
        ## Default is 5.  To turn off outlier removal, set this parameter to 0.
        if bool_removal == False:
            par += 'numoutlieriter: 0\n'
        ##
        ## outliersigmathresh: number of standard deviations which an individual must 
        ## exceed, along one of the top (numoutlierevec) principal components, in
        ## order for that individual to be removed as an outlier.  Default is 6.0.
##            par += 'outliersigmathresh: 6\n'
        ##
        ## outlieroutname: output logfile of outlier individuals removed. If not specified,
            ## smartpca will print this information to stdout, which is the default.
            par += 'outlieroutname: %s.EIGENSOFT.outlier\n' %(prefix)
        ##
        ## numoutlierevec: number of principal components along which to
        ## remove outliers during each outlier removal iteration.  Default is 10.
##            par += 'numoutlierevec: 2\n'
        ##
        ## qtmode: If set to YES, assume that there is a single population and that the
        ## population field contains information on real-valued phenotypes.
        ## The default value for this parameter is NO.
##        par += 'qtmode: YES\n'
        ##
        ## http://genepath.med.harvard.edu/~reich/eigensoftfaq.htm
        ## Question: How do I compute principal components using only a subset of populations and project other populations onto those principal components? 
        ## Answer: Use -w flag to smartpca.perl (see EIGENSTRAT/README), or use poplistname parameter to smartpca (see POPGEN/README).
        ##
        ## -w poplist       : compute eigenvectors using populations in poplist only,
        ##                   where poplist is an ASCII file with one population per line
        ##poplistname:   If wishing to infer eigenvectors using only individuals from a 
        ##  subset of populations, and then project individuals from all populations 
        ##  onto those eigenvectors, this input file contains a list of population names,
        ##  one population name per line, which will be used to infer eigenvectors.  
        ##  It is assumed that the population of each individual is specified in the 
        ##  indiv file.  Default is to use individuals from all populations.
        ##
        ## http://helix.nih.gov/Applications/README.convertf
        ##     6th column is case/control status (1 is control, 2 is case) OR
        ##      quantitative trait value OR population group label.
        par += 'poplistname: %s.EIGENSOFT.poplist\n' %(prefix)
        ##
        par += 'snpweightoutname: %s.EIGENSOFT.snpweight\n' %(prefix)

        fd = open('%s.par' %(prefix),'w')
        fd.write(par)
        fd.close()

        return


    def indep_pairwise_before(self,out_prefix,):

        l_cmds = []

        ## loop over chromosomes
        cmd = '\n##\n## loop 1 (plink execution)\n##\n'
        cmd += 'for chrom in {1..22}\ndo\n' ## magic number
        cmd += '## continue if output exists\n'
        cmd += 'if ['
        cmd += ' -f prune/%s.$chrom.prune.in -o' %(out_prefix)
        cmd += ' -s prune/%s.$chrom.log' %(out_prefix)
        cmd += ' ]\nthen continue\nfi\n\n'
        ## initiate command
        cmd += "cmd='\n"
        l_cmds += [cmd]

        return l_cmds


    def genome_before(self,in_prefix,out_prefix,):

        l_cmds = []

        if out_prefix == '%s.prehardy' %(self.o):
            fn_samples_remove = '%s.sampleQC.samples' %(self.o)
        elif out_prefix == '%s.posthardy' %(self.o):
            fn_samples_remove = '%s.sampleQC.IBD.samples' %(self.o)
        elif out_prefix == '%s.%s' %(self.o,self.fn1000g,):
            fn_samples_remove = '%s.sampleQC.IBD.samples' %(self.o)
        else:
            print(in_prefix)
            print(out_prefix)
            print(self.fn1000g)
            stop

        ##
        ## generate family lists
        ##

        ## sort sample removal file
        cmd = 'if ['
        cmd += ' ! -s %s.sorted' %(fn_samples_remove,)
        cmd += ' -a ! -s fam/%s.fam.00' %(out_prefix,)
        cmd += ' ]; then\n'
        cmd += 'parent=1\n' ## one time only
        cmd += 'cat %s | sort -k1,1 -k2,2 ' %(fn_samples_remove,)
        cmd += '> %s.sorted' %(fn_samples_remove,)
        cmd += '\nelse\nparent=0'
        cmd += '\nfi'
        l_cmds += [cmd]
        
        ## sort fam file
        cmd = 'if [ $parent -eq 1 ]; then\n' ## one time only
        cmd += "cat %s.fam | sort -k1,1 -k2,2" %(in_prefix,)
        cmd += " | awk '{print $1,$2}' > %s.fam.sorted" %(out_prefix,)
        cmd += '\nfi'
        l_cmds += [cmd]

        ## exclude samples
        cmd = 'if [ $parent -eq 1 ]; then\n' ## one time only
        cmd += "comm -23 %s.fam.sorted %s.sorted" %(out_prefix,fn_samples_remove,)
        cmd += " > %s.filtered.fam" %(out_prefix,)
        cmd += '\nfi'
        l_cmds += [cmd]

        ## parse FID and IID columns
        cmd = 'if [ $parent -eq 0 ]; then\n' ## one time only
        ## sleep while allowing parent to split fam files
        cmd += 'sleep 60'
        cmd += '\nelse\n'
        cmd += 'cat %s.filtered.fam' %(out_prefix,)
        cmd += ' | '
        cmd += " awk '{print $1,$2}' "
        ## and pipe them to individual fam files
        cmd += ' | split -d -a 2 -l %i - fam/%s.fam.' %(
            self.len_lists,out_prefix,
            )
        cmd += '\nfi'
        l_cmds += [cmd]

        ## clean up
        cmd = 'if [ $parent -eq 1 ]; then\n' ## one time only
        cmd += 'rm %s.sorted\n' %(fn_samples_remove)
        cmd += 'rm %s.fam.sorted\n' %(out_prefix)
        cmd += 'rm %s.filtered.fam\n' %(out_prefix)
        cmd += '\nfi'
        l_cmds += [cmd]

        ##
        ## count number of fam files generated
        ##
        cmd = 'n=$(cat %s.fam | wc -l)' %(in_prefix,)
        l_cmds += [cmd]
        if out_prefix != '%s.%s' %(self.o,self.fn1000g,):
            cmd = 'let n=$n-$(cat %s | wc -l)' %(fn_samples_remove,)
            l_cmds += [cmd]
        cmd = 'cnt=$(((${n}+%i-1)/%i))' %(
            self.len_lists, self.len_lists,
            )
        l_cmds += [cmd]

        ## loop over combinations of fam files
        cmd = '\n##\n## loop 1 (plink execution)\n##\n'
        cmd += 'for i in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
        cmd += 'for j in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
        cmd += 'if [ $j -gt $i ]\nthen continue\nfi\n\n'
        cmd += '## continue if output exists\n'
        cmd += 'if [ -s genome/%s.$i.$j.genome ]\nthen continue\nfi\n\n' %(
            out_prefix,
            )
        ## initiate command
        cmd += "cmd='\n"
        l_cmds += [cmd]

        return l_cmds


    def indep_pairwise_after(self,out_prefix,in_prefix,):
        l_cmds = []

        memMB = self.assign_memory(
            'indep-pairwise',
            memMB_per_1000_samples=150, ## based on chromosome 1
            )

        ##
        ## continuation of nested loop and command from indep_pairwise_before
        ##
        cmd = ''

        ## append to command
        if self.bool_run_all == True:
            cmd += ';'
            ## need to rerun at the end of each bsub
            cmd += self.cmd_rerun('indep-pairwise')
            pass
        ## terminate command
        cmd += "\n'\n\n"
        ## evaluate command
        cmd_LSF = self.append_LSF(
            'indep-pairwise',
            JOBID='%s.indep-pairwise.${chrom}' %(out_prefix),
            memMB = memMB,
            )
##        cmd += 'echo $cmd\n'
        cmd += '''cmd="%sbash -c '$cmd' $chrom"\n''' %(cmd_LSF)
        cmd += 'echo $cmd\n'
        cmd += 'eval $cmd\n'
        ## end loops
        cmd += '\ndone\n'
        ## exit before count if file out exists
        cmd += 'if [ -f %s.prune.in ];then\nexit\nfi\n' %(out_prefix)
        l_cmds += [cmd]

        ## count observed lines
        l_fn = [
            'prune/%s.$chrom.prune.in' %(out_prefix),
            'prune/%s.$chrom.prune.out' %(out_prefix),
            ]
        cmd = '\n##\n## loop 2 (count lines)\n##\n'
        cmd += 'nlines=0\n'
        cmd += 'for chrom in {1..22}; do\n\n' ## magic number
        for fn in l_fn:
            if fn == 'prune/%s.$chrom.prune.in' %(out_prefix): continue
            cmd += 'if [ ! -s %s ]; then\n' %(fn)
            cmd += 'break\nfi\n\n'
        cmd += 'let nlines=$nlines+$('
        cmd += 'cat'
        for fn in l_fn:
            cmd += ' %s' %(fn)
        cmd += ' | wc -l'
        cmd += ')\n'
        cmd += '\ndone\n'
        l_cmds += [cmd]

        ## calculate expected lines
        cmd = 'n=$('
        if out_prefix == '%s.prehardy' %(self.o):
            cmd += 'cat %s.preIBD.frq' %(self.o)
        elif out_prefix == '%s.posthardy' %(self.o):
            cmd += 'cat %s.postIBD.frq' %(self.o)
        elif out_prefix == '%s.%s' %(self.o,self.fn1000g,):
            cmd += 'cat %s.%s.frq' %(self.o,self.fn1000g,)
        else:
            print(out_prefix)
            stop
        cmd += " | awk 'NR>1{if($5>=0.05) print}'"
        cmd += ' | wc -l'
        cmd += ')'
        l_cmds += [cmd]

        ##
        ## concatenanate lines if they are all present
        ## and do some additional stuff if they are...
        ##
        cmd = '##\n## concatenate lines if they are all present\n##\n'
        cmd += 'echo actual $nlines expected $n\n'
        cmd += 'if [ $nlines -ne $n ];then\nexit\nfi\n'
        ## exit after count if file out exists
        cmd += 'if [ -f %s.prune.in ];then\nexit\nfi\n' %(out_prefix)
        l_cmds += [cmd]

        cmd = '\n##\n## loop 3 (concatenate files)\n##\n'
        ## loop .prune.in files
        cmd += 'for chrom in {1..22}\ndo\n' ## magic number
        ## append to .prune.in file
        cmd += "cat prune/%s.$chrom.prune.in >> %s.prune.in\n" %(
            out_prefix,out_prefix,
            )
        cmd += 'echo %s.prune.in >> %s.touch\n' %(out_prefix,self.o,)
        ## append to .prune.out file
        cmd += "cat prune/%s.$chrom.prune.out >> %s.prune.out\n" %(
            out_prefix,out_prefix,
            )
        cmd += '\ndone\n\n'

        if self.bool_run_all == True:
            cmd += self.cmd_rerun('indep-pairwise',)
        l_cmds += [cmd]

        cmd = ''
        if (
            out_prefix == '%s.posthardy' %(self.o)
            or
            out_prefix == '%s.%s' %(self.o,self.fn1000g)
            ):
            cmd_EIGENSOFT = self.EIGENSOFT(
                in_prefix,out_prefix,bool_removal=False,)
            cmd += '\n##\n## EIGENSOFT\n##\n'+cmd_EIGENSOFT+'\n\n'
        l_cmds += [cmd]

        return l_cmds


    def genome_after(self,in_prefix,out_prefix,):

        l_cmds = []

        memMB = self.assign_memory(
            'genome',
            memMB_per_1000_samples=760,
            )

        ##
        ## continuation of nested loop and command from genome_before
        ##
        cmd = ''

        ## append to command
        if self.bool_run_all == True:
            cmd += ';'
            ## need to rerun at the end of each bsub
            cmd += self.cmd_rerun('genome')
        ## terminate command
        cmd += "\n'\n\n"
        ## evaluate command
        cmd_LSF = self.append_LSF(
            'genome',JOBID='genome.${i}.${j}',memMB=memMB,)
##        cmd += 'echo $cmd\n'
        cmd += '''cmd="%sbash -c '$cmd' $i $j"\n''' %(cmd_LSF)
        cmd += 'echo $cmd\n'
        cmd += 'eval $cmd\n'
        ## end loops
        cmd += '\ndone\ndone\n'
        ## exit before count if file out exists
        cmd += 'if [ -f %s.genome ];then\nexit\nfi\n' %(out_prefix)
        l_cmds += [cmd]

        ## count lines
        cmd = '\n##\n## loop 2a (count lines)\n##\n'
        cmd += 'nlines=0\n'
        cmd += 'break=0\n\n'
        cmd += 'for i in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
        cmd += 'if [ $break -eq 1 ]\nthen\nbreak\nfi\n'
        cmd += 'for j in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
        cmd += 'if [ $j -gt $i ]\nthen continue\nfi\n\n'
        cmd += 'if [ ! -s genome/%s.$i.$j.genome ]\nthen\n' %(out_prefix)
        cmd += 'break=1\nbreak\nfi\n\n'
        cmd += 'let nlines=$nlines+$('
        cmd += 'cat genome/%s.$i.$j.genome | wc -l' %(out_prefix)
        cmd += ')-1\n'
        cmd += '\ndone\ndone\n'
        l_cmds += [cmd]

        ##
        ## concatenanate lines if they are all present
        ## and do some additional stuff if they are...
        ##
        cmd = '##\n## concatenate lines if they are all present\n##\n'
        cmd += 'echo actual $nlines expected $(($n*($n-1)/2))\n'
        cmd += 'if [ $nlines -ne $(($n*($n-1)/2)) ];then\nexit\nfi\n'
        ## exit after count if file out exists
        cmd += 'if [ -f %s.genome ];then\nexit\nfi\n' %(out_prefix)
        l_cmds += [cmd]

        ## initiate .genome file with header
        s = '                   FID1                   IID1                   FID2                   IID2 RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO'
        cmd = '\n##\n## loop 2b (concatenate genome files)\n##\n'
        cmd += 'echo "%s" > %s.genome\n' %(s,out_prefix,)
        ## loop .genome files
        cmd += 'for i in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
        cmd += 'for j in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
        cmd += 'if [ $j -gt $i ]\nthen continue\nfi\n\n'
        ## append to .genome file
        cmd += "cat genome/%s.$i.$j.genome | awk 'NR>1' >> %s.genome\n" %(
            out_prefix,out_prefix,
            )
        cmd += '\ndone\ndone\n\n'
        cmd += 'echo %s.genome >> %s.touch\n' %(out_prefix,self.o,)

##        ## sort to always get the same result from IBD prune irrrespective of fam file sizes
##        ## for now using the same fam file size of 400 samples per fam file
##        cmd += 'sort -k2,2 -k3,3 %s.genome -o %s.genome\n\n' %(
##            out_prefix,out_prefix,)

        if out_prefix == '%s.prehardy' %(self.o):

            ##
            ## run IBD prune
            ##

            ## disallow duplicates or non-founders (e.g. IBD>0.90 or IBD>0.05)
            cmd += 'python %s/QC_IBD_prune.py' %(
                os.path.dirname(sys.argv[0]),
                )
            if self.args.fn_twins:
                cmd += ' --twins %s' %(self.args.fn_twins)
            ## do not include twins in sample exclusion list
            cmd += ' --pi_hat_max %.2f --genome %s.prehardy --imiss %s' %(
                float(self.args.threshold_pi_hat_max_postHWE),self.o,self.o,
                )
            cmd += ' --out %s.excltwins\n\n' %(self.o)

            ## disallow related samples for HWE (e.g. IBD>0.05; also twins)
            ## write minimal sample exclusion list
            f1 = float(self.args.threshold_pi_hat_max_preHWE)
            f2 = float(self.args.threshold_pi_hat_max_postHWE)
##            if f1 != f2 or self.args.fn_twins == None:
            if True:
                cmd += 'python %s/QC_IBD_prune.py' %(
                    os.path.dirname(sys.argv[0]),
                    )
                cmd += ' --pi_hat_max %.2f' %(
                    float(self.args.threshold_pi_hat_max_preHWE))
                cmd += ' --genome %s.prehardy --imiss %s --out %s\n\n' %(
                    self.o,self.o,self.o)
            cmd += 'echo %s.genome.%.2f.samples >> %s.touch\n' %(
                self.o,float(self.args.threshold_pi_hat_max_preHWE),self.o,)

            ##
            ## concatenate sample removal lists
            ##

            cmd += '\nif [ ! -f %s.excltwins.genome.%.2f.samples ]; then exit; fi\n' %(
                self.o,float(self.args.threshold_pi_hat_max_postHWE))

            ## after HWE
            ## concatenate sample exclusion lists (i.e. IBD>low threshold)
            cmd += 'cat %s.sampleQC.samples %s.excltwins.genome.%.2f.samples' %(
                self.o,self.o,float(self.args.threshold_pi_hat_max_postHWE),
                )
            cmd += ' > %s.sampleQC.IBD.samples\n' %(self.o,)
            cmd += 'echo %s.sampleQC.IBD.samples >> %s.touch\n\n' %(self.o,self.o,)

            ## before HWE
            cmd += 'cat %s.sampleQC.samples %s.genome.%.2f.samples' %(
                self.o,self.o,float(self.args.threshold_pi_hat_max_preHWE),
                )
            cmd += ' > %s.sampleQC.nonfounders.samples\n' %(self.o,)
            cmd += 'echo %s.sampleQC.nonfounders.samples >> %s.touch\n\n' %(self.o,self.o,)

        if self.bool_run_all == True:
            cmd += self.cmd_rerun('genome',)

##        cmd += self.mail('genome',out_prefix,'finished',)

        l_cmds += [cmd]

        return l_cmds


    def recodeA_after(self):

        l_cmds = []

        cmd = 'cat %s.raw' %(self.o)
        ## init awk
        cmd += " | awk '"
        ## loop rows
        cmd += ' NR>1 {cntNM=0;cntHET=0;'
        ## loop fields
        cmd += ' for(i=7;i<=NF;i++)'
        cmd += ' if($i==1) {cntHET++; cntNM++} else if($i!="NA") cntNM++;'
        ## loop rows continued
        cmd += ' het=cntHET/cntNM; print $1,$2,cntHET,cntNM,het;'
        cmd += ' }'
        ## term awk
        cmd += "'"
        cmd += ' > %s.het' %(self.o)
        l_cmds += [cmd]
        cmd = 'echo %s.het >> %s.touch\n' %(self.o,self.o,)
        l_cmds += [cmd]

        ## append header
        cmd = "sed -i '1i FID IID N(HET) N(NM) HET' %s.het" %(self.o)
##        cmd = 'echo "FID IID N(HET) N(NM) HET" > %s.het' %(bfile)
        l_cmds += [cmd]
        
        ## a) calculate heterozygosity mean and stddev
        l, mean,average, het_min, het_max = self.het2stddev(
            float(self.args.threshold_imiss),bool_execute=False,)
        l_cmds += l

        l_cmds += ['mean=$(echo $s | cut -d " " -f1)']
        l_cmds += ['stddev=$(echo $s | cut -d " " -f2)']

        l_cmds += ['echo mean $mean stddev $stddev']

        ## b) write list of samples outside mean+-3stddev
        cmd = 'awk -v mean=$mean -v stddev=$stddev '
        cmd += " 'NR>1 {if( "
        cmd += '$5<mean-%i*stddev || $5>mean+%i*stddev ' %(
            int(self.args.threshold_het_stddev),
            int(self.args.threshold_het_stddev),)
        cmd += ") print $1,$2}' "
        cmd += ' %s.het > %s.het.samples' %(self.o, self.o,)
        l_cmds += [cmd]
        cmd = 'echo %s.het.samples >> %s.touch\n' %(self.o,self.o,)
        l_cmds += [cmd]

##        ## remove huge .raw file to save disk space
        ## don't remove, otherwise it will be continuously regenerated
##        l_cmds += ['rm %s.raw' %(bfile)]

        return l_cmds


    def missing_after(self):

        l_cmds = []

        ##
        ## imiss.samples
        ##
        cmd = 'if [ -s %s.imiss -a ! -f %s.imiss.samples ]\n' %(self.o,self.o,)
        cmd += 'then\n'
        ## exclusion list
        cmd += "awk 'NR>1 {if($6>%f) print $1,$2}' %s.imiss " %(
            1-float(self.args.threshold_imiss),
            self.o,)
        cmd += ' > %s.imiss.samples\n' %(self.o)
        cmd += 'echo %s.imiss.samples >> %s.touch\n' %(self.o,self.o,)
        cmd += 'echo %s.imiss >> %s.touch\n' %(self.o,self.o,)
##        cmd += 'rm %s.lmiss\n' %(bfile)
        cmd += self.concatenate_sampleQC_remove_lists()
        cmd += 'fi\n'
        l_cmds += [cmd]

        ##
        ## lmiss.SNPs
        ##
        cmd = 'if [ -s %s.SNPQC.lmiss -a ! -f %s.lmiss.SNPs ]\n' %(self.o,self.o,)
        cmd += 'then\n'

        cmd += 'echo %s.SNPQC.lmiss >> %s.touch\n' %(self.o,self.o,)

        cmd += 'cat %s.SNPQC.lmiss' %(self.o)
        cmd += ' | '
        cmd += "awk 'NR>1 {if ((1-$5)<%.2f) print $2}' " %(
            float(self.args.threshold_lmiss),)
        cmd += ' > %s.lmiss.SNPs\n' %(self.o)
        cmd += 'echo %s.lmiss.SNPs >> %s.touch\n' %(self.o,self.o,)

        cmd += 'cat %s.lmiss.SNPs %s.ldregions.SNPs > %s.lmiss.exclLRLD.SNPs\n' %(
            self.o,self.o,self.o)
        cmd += 'echo %s.lmiss.exclLRLD.SNPs >> %s.touch\n' %(self.o,self.o,)
##        cmd += 'rm %s.SNPQC.imiss\n' %(bfile)
        cmd += 'fi\n'
        l_cmds += [cmd]

        ##
        ## lmiss.X.SNPs
        ##
        for sex in ['males','females',]:
            if sex == 'males' and self.bool_filter_females == True: continue
            cmd = 'if [ -s %s.X.%s.lmiss -a ! -f %s.X.%s.lmiss.SNPs ]\n' %(
                self.o,sex,self.o,sex,)
            cmd += 'then\n'
            cmd += 'cat %s.X.%s.lmiss' %(self.o,sex,)
            cmd += ' | '
            cmd += "awk 'NR>1 {if ((1-$5)<%.2f) print $2}' " %(
                float(self.args.threshold_lmiss),)
            cmd += ' > %s.X.%s.lmiss.SNPs\n' %(self.o,sex,)
            cmd += 'rm %s.X.%s.imiss\n' %(self.o,sex,)
            cmd += 'fi\n'
            l_cmds += [cmd]

        if self.bool_filter_females == True:
            cmd = 'if [ -f %s.X.females.lmiss.SNPs ]; then\n' %(
                self.o,)
            cmd += 'cat %s.X.females.lmiss.SNPs | sort -u > %s.X.lmiss.union.SNPs\n' %(
                self.o,self.o,)
        else:
            cmd = 'if [ -f %s.X.males.lmiss.SNPs -a -f %s.X.females.lmiss.SNPs ]; then\n' %(
                self.o,self.o,)
            cmd += 'cat %s.X.males.lmiss.SNPs %s.X.females.lmiss.SNPs | sort -u > %s.X.lmiss.union.SNPs\n' %(
                self.o,self.o,self.o,)
        cmd += 'echo %s.X.lmiss.union.SNPs >> %s.touch\n' %(self.o,self.o,)
        cmd += 'fi\n'
        l_cmds += [cmd]

        return l_cmds


    def concatenate_sampleQC_remove_lists(self):

        cmd = 'if [ '
        ## check file out
        cmd += '! -s %s.sampleQC.samples ' %(self.o)
        ## check file in
        cmd += '-a -s %s.imiss ' %(self.o)
        cmd += '-a -s %s.het ' %(self.o)
        if self.bool_no_sexcheck == False:
            cmd += '-a -s %s.sexcheck ' %(self.o)
        cmd += '-a -f %s.imiss.samples ' %(self.o)
        cmd += '-a -f %s.het.samples ' %(self.o)
        if self.bool_no_sexcheck == False:
            cmd += '-a -f %s.sexcheck.samples ' %(self.o)
        cmd += ']\n'
        cmd += 'then\n'
        ## concatenate sample files
        cmd += 'cat '
        cmd += '%s.imiss.samples ' %(self.o,)
        cmd += '%s.het.samples ' %(self.o,)
        if self.bool_no_sexcheck == False:
            cmd += '%s.sexcheck.samples ' %(self.o,)
        cmd += '| sort -u > %s.sampleQC.samples\n' %(self.o,)
        cmd += 'echo %s.sampleQC.samples >> %s.touch\n' %(self.o,self.o,)
        cmd += 'fi\n'

        return cmd
    

    def histogram_het(self,bool_with_stddev=True,):

        if os.path.isfile('%s.het.png' %(self.o)):
            return

        n_samples = int(os.popen('cat %s.het | wc -l' %(self.o)).read())-1
        cmd = 'cat %s.imiss.samples | wc -l' %(self.o)
        n_samples -= int(os.popen(cmd).read())
        cmd = 'cat %s.autosomes.SNPs | wc -l' %(self.o)
        n_SNPs = int(os.popen(cmd).read())

        x_step=0.0005

        if bool_with_stddev == True:

            fn_plot = '%s.het.joined' %(self.o)

            l_cmds,average,stddev,het_min,het_max = self.het2stddev(
                float(self.args.threshold_imiss),
                bool_execute=True,bool_remove=False,)

            l_arrows = [
                'set arrow from %f,0 to %f,graph(1) nohead lc 0\n' %(
                    average-self.args.threshold_het_stddev*stddev,
                    average-self.args.threshold_het_stddev*stddev,
                    ),
                'set arrow from %f,0 to %f,graph(1) nohead lc 0\n' %(
                    average+self.args.threshold_het_stddev*stddev,
                    average+self.args.threshold_het_stddev*stddev,
                    ),
                ]

            s_plot = 'plot [%f:%f]"%s" ' %(
                min(het_min,average-self.args.threshold_het_stddev*stddev,)-2*x_step,
                max(het_max,average+self.args.threshold_het_stddev*stddev,)+2*x_step,
                fn_plot,
                )
            s_plot += 'u (hist($5,width)):(1.0) smooth freq w boxes '
            s_plot += 'lc rgb"red" notitle\n'

            s_average = '%.4f' %(average)
            s_stddev = '%.4f' %(stddev)

            if int(os.popen('cat %s | wc -l' %(fn_plot)).read()) != n_samples:
                print(n_samples)
                print(int(os.popen('cat %s | wc -l' %(fn_plot)).read()))
                stop

            pass

        else:

            fn_plot = '%s.het' %(self.o)

            s_plot = 'plot "%s" ' %(fn_plot)
            s_plot += 'u (hist($5,width)):(1.0) smooth freq w boxes '
            s_plot += 'lc rgb"red" notitle\n'
            l_arrows = None

            s_average = 'N/A'
            s_stddev = 'N/A'

        if self.bool_verbose == True:
            bool_timestamp = True

        self.gnuplot_histogram(
            '%s.het' %(self.o),
##            x_min=het_min,x_max=het_max,#tic_step=0.05,
            x_step = x_step,
            s_plot = s_plot,
            xlabel='heterozygosity',
            title='%s (n_{samples}=%i, n_{SNPs}=%i, mean=%s, SD=%s)' %(
                self.o.replace('_','\\\\_'),n_samples,n_SNPs,
                s_average, s_stddev,
                ),
            lines_extra = l_arrows,
            color = 'red',
            bool_timestamp = True,
            )

        if bool_with_stddev == True:
            os.remove(fn_plot)

        return


    def histogram_genome(self):

##awk 'NR>1{if($10!="PI_HAT") print $0;}' $chip.genome > $chip.genome2

        if os.path.isfile('%s.prehardy.genome1.png' %(self.o)):
            return

        print('histogram PI_HAT', self.o)

        ## no data to plot
        if int(os.popen('head %s.prehardy.genome | wc -l' %(self.o,)).read()) == 1:
            return

##        cmd = 'cat %s.prehardy.genome | wc -l' %(bfile)
##        if self.verbose == True:
##            print cmd
##        n = int(os.popen(cmd).read())-1
##        ## http://simple.wikipedia.org/wiki/Quadratic_equation
##        n_samples1 = int(((--.5)+math.sqrt((-.5)**2-4*.5*-n))/(2*.5))

        n_samples2 = int(os.popen('cat %s.fam | wc -l' %(self.i)).read())
        n_samples2 -= int(os.popen('cat %s.sampleQC.samples | wc -l' %(self.o)).read())

##        if n_samples1 != n_samples2:
##            print n_samples1
##            print n_samples2
##            stop
##        else:
##            n_samples = n_samples1 = n_samples2

        n_samples = n_samples2

        n_SNPs = int(os.popen('cat %s.prehardy.prune.in | wc -l' %(self.o)).read())

        ##
        ## find highest pairwise IBD per sample
        ##
        d_IBD = {}
        ## populate dictionary with keys
        fd = open('%s.fam' %(self.i))
        lines = fd.readlines()
        fd.close()
        for line in lines:
            IID = line.split()[0]
            d_IBD[IID] = -1

        fd = open('%s.prehardy.genome' %(self.o),'r')
        ## skip header
        for line in fd: break

        for line in fd:
            l = line.split()
            IID1 = l[1]
            IID2 = l[3]
            PI_HAT = float(l[9])
##            if not IID1 in d_IBD.keys():
##                d_IBD[IID1] = 0
##            if not IID2 in d_IBD.keys():
##                d_IBD[IID2] = 0
            if PI_HAT > d_IBD[IID1]:
                d_IBD[IID1] = PI_HAT
            if PI_HAT > d_IBD[IID2]:
                d_IBD[IID2] = PI_HAT
        fd.close()

        ## remove IIDs that weren't in the .genome file
        for IID in list(d_IBD.keys()):
            if d_IBD[IID] == -1:
                del d_IBD[IID]

        l = ['%s\n' %(v) for v in list(d_IBD.values())]
        fd = open('%s.genome.max' %(self.o),'w')
        fd.writelines(l)
        fd.close()

        if self.bool_verbose == True:
            bool_timestamp = True

        for prefix,column,xlabel,x_min,x_max,x_step,prefix_out in [
            [
                '%s.prehardy.genome' %(self.o),
                '$10','PI CIRCUMFLEX',0,0.1,0.002,
                '%s.prehardy.genome1' %(self.o),
                ],
            ['%s.genome.max' %(self.o),'$1','PI CIRCUMFLEX MAX',0,1,0.01,'%s.genome.max' %(self.o),],
            [
                '%s.prehardy.genome' %(self.o),
                '$10','PI CIRCUMFLEX',0.1,1,0.002,
                '%s.prehardy.genome2' %(self.o),
                ],
            ]:

            self.gnuplot_histogram(
                prefix,
                x_min=x_min,
                x_step=x_step,
                x_max=x_max,
                column=column,
                xlabel=xlabel,
                title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                    self.o.replace('_','\\\\_'),n_samples,n_SNPs,
                    ),
                color = 'green',
                prefix_out=prefix_out,
                bool_timestamp = True,
                )

        os.remove('%s.genome.max' %(self.o))

        return


    def histogram_frq(self):

        if os.path.isfile('%s.postIBD.frq.png' %(self.o)):
            return

        cmd = 'cat %s.fam | wc -l' %(self.i)
        n_samples = int(os.popen(cmd).read())
        cmd = 'cat %s.sampleQC.IBD.samples | wc -l' %(self.o,)
        n_samples -= int(os.popen(cmd).read())
        cmd = 'cat %s.postIBD.frq | wc -l' %(self.o)
        n_SNPs = int(os.popen(cmd).read())-1

        if n_samples < 200:
            x_step=.01
        else:
            x_step=.005

        if self.bool_verbose == True:
            bool_timestamp = True

        self.gnuplot_histogram(
            '%s.postIBD.frq' %(self.o),
            x_step=x_step,
            x_max=0.5,
            column='$5',
##            s_plot = 'plot "%s.frq" u (hist($5,width)):(1.0) smooth freq w boxes lc rgb"blue" notitle\n' %(
##                bfile,
##                ),
            xlabel='MAF',
            title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                self.o.replace('_','\\\\_'),n_samples,n_SNPs,
                ),
            color = 'blue',
            bool_timestamp = True,
            )

        return


    def histogram_lmiss(self):

        if os.path.isfile('%s.lmiss.png' %(self.o)):
            return

        n_samples = int(os.popen('cat %s.fam | wc -l' %(self.i)).read())
        n_samples -= int(os.popen('cat %s.sampleQC.samples | wc -l' %(self.o,)).read())
        n_SNPs = int(os.popen('cat %s.SNPQC.lmiss | wc -l' %(self.o)).read())
        n_SNPs -= 1

        if n_samples > 200:
            x_step = .0002
        elif n_samples < 20:
            x_step = .1
        else:
            x_step = .01

        if self.bool_verbose == True:
            bool_timestamp = True

        cmd = "awk 'NR>1{print 1-$5}' %s.lmiss" %(self.o)
        x_min = min(
            0.95,
            float(self.args.threshold_lmiss),
            min([float(s) for s in os.popen(cmd).readlines()]),)
    
        self.gnuplot_histogram(
            '%s.SNPQC.lmiss' %(self.o),
            prefix_out='%s.lmiss' %(self.o),
            column='(1-$5)',
            x_step=x_step,
            x_min=x_min,
            x_max=1,
            xlabel='call rate, SNPs',
            title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                self.o.replace('_','\\\\_'),n_samples,n_SNPs,
                ),
            color = 'purple',
            bool_timestamp = True,
            )

        return


    def histogram_imiss(self):

        if os.path.isfile('%s.imiss.png' %(self.o)):
            return

        n_samples = int(os.popen('cat %s.imiss | wc -l' %(self.o)).read())-1
        if n_samples < 200:
            x_step=.005
        else:
            x_step=.002

        n_SNPs = int(os.popen('cat %s.autosomes.SNPs | wc -l' %(self.o)).read())

        cmd = "awk 'NR>1{print 1-$6}' %s.imiss" %(self.o)
        x_min = min(
            0.95,
            float(self.args.threshold_imiss),
            min([float(s) for s in os.popen(cmd).readlines()]),)

        if self.bool_verbose == True:
            bool_timestamp = True
   
        self.gnuplot_histogram(
            '%s.imiss' %(self.o),
            column='(1-$6)',
##            x_step=.0001,
            x_step=x_step,
            x_max=1,
            x_min=x_min,
            xlabel='call rate, samples',
            title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                self.o.replace('_','\\\\_'),n_samples,n_SNPs,
                ),
            color = 'orange',
            bool_timestamp = True,
            bool_remove = False,
            )

        return


    def init(self):

        ##
        ## check file existence
        ##
        if not os.path.isfile('%s.bed' %(self.i)):
            print('bed not found:', self.i)
            sys.exit(0)

        ##
        ## check X chromosome genotype data existence
        ##
        if self.bool_no_sexcheck == False:
            cmd = "cat %s.bim | awk '{if($1==23) print}' | wc -l" %(self.i)
            if int(os.popen(cmd).read()) == 0:
                print('No X chromosome SNPs.')
                print('Use --no-sex-check to avoid this error message')
                print('and skip the sex check')
                sys.exit(0)
        else:
            for fn in ['%s.sexcheck' %(self.o),'%s.sexcheck.samples' %(self.o),]:
                if os.path.isfile(fn): continue
                self.execmd('touch %s' %(fn))

        ##
        ## check that FIDs and IIDs in fam file are unique
        ## otherwise this script needs to be rewritten
        ## in particular when fgrep is used
        ##
        cmd_FID = "cat %s.fam | awk '{print $1}' | sort | uniq -d" %(self.i)
        FID_dup = os.popen(cmd_FID).read().strip()
        cmd_IID = "cat %s.fam | awk '{print $2}' | sort | uniq -d" %(self.i)
        IID_dup = os.popen(cmd_IID).read().strip()
        if FID_dup != '' and IID_dup != '':
            print('duplicate FIDs or IIDs')
            print(FID_dup)
            print(IID_dup)
            sys.exit()

        ##
        ## create touch file *after* initial checks
        ##
        if not os.path.isfile('%s.touch' %(self.o)):
            s = ''
            for extension in ['bed','bim','fam',]:
                s += '%s.%s\n%s.%s\n' %(self.i,extension,self.fp1000g,extension)
            fd = open('%s.touch' %(self.o),'w')
            fd.write(s)
            fd.close()

        ##
        ## create dirs
        ##
        for dn in ['stdout','stderr','fam','genome','prune',]:
            if not os.path.isdir(dn):
                os.mkdir(dn)
                pass
            continue

        ##
        ## write options in use to file
        ##
        l = []
        for k,v in list(vars(self.args).items()):
            l += [[k,v,]]
            continue
        l.sort()
        s = ''
        for k,v in l:
            s += '%-20s\t%s\n' %(k,v,)
            continue
        fd = open('%s.arguments' %(self.o),'w')
        fd.write(s)
        fd.close()

        ##
        ## write lists of autosomal SNPs and X chromosome SNPs
        ## the latter is used for --check-sex
        ## subtract concatenated SNP exclusion list
        ##
        self.write_list_of_autosomal_and_X_SNPs()

        ##
        ## convert range exclusion list to SNP exclusion list
        ## to allow --exclude and --extract to be combined
        ## which they cannot be, when --range is submitted for one of them
        ##
        self.convert_range_exclusion_to_SNP_exclusion()

        ## skip if output exists...
        self.d_out_suffix = {
            'missing':['imiss','lmiss',],
            'het':['het'],
            'cluster':['mds','cluster0','cluster1','cluster2','cluster3',],
            'freq':['frq'],
            'indep-pairwise':['prune.in','prune.out',],
            'check-sex':['sexcheck',],
            'genome':['genome',],
            'bmerge':['bed','fam','bim',],
            'hardy':['hwe',],
            'make-bed':['bed','fam','bim',],
            'recodeA':['raw',],
            }

        return


    def convert_range_exclusion_to_SNP_exclusion(self):

        fn_out = '%s.ldregions.SNPs' %(self.o)

        if os.path.isfile(fn_out):
            return

        fd = open(self.cluster_exclude,'r')
        lines = fd.readlines()
        fd.close()
        d_ranges = {}
        for line in lines:
            if line.strip() == '':
                continue
            l = line.split()
            chrom = l[0]
            i1 = int(l[1])
            i2 = int(l[2])
            if not chrom in list(d_ranges.keys()):
                d_ranges[chrom] = []
            d_ranges[chrom] += [[i1,i2,]]

        fd = open('%s.bim' %(self.i),'r')
        lines = []
        for line in fd:
            ## split line
            l = line.split()
            ## parse chromosome
            chrom = l[0]
            if not chrom in list(d_ranges.keys()):
                continue
            ## parse position
            pos = int(l[3])
            ## determine if within range
            bool_within_range = False
            for pos1,pos2 in d_ranges[chrom]:
                if pos > pos1 and pos < pos2:
                    bool_within_range = True
                    break
                continue
            if bool_within_range == False:
                continue
            SNP = l[1]
            lines += ['%s\n' %(SNP)]
        fd.close()

        fd = open(fn_out,'w')
        fd.writelines(lines)
        fd.close()

        return


    def execmd(self,cmd,):

        l_indexes = []
        l_cmd = cmd.split()

        if '%' in cmd and ".hwe | awk '{if(NR" not in cmd:
            print(cmd)
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
                if '<' not in cmd:
                    print(cmd)
                    print(inspect.stack()[1][3])
                    print('does not exist:', cmd.split()[index])
                    sys.exit(0)

        if self.bool_verbose == True:
            print()
            print(inspect.stack()[1][3])
            print(cmd)
        os.system(cmd)

        return


    def write_list_of_autosomal_and_X_SNPs(self):

        '''somehow I can never remember what this function does.
maybe I should rename it.'''

        for out_prefix, condition in [
            ['X_and_XY','$1==23||$1==25',],
            ['X_not_XY','$1==23',],
            ['autosomes','$1>=1 && $1<=22',], ## magic number
            ]:

            fn_out = '%s.%s.SNPs' %(self.o,out_prefix,)

            ## continue if file already exists
            if os.path.isfile(fn_out):
                continue

            ## generate inclusion list
            cmd = "cat %s.bim | awk '{if(%s) print $2}'  > %s" %(
                self.i,condition,fn_out,
                )
            self.execmd(cmd)

            fd = open('%s.touch' %(self.o),'a')
            fd.write('%s\n' %(fn_out))
            fd.close()

        return


    def parse_arguments(self,):

        parser = argparse.ArgumentParser()

        parser.add_argument(
            '-p', '--project', dest='project',
            metavar='STR',
            )

        parser.add_argument(
            '--threshold_pi_hat_max_postHWE', '--pi_hat_max_postHWE',
            dest='threshold_pi_hat_max_postHWE',
            help='Use to exclude non-founders (e.g. 0.05) or duplicates (e.g. 0.90)',
            metavar='FLOAT',default=0.90, type=float,
            )

        parser.add_argument(
            '--threshold_pi_hat_max_preHWE','--pi_hat_max_preHWE',
            dest='threshold_pi_hat_max_preHWE',
            help='The HWE step is dependent on the chosen PI CIRCUMFLEX MAX. Use to exclude non-founders before HWE check (e.g. 0.05)',
            metavar='FLOAT',default=0.05, type=float,
            )

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune
        parser.add_argument(
            '--indepWindow','--threshold_indepWindow',
            dest='threshold_indepWindow',
            help='The SNP pruning step (indep-pairwise) is dependent on the chosen value.',
            metavar='INT',default=50, type=int,
            )

        parser.add_argument(
            '--indepShift','--threshold_indepShift',
            dest='threshold_indepShift',
            help='The SNP pruning step (indep-pairwise) is dependent on the chosen value.',
            metavar='INT',default=5, type=int,
            )

        parser.add_argument(
            '--Rsquared','--threshold_indepRsquared',
            dest='threshold_indepRsquared',
            help='The SNP pruning step (indep-pairwise) is dependent on the chosen value.',
            metavar='FLOAT',default=.2, type=float,
            )

        parser.add_argument(
            '--threshold_het_stddev',
            dest='threshold_het_stddev',
            help='Used to determine, which samples fail the heterozygosity check.',
            metavar='INT',default=3, type=int,
            )

        ## http://pngu.mgh.harvard.edu/~purcell/plink/thresh.shtml
        parser.add_argument(
            '-i', '--mind', '--imiss',
            '--threshold_imiss',
            dest='threshold_imiss',
            help='The sample call rate calculation (--missing) is dependent on this value.',
            metavar='FLOAT',default=.97, type=float,
            )

        ## http://pngu.mgh.harvard.edu/~purcell/plink/thresh.shtml
        parser.add_argument(
            '-l', '--geno', '--lmiss',
            '--threshold_lmiss', '--threshold_lmiss',
            '--threshold_callrate_SNPs', '--threshold_callrate_SNP', '--threshold_missing_SNP',
            dest='threshold_lmiss',
            help='The SNP call rate calculation (--missing) is dependent on this value.',
            metavar='FLOAT',default=.97, type=float,
            )

        ## http://pngu.mgh.harvard.edu/~purcell/plink/thresh.shtml
        parser.add_argument(
            '--maf','--threshold_indepMAF',
            dest='threshold_indepMAF',
            help='Allele frequency threshold',
            metavar='FLOAT',default=.05, type=float,
            )

        ## http://pngu.mgh.harvard.edu/~purcell/plink/thresh.shtml
        parser.add_argument(
            '--hwe','--hwe_min', '--threshold_hwe_min', '--threshold_hwe',
            dest='threshold_hwe_min',
            help='HWe threshold used to determine failing SNPs after the HWE step.',
            metavar='FLOAT',default=10**-8, type=float,
            )

        parser.add_argument(
            '-b', '--bfile',
            dest='bfile',
            help='Specify an input file.',
            metavar='FILE', required=True,
            )
        
        parser.add_argument(
            '--fn_options', '-o', '--options',
            dest='fn_options',
            help='Specify an options input file.',
            metavar='FILE',
            )

        s_help = 'File used to label PCA/MDS plots;'
        s_help += 'column 1 = sample/individual ID, column2 = population/tribe'
        parser.add_argument(
            '--fn_pops','--pops','--dic',
            dest='fn_pops',
            help=s_help,
            metavar='FILE',default=None,
            )

        parser.add_argument(
            "-v", '--bool_verbose', '--verbose', action="store_true",
            dest="bool_verbose", default=True)
        parser.add_argument(
            "-q", '--silent', '--quiet', action="store_false",
            dest="bool_verbose")

        parser.add_argument(
            '--bool_filter_females', '--filter-females',
            action="store_true", dest="bool_filter_females", default=False)

        parser.add_argument(
            '--bool_no_sexcheck', '--no-sex-check',
            action="store_true", dest="bool_no_sexcheck", default=False)

        parser.add_argument(
            '--fn_twins', '--twins', metavar="FILE", type=str)

        parser.add_argument('--path_plink','--plink',required=True)

        ## parse arguments
        namespace_args = parser.parse_args()

        if namespace_args.bfile is None and namespace_args.arguments is None:
            parser.error('--bfile or --options not specified')
        if namespace_args.project is None and namespace_args.arguments is None:
            parser.error('--project or --options not specified')

        return namespace_args


    def __init__(self,):

        ##
        ## options
        ##

        args = self.parse_arguments()

        ##
        ## options file
        ##
        if args.fn_options not in [None,'None',False,'False',]:
            if not os.path.isfile(args.fn_options):
                print('does not exist:', args.fn_options)
                sys.exit(0)
            fd = open(args.arguments,'r')
            lines = fd.readlines()
            fd.close()
            d = {}
            for line in lines:
                l = line.strip().split()
                k = l[0]
                v = l[1]
                d[k] = v
                continue
        ##
        ## command line options
        ##
        else:
            d = vars(args)

        for k,v in list(d.items()):
            if k[:5] == 'bool_':
                if v in ['False',False,0,'0','no','n','false',]:
                    v = False
                elif v in ['True',True,1,'1','yes','y','true',]:
                    v = True
                else:
                    raise Exception('Unexpected value for %s: %s' %(k,v))
            setattr(self,k,v)
            setattr(args,k,v)

        self.verbose = args.bool_verbose

        self.i = self.bfile
        self.o = os.path.basename(self.bfile)

        ##
        ## other options
        ##

        ## length of lists used for pairwise IBD estimation run in parallel
        self.len_lists = 400

        self.bool_run_all = False
        if '--run-all' in sys.argv:
            self.bool_run_all = True
        self.bool_run_all = True
##        self.bool_run_all = False

        self.bool_sequential = False

        ##
        ## paths ("magic numbers")
        ##

        self.dn1000g = '/lustre/scratch113/projects/uganda_gwas/users/tc9/QC'
        self.fn1000g = '1000G_True' ## quad, build37, excl "chip effect" SNPs
        self.dn1000g = '/lustre/scratch113/projects/agv/AGV_datasets/global'
        self.fn1000g = 'globalagv_1.4M'

        self.dn1000g = '/lustre/scratch113/projects/uganda_gwas/users/tc9/PCA'
        self.fn1000g = 'globalagv_1.4M_exclug'
        self.dn1000g = '/lustre/scratch113/projects/uganda_gwas/users/tc9/PCA'
        self.fn1000g = 'Africaagv_1.6M_without_Uganda'

        self.dn1000g = '/lustre/scratch113/projects/uganda_gwas/users/tc9/PCA'
        self.fn1000g = 'globalagv_1.4M_excl_uganda_gwas_agv_intersection'
        self.dn1000g = '/lustre/scratch113/projects/uganda_gwas/users/tc9/PCA'
        self.fn1000g = 'Africaagv_1.6M_excl_uganda_gwas_agv_intersection'

        self.fp1000g = os.path.join(self.dn1000g,self.fn1000g,)

        ## dg11 2012oct09 06:32
        ## build 37
        self.cluster_exclude = '/nfs/t149_influenza_exomes/working/analysis2/ldregions.txt'

        ## http://genetics.med.harvard.edu/reich/Reich_Lab/Software.html
        ## http://computing.bio.cam.ac.uk/local/doc/eigenstrat.txt        
        ##Estimated running time of the smartpca program is 
        ##  2.5e-12 * nSNP * NSAMPLES^2 hours            if not removing outliers.
        ##  2.5e-12 * nSNP * NSAMPLES^2 hours * (1+m)    if m outlier removal iterations.
        ##Thus, under the default of up to 5 outlier removal iterations, running time is 
        ##  up to 1.5e-11 * nSNP * NSAMPLES^2 hours.
        self.eigensoft = '/nfs/team149/Software/EIG4.2/bin/smartpca'

        for file in (self.eigensoft,self.path_plink,self.cluster_exclude):
            if not os.path.isfile(file):
                print('does not exist:', file)
                sys.exit()

        self.args = args

        return

if __name__ == '__main__':
    instance_main = main()
    instance_main.main()
