#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

## built-ins
import os,sys,time,math,inspect,optparse
## add-ons
import numpy
import pwd
## http://docs.scipy.org/doc/scipy/reference/stats.html
from scipy import stats
## other
sys.path.append('/nfs/users/nfs_t/tc9/github/sandbox')
import gnuplot

## todo 2013-01-18: add option to do --options bfile.options

class main:

    ## Essential document describing the order of operations in PLINK:
    ## http://pngu.mgh.harvard.edu/~purcell/plink/flow.shtml

    '''This script assumes that FID and IID are always identical'''

    def main(self,bfile,):

        self.init(bfile,)

        self.check_logs() ## tmp!!!

        if self.verbose == True: print '############ execute ############'
        self.plink_execution(bfile,)

        if self.verbose == True: print '############ plot ############'
        self.plink_plots(bfile,)

        if self.verbose == True: print '############ tabulate ############'
        self.plink_tables(bfile,)

        return


    def check_logs(self,):

        print 'checking logs - temporary function!!!'

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
                print s
                print fn
####                os.remove(fn)
                stop_error
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
                print s
                print 'stderr/%s' %(fn)
####                os.remove('stderr/%s' %(fn))
                stop_error
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
                print s
                print fn
####                os.remove('stdout/%s' %(fn))
                stop_error

        print 'checked logs'

        return


    def scatter_hwe(self,bfile,):

        if os.path.isfile('%s.qq.hwe.png' %(bfile)):
            return

        n_samples = int(os.popen('cat %s.fam | wc -l' %(bfile)).read())
        cmd = 'cat %s.genome.prehardy.EIGENSOFT.samples | wc -l' %(bfile)
        n_samples -= int(os.popen(cmd).read())

        n_SNPs = (int(os.popen('cat %s.hwe | wc -l' %(bfile)).read())-1)/3

        cmd = 'cat %s.hwe' %(bfile)
        cmd += ' | '
        cmd += "awk 'NR>1{"
        cmd += 'if ($9!="NA") {logp=-log($9)/log(10); print logp} '
        cmd += "}'"
        cmd += ' | sort -r -n '
        ## http://gettinggeneticsdone.blogspot.co.uk/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html
        ## Thanks to Daniel Shriner at NHGRI for providing this code for creating expected and observed values
        cmd += " | awk '{print $1,-log(NR/%i)/log(10)}' " %(n_SNPs,)
        cmd += '> %s.qq.hwe.dat' %(bfile)
        self.execmd(cmd)

        gnuplot.scatter_plot_2d(
            '%s.qq.hwe' %(bfile),
##            s_plot = '"< paste %s.het %s.imiss" u (($5-$3)/$5):(1-$12)' %(bfile,bfile,),
##            line_plot = line_plot,
            column1 = 1, column2 = 2,
            xlabel = '-log(p_{HWE,observed})',
            ylabel = '-log(p_{HWE,expected})',
##            xmin=0,xmax=8,
            title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                bfile.replace('_','\\\\_'),n_samples,n_SNPs,
                ),
            )

        ##
        ## 2nd plot...
        ##
        line_plot = 'plot '
        for s_ineq, suffix,color in [
            ['<','out',1,],
            ['>','in',0,],
            ]:
            cmd = 'cat %s.hwe' %(bfile)
            cmd += " | awk '"
            cmd += 'NR>1{if ($9!="NA" && $9 %s %.1e) print}' %(
                s_ineq,self.hwe_min)
            cmd += "'"
            cmd += ' > %s.hwe.%s' %(bfile,suffix)
            self.execmd(cmd)
            line_plot += '"%s.hwe.%s" u 7:8 ps 1 pt 7 lc %i t "",' %(
                bfile,suffix,color,)
        line_plot = line_plot[:-1]+'\n'

        gnuplot.scatter_plot_2d(
            '%s.hwe.scatter' %(bfile),
            column1 = '$7', column2 = '$8',
            line_plot = line_plot,
            xlabel = 'O(HET)',
            ylabel = 'E(HET)',
##            xmin=0,xmax=8,
            title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                bfile.replace('_','\\\\_'),n_samples,n_SNPs,
                ),
            prefix_out = '%s.hwe.scatter' %(bfile),
            )

        os.remove('%s.hwe.in' %(bfile))
        os.remove('%s.hwe.out' %(bfile))

        return


    def histogram_hwe(self,bfile,):

        if os.path.isfile('%s.hwe.png' %(bfile)):
            return
        
##        cmd = '''awk 'NR>1{if ($9!="NA") print}' %s.hwe > %s.hwe.dat''' %(
##            bfile,bfile,
##            )
        cmd = 'cat %s.hwe' %(bfile)
        cmd += ' | '
        cmd += "awk 'NR>1{"
        cmd += 'if ($9!="NA") {logp=-log($9)/log(10); if(logp<=10) print logp}'
        cmd += "}'"
        cmd += '> %s.hwe.dat' %(bfile)
        self.execmd(cmd)

        n_samples = int(os.popen('cat %s.fam | wc -l' %(bfile)).read())
        cmd = 'cat %s.genome.prehardy.EIGENSOFT.samples | wc -l' %(bfile)
        n_samples -= int(os.popen(cmd).read())

        n_SNPs = (int(os.popen('cat %s.hwe | wc -l' %(bfile)).read())-1)/3

        gnuplot.histogram2(
##                '%s.hwe' %(bfile),
            '%s.hwe.dat' %(bfile),
            prefix_out = '%s.hwe' %(bfile),
##            x_step=0.01,
            x_step=0.1,
##            x_max=0.1,
##            column='$9',
            column='$1',
            xlabel='-log_{10}({/Helvetica-Italic p}_{HWE})',
            title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                bfile.replace('_','\\\\_'),n_samples,n_SNPs,
                ),
            color = 'yellow',
            )

        os.remove('%s.hwe.dat' %(bfile))

        return


##    def scatter_lmiss_frq(self,bfile,):
##
##        if os.path.isfile('%s.lmiss.frq.png' %(bfile)):
##            return
##
##        ##
##        ## join
##        ##
##        cmd = "sort -k2,2 %s.SNPQC.lmiss > %s.SNPQC.lmiss.sorted" %(bfile,bfile,)
##        self.execmd(cmd)
##        cmd = "sort -k2,2 %s.frq > %s.frq.sorted" %(bfile,bfile,)
##        self.execmd(cmd)
##        cmd = 'join -1 2 -2 2 %s.SNPQC.lmiss.sorted %s.frq.sorted > %s.lmiss.frq.joined' %(
##            bfile,bfile,bfile,
##            )
##        self.execmd(cmd)
##
##        os.remove('%s.SNPQC.lmiss.sorted' %(bfile))
##        os.remove('%s.frq.sorted' %(bfile))
##
##        line_plot = 'plot [.9:1][0:.5]'
##        line_plot += '"%s.lmiss.frq.joined" u (1-$5):10 lc 0 ps 1 pt 7 t ""' %(bfile,)
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
##        n_samples = int(os.popen('cat %s.lmiss.frq.joined | wc -l' %(bfile)).read())
##                            
##        cmd = '''cat %s.bim | awk '{if ($1<=22 && $1>=1) print}' | wc -l''' %(bfile)
##        cmd = '''cat %s.SNPQC.lmiss | wc -l''' %(bfile)
##        n_SNPs = int(os.popen(cmd).read())
##        ## subtract header
##        n_SNPs -= 1
##
##        gnuplot.scatter_plot_2d(
##            '%s.lmiss.frq.joined' %(bfile),
####            s_plot = '"< paste %s.lmiss %s.frq" u (1-$5):$10' %(bfile,bfile,),
####            s_plot = s_plot,
##            line_plot = line_plot,
##            column1 = '(1-$5)', column2 = 9,
##            xlabel = 'F MISS',
##            ylabel = 'MAF',
##            title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
##                bfile.replace('_','\\\\_'),n_samples,n_SNPs,
##                ),
####            ymax=.5,ymin=0,
######            xmin=0,xmax=1,
####            xmin=0.9,xmax=1,
##            lines_extra = l_arrows,
##            prefix_out='%s.lmiss.frq' %(bfile),
##            )
##
##        return


    def plink_tables(self,bfile,):

        for l_suffixes_in,function in [
            [['imiss'],self.table_imiss,],
            [['het'],self.table_het,],
            [['sexcheck'],self.table_sexcheck,],
            [['sampleQC.frq'],self.table_frq,],
            [['SNPQC.lmiss'],self.table_lmiss,],
            [['prehardy.genome'],self.table_genome,],
            [['hwe'],self.table_hwe,],
            ]:
            bool_continue = False
            for in_suffix in l_suffixes_in:
                fp_in = '%s.%s' %(bfile,in_suffix,)
                out_suffix = in_suffix
                fp_out = '%s.table' %(out_suffix)
                if not os.path.isfile(fp_in):
                    bool_continue = True
                    break
                if not os.path.getsize(fp_in) > 0:
                    bool_continue = True
                    break
                if not time.time()-os.path.getmtime(fp_in) > 5*60:
                    bool_continue = True
                    break
                if os.path.isfile(fp_out):
                    bool_continue = True
                    break
                continue
            if bool_continue == True:
                continue
            function(bfile,)

        ##
        ## combined
        ##
        print '\ncounts'
        for suffix in ['imiss','het','sexcheck',]:
            if not os.path.isfile('%s.%s.samples' %(bfile,suffix,)):
                continue
            cmd = 'cat %s.%s.samples | wc -l' %(bfile,suffix,)
            print suffix, os.popen(cmd).read(), bfile
        if (
            os.path.isfile('%s.imiss' %(bfile))
            and
            os.path.isfile('%s.het' %(bfile))
            and
            os.path.isfile('%s.sexcheck' %(bfile))
            ):
            set_combined = set()
            for suffix in ['imiss','het','sexcheck',]:
                l = os.popen('cat %s.%s.samples' %(bfile,suffix,)).readlines()
##                print len(l)
                set_combined |= set(l)
            print 'combined', len(set_combined)

        return


    def table_flip(self,bfile,):

        '''orphant function'''

##        cmd = 'sort %s.flip > %s.flip.sorted' %(bfile,bfile,)
##        self.execmd(cmd)
##
##        col_id = 2
##
##        cmd = 'cat %s.bim | sort -k%i,%i > %s.bim.sorted' %(bfile,col_id,col_id,bfile,)
##        self.execmd(cmd)
##        cmd = 'join -1 1 -2 %i ' %(col_id)
##        cmd += '%s.flip.sorted %s.bim.sorted' %(bfile,bfile,)
##        cmd += '> %s.bim.flip.joined' %(bfile)
##        self.execmd(cmd)
##
##        cmd = 'cat %s.flipped.bim | sort -k%i,%i > %s.flipped.bim.sorted' %(bfile,col_id,col_id,bfile,)
##        self.execmd(cmd)
##        cmd = 'join -1 1 -2 %i ' %(col_id)
##        cmd += '%s.flip.sorted %s.flipped.bim.sorted' %(bfile,bfile,)
##        cmd += '> %s.flipped.bim.flip.joined' %(bfile)
##        self.execmd(cmd)

##        cmd = 'cat %s.bim.flip.joined' %(bfile)
##        cmd = 'cat %s.flipped.bim.flip.joined' %(bfile)
##        cmd = 'cat %s.flipped.bim' %(bfile)
        cmd = 'cat %s.bim' %(bfile)
        cmd += " | awk '"
        ## BEGIN
        cmd += 'BEGIN{sumAC=0; sumAG=0; sumAT=0; sumCG=0; sumCT=0; sumGT=0; sumDI=0; sumother=0}'
        ## count
        cmd += ' '
        cmd += '{'
        cmd += 'if(($5=="A" && $6=="C") || ($5=="C" && $6=="A")) sumAC++; '
        cmd += 'if(($5=="A" && $6=="G") || ($5=="G" && $6=="A")) sumAG++; '
        cmd += 'if(($5=="A" && $6=="T") || ($5=="T" && $6=="A")) sumAT++; '
        cmd += 'if(($5=="C" && $6=="G") || ($5=="G" && $6=="C")) sumCG++; '
        cmd += 'if(($5=="C" && $6=="T") || ($5=="T" && $6=="C")) sumCT++; '
        cmd += 'if(($5=="G" && $6=="T") || ($5=="T" && $6=="G")) sumGT++; '
        cmd += 'if(($5=="D" && $6=="I") || ($5=="I" && $6=="D")) sumDI++; '
        cmd += 'if($5=="0" || $6=="0") sumother++ '
        cmd += '}'
        ## END
        cmd += ' '
        cmd += 'END{'
        cmd += 'print sumAC; print sumAG; print sumAT; '
        cmd += 'print sumCG; print sumCT; print sumGT; '
        cmd += 'print sumDI; print sumother}'
        cmd += "'"

##        os.remove('%s.bim.sorted' %(bfile))
##        os.remove('%s.flipped.bim.sorted' %(bfile))
##        os.remove('%s.flip.sorted' %(bfile))
##        os.remove('%s.bim.flip.joined' %(bfile))
##        os.remove('%s.flipped.bim.flip.joined' %(bfile))

        if self.verbose == True:
            print cmd
        print os.popen(cmd).read()

        return


    def table_frq(self,bfile,):

        fn_table = fn = 'frq.sampleQC.table'
        
        if os.path.isfile(fn_table): return

        cmd = '''awk 'NR>1{print $5}' %s.sampleQC.frq''' %(bfile,)
        if self.verbose == True:
            print cmd
        l_frq = [float(s) for s in os.popen(cmd).readlines()]

##        l_scores = [.0,.1,.2,.3,.4,.5,1.,2.,3.,4.,5.,10.,20.,30.,40.,50.,]
        l_scores = [.0,.1,.2,.5,1.,2.,5.,10.,20.,50.,]

        s = '%s' %(bfile)
        for score in l_scores:
            score /= 100
            try:
                n = int(len(l_frq)*stats.percentileofscore(l_frq, score, 'weak',))/100
            except:
                print score
                for frq in l_frq:
                    if type(frq) == str:
                        print frq
                n = int(len(l_frq)*stats.percentileofscore(l_frq, score, 'weak',))/100
                sys.exit(0)
            if self.verbose == True:
                print score, n
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

        os.system('sort %s -k1,1 -u -o %s' %(fn,fn))

        return


    def table_sexcheck(self,bfile,):

        fn_table = fn_out = 'sexcheck.table'

        if os.path.isfile(fn_table):
            return

        cmd = '''awk '{if($5=="PROBLEM") print $3,$4,$6}' %s.sexcheck''' %(
            bfile,
            )
        l_sexcheck = [s.split() for s in os.popen(cmd).readlines()]
        if self.verbose == True:
            print
            print 'sexcheck'
            print 'PEDSEX SNPSEX F'
            print l_sexcheck

##        ## opposite sex
##        cmd = 'cat %s.sexcheck' %(bfile,)
##        cmd += ''' | awk '{if($5=="PROBLEM" && $4!="0") print}' '''
##        cmd += '| wc -l'
##        n_opposite = int(os.popen(cmd).read())
##
##        ## unknown sex
##        cmd = 'cat %s.sexcheck' %(bfile,)
##        cmd += ''' | awk '{if($5=="PROBLEM" && $4=="0") print $1}' '''
##        cmd += '| wc -l'
##        n_unknown = int(os.popen(cmd).read())
##
##        n_total = n_opposite+n_unknown
##
##        if not os.path.isfile('sexcheck.table'):
##            fd = open('sexcheck.table','w')
##            fd.write('# total opposite_sex unknown_sex \n')
##            fd.close()
##            
##        fd = open('sexcheck.table','a')
##        fd.write('%s\t%i\t%i\t%i\n' %(bfile, n_total, n_opposite, n_unknown,))
##        fd.close()
##
##        os.system('sort -k1,1 -u -o sexcheck.table')

        cmd = 'cat %s.sexcheck' %(bfile,)
        cmd += ''' | awk '{if($5=="PROBLEM") print $2,$3,$4,$6}' '''
        lines = os.popen(cmd).readlines()
        lines = ['%s\t%s\n' %(bfile,'\t'.join(line.split())) for line in lines]

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


    def table_hwe(self,bfile,):

        l_range = range(4,8+1)

        for suffix in ['','.X.females',]:

            ## do loop over thresholds even though it's not optimal...
            s = '%s' %(bfile)
            for i in l_range:
                f = 10**-i
                cmd = 'grep ALL %s%s.hwe' %(bfile,suffix,)
                cmd += " | awk '{if ($9 < %.1e) print $2}' " %(f,)
                cmd += ' | wc -l'
                n = int(os.popen(cmd).read())
                print i, n
                s += '\t%i' %(n)
            s += '\n'

            fn_table = 'hwe%s.table' %(suffix)
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

            os.system('sort %s -k1,1 -u -o %s' %(fn_table,fn_table,))

        return


    def table_lmiss(self,bfile,):

        for suffix in ['SNPQC','X.males','X.females',]:

            if not os.path.isfile('%s.%s.lmiss' %(bfile,suffix,)):
                continue

            fn_table = 'lmiss.%s.table' %(suffix)

            if os.path.isfile(fn_table):
                continue

            s = ''
            cmd = '''awk 'NR>1{print 1-$5}' %s.%s.lmiss''' %(
                bfile,suffix,
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
                print 'cumfreq+extrapoints', l_cumfreq[0]
                print 'lowerreallimit', l_cumfreq[1]
                print 'binsize', l_cumfreq[2]
                print 'extrapoints', l_cumfreq[3]

            self.table_header_cumfreq(fn_table,lowerreallimit,binsize,)

            s = self.cumfreq2string(l_cumfreq,)
            line = '%s\t%s\n' %(bfile,s,)
            fd = open(fn_table,'a')
            fd.write(line)
            fd.close()

            os.system('sort %s -k1,1 -u -o %s' %(fn_table,fn_table,))

        return


    def table_header_cumfreq(self,fn_table,lowerreallimit,binsize,):

        if not os.path.isfile(fn_table):
            l = ['#']
            for i in xrange(int((1-lowerreallimit)/binsize)+1):
                l += ['%.3f' %(lowerreallimit+i*binsize)]
            s = '\t'.join(l)+'\n'
            fd = open(fn_table,'w')
            fd.write(s)
            fd.close()

        return
    

    def table_genome(self,bfile,):

        if (
            os.path.isfile('genome.table')
            and
            os.path.isfile('genome_outliers.table')
            ):
            return

        ##
        ## define pi hat range
        ##
        numbins = 4
        l_pi_hat_max = []
        for i in xrange(numbins):
            pi_hat_max = (i+1)*((.20-0.)/numbins)
            l_pi_hat_max += [pi_hat_max]
        l_pi_hat_max += [0.5,0.9,]

        ##
        ## check that output doesn't already exist
        ##
        bool_all = True
        for pi_hat_max in l_pi_hat_max:
            fn = '%s.genome.%.2f.samples' %(bfile,pi_hat_max,)
            if not os.path.isfile(fn):
                bool_all = False
                break
        if bool_all == True:
            return

        print '\n\ngenome'

        ##
        ## related samples at different thresholds
        ##

        n_samples = int(os.popen('cat %s.fam | wc -l' %(bfile)).read())
        n_samples -= int(os.popen('cat %s.sampleQC.samples | wc -l' %(bfile)).read())

        l = [bfile]
        l_header = ['#']
        for pi_hat_max in l_pi_hat_max:
            cmd = 'python %s/QC_IBD_prune.py ' %(os.path.dirname(sys.argv[0]))
            cmd += '--pi_hat_max %.02f --genome %s.prehardy --imiss %s --out %s\n\n' %(
                pi_hat_max, bfile, bfile, bfile,
                )
            self.execmd(cmd)
            cmd = 'cat %s.genome.%.2f.samples | wc -l' %(bfile,pi_hat_max,)
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

        os.system('sort %s -k1,1 -u -o %s' %(fn_table,fn_table,))

        ##
        ## highly related samples
        ##
        if self.project == 'uganda_gwas':
            cmd = '''awk 'NR>1{if($10>0.9) print $1,$3,$10}' %s.prehardy.genome''' %(bfile)
        elif self.project == 'agv':
            cmd = '''awk 'NR>1{if($10>0.1) print $1,$3,$10}' %s.prehardy.genome''' %(bfile)
        else:
            stop
        lines = os.popen(cmd).readlines()
        s = ''
        for line in lines:
            s += '%s\t%s\n' %(bfile,'\t'.join(line.split()))
        fn_out = 'genome_outliers.table'
        fd = open(fn_out,'a')
        fd.write(s)
        fd.close()
        os.system('sort -u %s | sort -k1,1 -o %s' %(fn_out,fn_out,))

        ##
        ## highly related samples
        ##
        fn = 'genome.high.table'
        cmd = '''awk 'NR>1{if($10>0.90) print $1,$3,$10}' '''
        cmd += '%s.prehardy.genome > %s.genome.high' %(bfile,bfile,)
        self.execmd(cmd)
        ## sort
        cmd = 'sort -k2,2 %s.imiss > %s.imiss.sorted' %(bfile,bfile,)
        self.execmd(cmd)
        ## sort and join 1
        cmd = 'sort -k1,1 %s.genome.high > %s.genome.high.sorted' %(bfile,bfile)
        self.execmd(cmd)
        cmd = 'join -1 2 -2 1 -o 0,2.3,1.6 '
        cmd += '%s.imiss.sorted %s.genome.high.sorted >> %s' %(
            bfile,bfile,fn,)
        self.execmd(cmd)
        ## sort and join 2
        cmd = 'sort -k2,2 %s.genome.high > %s.genome.high.sorted' %(bfile,bfile)
        self.execmd(cmd)
        cmd = 'join -1 2 -2 2 -o 0,2.3,1.6 '
        cmd += '%s.imiss.sorted %s.genome.high.sorted >> %s' %(
            bfile,bfile,fn,)
        self.execmd(cmd)
        ## get rid of duplicates
        cmd = 'sort -u %s -o %s' %(fn,fn)
        self.execmd(cmd)
        ## clean up
        os.remove('%s.genome.high' %(bfile))
        os.remove('%s.genome.high.sorted' %(bfile))
        os.remove('%s.imiss.sorted' %(bfile))

        return


    def table_imiss(self,bfile,):

        fn_table = 'imiss.table'
        if os.path.isfile(fn_table): return

        cmd = '''awk 'NR>1{print 1-$6}' %s.imiss''' %(bfile)
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
            print 'cumfreq', l_cumfreq[0]
            print 'lowerreallimit', l_cumfreq[1]
            print 'binsize', l_cumfreq[2]
            print 'extrapoints', l_cumfreq[3]

        self.table_header_cumfreq(fn_table,lowerreallimit,binsize,)

        s = self.cumfreq2string(l_cumfreq,)
        line = '%s\t%s\n' %(bfile,s,)
        fd = open(fn_table,'a')
        fd.write(line)
        fd.close()

        os.system('sort %s -k1,1 -u -o %s' %(fn_table,fn_table))

        return l_imiss


    def cumfreq2string(self,l_cumfreq,):

        s = '\t'.join(
            [str(l_cumfreq[3])]+[str(int(v)) for v in l_cumfreq[0]+l_cumfreq[3]]
            )

        return s


    def table_het(self,bfile,):

        fn_table = 'het.table'

        if os.path.isfile(fn_table):
            return

        l_cmds,average,stddev,het_min,het_max = self.het2stddev(
            bfile,self.threshold_imiss,bool_execute=True,bool_remove=False,)

        ## parse heterozygosities for stats.histogram
        cmd = 'cat %s.het.joined ' %(bfile,)
##        cmd += " | awk '{het=($5-$3)/$5; print het}'"
        cmd += " | awk '{print $5}'"
        l_het = [float(line) for line in os.popen(cmd).readlines()]
        os.remove('%s.het.joined' %(bfile))

        ## write stats to file
        fd = open('het.stats','w')
        fd.write('%s %f %f\n' %(bfile,average,3*stddev,))
        fd.close()

        ## histogram
        l_histogram = stats.histogram(
            l_het, numbins=2*self.threshold_het_stddev, defaultlimits=(
                average-self.threshold_het_stddev*stddev,
                average+self.threshold_het_stddev*stddev,
                ),
            )

        ## cumulated frequency...
        l_cumfreq = stats.cumfreq(
            l_het,numbins=1,defaultreallimits=(
                0,
                average+self.threshold_het_stddev*stddev,
                ),
            )

        histogram_extrapoints = l_histogram[3]
        cumfreq_extrapoints = l_cumfreq[3]

        if self.bool_verbose == True:
            print 'histogram', l_histogram[0]
            print 'less than lower', int(
##                sum(l_cumfreq[0])-sum(l_histogram[0][:])
                histogram_extrapoints-cumfreq_extrapoints
                )
            print 'greater than upper', int(
##                len(l_het)-sum(l_cumfreq[0])
                cumfreq_extrapoints
                )

        ## less than lower
        l = [str(int(histogram_extrapoints-cumfreq_extrapoints))]
        l += [str(int(v)) for v in l_histogram[0]]
        ## greater than upper
        l += [str((cumfreq_extrapoints))]
        s = '\t'.join(l)
        line = '%s\t%s\n' %(bfile,s,)
        fd = open(fn_table,'a')
        fd.write(line)
        fd.close()

        os.system('sort %s -k1,1 -u -o %s' %(fn_table,fn_table))

        ##
        ## table for appendix
        ##
##        ## 1) sort imiss
##        cmd = 'sort -k2,2 %s.imiss > %s.imiss.sorted' %(bfile,bfile,)
##        self.execmd(cmd)
##        ## 2) sort het
##        cmd = 'sort -k2,2 %s.het > %s.het.sorted' %(bfile,bfile,)
##        self.execmd(cmd)
##        ## 3) join
##        cmd = 'join -1 2 -2 2 -o 0,1.6,2.3,2.5 '
##        cmd + = '%s.imiss.sorted %s.het.sorted ' %(bfile,bfile,)
##        cmd += '> %s.imiss.het.joined' %(bfile)
##        self.execmd(cmd)
##        ## 4) fgrep -w -v -f imiss.samples
##        cmd = 'fgrep -w -v -f %s.imiss.1column.samples ' %(bfile,)
##        cmd += '%s.imiss.het.joined' %(bfile,)
##        cmd += '''| awk '{print $1,(($4-$3)/$4-%f)/%f}' ''' %(mean,stddev,)
##        cmd += '> %s.appendix
##        print cmd
##        stop
##        self.execmd()

##        fn_table = 'het.appendix.table'
##        cmd = 'fgrep -w -v -f %s.imiss.1column.samples ' %(bfile,)
##        cmd += '%s.het ' %(bfile,)
##        cmd += '''| awk 'BEGIN{OFS="\t"} NR>1{'''
##        cmd += 'if(($5-$3)/$5<%f' %(average-3*stddev)
##        cmd += '||'
##        cmd += '($5-$3)/$5>%f) ' %(average+3*stddev)
##        cmd += "print $2,($5-$3)/$5,(($5-$3)/$5-%f)/%f}' " %(average,stddev,)
##        cmd += '>> %s' %(fn_table)
##        self.execmd(cmd)
##
##        cmd = 'cat %s | sort -k1,1 -u -o %s' %(fn_table,fn_table)
##        self.execmd(cmd)
##
##        os.remove('%s.imiss.1column.samples' %(bfile))

        return l_het


    def xfrange(self, start, stop, step):
        while start < stop:
            yield start
            start += step


    def plink_plots(self,bfile,):

        for l_suffixes_in,function in [
            [['imiss'],self.histogram_imiss,],
            [['SNPQC.lmiss'],self.histogram_lmiss,],
            [['het'],self.histogram_het,],
            [['sampleQC.frq'],self.histogram_frq,],
            [['prehardy.genome'],self.histogram_genome,],
            [
                [
                    'imiss','het',
                    'sexcheck', ## red dots...
                    ],
                self.scatter_het_call,
                ],
            [['hwe',],self.histogram_hwe,],
            [['hwe',],self.scatter_hwe,],
##            [['SNPQC.lmiss','frq',],self.scatter_lmiss_frq,],
##            [['mds',],self.scatter_mds,],
            [['posthardy.mds',],self.scatter_mds_excl_1000g,],
            [['%s.mds' %(self.fn1000g),],self.scatter_mds_incl_1000g,],
            ]:
            bool_continue = False
            for in_suffix in l_suffixes_in:
                fp_in = '%s.%s' %(bfile,in_suffix,)
                if not os.path.isfile(fp_in):
                    bool_continue = True
                    break
                if not os.path.getsize(fp_in) > 0:
                    bool_continue = True
                    break
                if not time.time()-os.path.getmtime(fp_in) > 5*60:
                    bool_continue = True
                    break
            if bool_continue == True:
                continue
            function(bfile,)

        if (
            not os.path.isfile('venn3_%s.png' %(bfile))
            and
            os.path.isfile('%s.%s.samples' %(bfile,'imiss',))
            and
            os.path.isfile('%s.%s.samples' %(bfile,'het',))
            and
            os.path.isfile('%s.%s.samples' %(bfile,'sexcheck',))
            ):
            cmd = 'cat %s.het.samples' %(bfile)
            l_het = os.popen(cmd).read().split()
            cmd = 'cat %s.imiss.samples' %(bfile)
            l_imiss = os.popen(cmd).read().split()
            cmd = 'cat %s.sexcheck.samples' %(bfile)
            l_sexcheck = os.popen(cmd).read().split()
##            print len(l_het)
##            print len(l_imiss)
##            print len(l_sexcheck)
##            print 'imiss only', len(set(l_imiss)-set(l_sexcheck+l_het))
##            print 'het only', len(set(l_het)-set(l_imiss+l_sexcheck))
##            print 'sexcheck only', len(set(l_sexcheck)-set(l_imiss+l_het))
##            print 'imiss and het', len(set(set(l_imiss)&set(l_het))-set(l_sexcheck))
##            print 'het and sexcheck', len(set(set(l_sexcheck)&set(l_het))-set(l_imiss))
##            print 'imiss and sexcheck', 0, len(set(set(l_imiss)&set(l_sexcheck))-set(l_het))
##            print 'combined', len(set(l_sexcheck)|set(l_imiss)|set(l_het))
            if l_het != [] and l_imiss != [] and l_sexcheck != []:
                print 'venn3', len(l_het), len(l_imiss), len(l_sexcheck)
                try:
                    gnuplot.venn3(
                        l1=l_imiss,l2=l_het,l3=l_sexcheck,
                        suffix=bfile,
                        text1='call rate',text2='heterozygosity',text3='sex',
                        )
                except:
                    print 'need to fix venn3 function when no overlap...'

        return


    def scatter_mds_excl_1000g(self,bfile,):

        if os.path.isfile('%s.posthardy.mds.pc1.pc2.mds.png' %(bfile)):
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
                'IBO',
                ]
        else:
            l_tribes = [bfile,]

        fn = '%s.posthardy.mds' %(bfile)

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
                if FID in d_IID2tribe.keys():
                    prefix = d_IID2tribe[FID]
                else:
                    prefix = 'other'
                    print 'no tribe assigned', FID, os.popen('grep %s samples/*' %(FID)).read()
                if FID == 'APP5339441':
                    print tribe
                    stop
                fd = open('mds/%s.mds' %(prefix),'a')
                fd.write(line)
                fd.close()
                continue
        else:
            os.system('cp %s mds/%s.mds' %(fn,bfile))

        ##
        ## count samples and SNPs
        ##
        cmd = 'cat %s | wc -l' %(fn)
        n_samples = int(os.popen(cmd).read().strip())-1
        n_SNPs = int(os.popen('cat %s.posthardy.prune.in | wc -l' %(bfile)).read())

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

        for pc1 in xrange(1,5):
            for pc2 in xrange(1,5):
                if (
                    (pc1 != 1 or pc2 != 2)
                    and
                    (pc1 != 3 or pc2 != 4)
                    ):
                    continue
                if self.verbose == True:
                    print 'mds', pc1, pc2
                if pc2 <= pc1: continue
        
                line_plot = 'plot '
        ##        line_plot += '"%s.mds" u 4:5 lc 0 ps 2 pt 7 t ""\n' %(bfile)
                for i_tribe in xrange(len(l_tribes)):
                    tribe = l_tribes[i_tribe]
                    color = l_colors[i_tribe]
                    line_plot += '"mds/%s.mds" ' %(tribe)
                    line_plot += 'u %i:%i ' %(d_columns[pc1],d_columns[pc2],)
                    line_plot += 'lc rgb "#%s" ps 2 pt %i t "%s",' %(
                        "".join(map(chr, color)).encode('hex'),
                        l_pt[i_tribe%len(l_pt)],
                        tribe,
                        )
##                line_plot = line_plot[:-1]+'\n'

                ##
                ## add labels for MD outliers
                ##
                line_plot = self.add_MD_labels(
                    bfile,array_components,l_IIDs,pc1,pc2,line_plot,)

                ## finalize plot line
                line_plot = line_plot[:-1]+'\n'

                gnuplot.scatter_plot_2d(
                    '%s.mds' %(bfile),
        ##            s_plot = '"< paste %s.het %s.imiss" u (($5-$3)/$5):(1-$12)' %(bfile,bfile,),
        ##            s_plot = s_plot,
                    line_plot = line_plot,
##                    column1 = 4, column2 = 5,
                    xlabel = 'C%i' %(pc1),
                    ylabel = 'C%i' %(pc2),
                    title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                        bfile.replace('_','\\\\_'),n_samples,n_SNPs,
                        ),
                    prefix_out='%s.pc%i.pc%i.mds' %(fn,pc1,pc2,),
                    lines_extra=['set key out\n'],
##                    bool_remove=False,
                    )

        return


    def add_MD_labels(self,bfile,array_components,l_IIDs,pc1,pc2,line_plot,):

##        sys.path.append('/nfs/users/nfs_t/tc9/github/tc9/math')
##        import statistics
##        instance_tests = statistics.tests()

        u = array_components[:,pc1-1]
        v = array_components[:,pc2-1]
        ## http://docs.scipy.org/doc/numpy/reference/generated/numpy.cov.html
##                cov_pooled = np.cov(u,v,ddof=-1)
        cov_pooled = numpy.cov(u,v)
        cov_pooled_inverted = numpy.linalg.inv(cov_pooled)
        lines = ['#PC1\tPC2\tIID\tMD F p\n']
        for i in xrange(len(u)):
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
        fd = open('%s.mds.labels' %(bfile),'w')
        fd.writelines(lines)
        fd.close()
        
        line_plot += '"%s.mds.labels" u 1:2:3 w labels ' %(
            bfile,
            )
        line_plot += 'font "Helvetica,20" noenhanced t ""'

        return line_plot


    def scatter_mds_incl_1000g(self,bfile,):

        '''this function needs to be rewritten.
it's ugly and I will not understand it 1 year form now.'''

        if os.path.isfile('%s.%s.mds.1.2.mds.png' %(bfile,self.fn1000g,)):
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
                if s in d_sample2pop.keys():
                    if ethnicity == 'Unknown':
                        continue
                    print s, d_sample2pop[s], ethnicity
                    stop
                d_sample2pop[s] = ethnicity
                d_pops[ethnicity] += [s]
            fd = open('mds/%s_%s.mds' %(bfile,ethnicity),'w')
            fd.close()
            continue

        fn_mds = fn = '%s.%s.mds' %(bfile,self.fn1000g,)

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
            fd = open('mds/%s_%s.mds' %(bfile,prefix),'a')
            fd.write(line)
            fd.close()
##            print FID, os.popen('grep %s samples/*' %(FID)).read()
##        stop

        i = 0
        d_colors = {}
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
            'IBO',
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
            for i_pop in xrange(len(l_pops)):
                pop = l_pops[i_pop]
                fn_mds_pop = 'mds/%s_%s.mds' %(bfile,pop)
                if not os.path.isfile(fn_mds_pop):
                    print 'skipping', pop
                    continue
                if os.path.getsize(fn_mds_pop) == 0:
                    print 'skipping', pop
                    continue
                if pop in d_colors.keys():
                    color = d_colors[pop]
                    ps = 2
                    pt = 7
                else:
                    print pop
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
                if pop in d_colors.keys():
                    d_plot['uganda'] += s
                else:
                    d_plot['1000g'] += s
            line_plot += d_plot['uganda']
            line_plot += d_plot['1000g']
            line_plot = line_plot[:-1]+'\n'
            n_samples = int(os.popen('cat %s | wc -l' %(fn_mds)).read())
            n_samples -= 1
            cmd = 'cat %s.%s.prune.in | wc -l' %(bfile,self.fn1000g,)
            n_SNPs = int(os.popen(cmd).read())
            gnuplot.scatter_plot_2d(
                '%s.%s.mds' %(bfile,self.fn1000g,),
    ##            s_plot = '"< paste %s.het %s.imiss" u (($5-$3)/$5):(1-$12)' %(bfile,bfile,),
    ##            s_plot = s_plot,
                line_plot = line_plot,
                column1 = 4-1+c1, column2 = 4-1+c2,
                xlabel = 'C%i' %(c1),
                ylabel = 'C%i' %(c2),
                title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                    bfile.replace('_','\\\\_'),n_samples,n_SNPs,
                    ),
                prefix_out = '%s.%s.mds.%i.%i' %(bfile,self.fn1000g,c1,c2,),
                )

        return


    def scatter_het_call(self,bfile,bool_with_stddev=True,):

        '''this function is a mess. i was in a hurry, when I wrote it'''

        if os.path.isfile('%s.het.call.png' %(bfile)):
            return

        ## check that number of samples are equal
        ## which they should be if checks are run in *parallel*
        if (
            os.popen('cat %s.het | wc -l' %(bfile)).read()
            !=
            os.popen('cat %s.imiss | wc -l' %(bfile)).read()
            ):
            stop
        ## sort and join lists of samples in case they are not of equal length
        cmd = "sed '1d' %s.het | sort > %s.het.sorted" %(bfile,bfile,)
        self.execmd(cmd)
        cmd = "sed '1d' %s.imiss | sort > %s.imiss.sorted" %(bfile,bfile,)
        self.execmd(cmd)
        cmd = 'join %s.imiss.sorted %s.het.sorted > %s.imiss.het.joined' %(
            bfile,bfile,bfile,
            )
        self.execmd(cmd)

        os.remove('%s.het.sorted' %(bfile))
        os.remove('%s.imiss.sorted' %(bfile))

        if (
            os.popen('cat %s.het | wc -l' %(bfile)).read()
            !=
            os.popen('cat %s.sexcheck | wc -l' %(bfile)).read()
            ):
            stop

        ##
        ## join het and imiss (sexcheck samples only)
        ##
        cmd = "cat %s.sexcheck.samples | awk '{print $1}' " %(bfile)
        cmd += '> %s.sexcheck.samples.1column' %(bfile)
        self.execmd(cmd)
        ## use comm and join instead of grep... safer if IIDs are not unique...
        ## join1 (imiss)
        cmd = 'grep -f %s.sexcheck.samples.1column %s.imiss | ' %(bfile,bfile,)
        cmd += "awk '{print $1, 1-$6}' > %s_join1.txt" %(bfile,)
        self.execmd(cmd)
        ## join2 (het)
        cmd = 'grep -f %s.sexcheck.samples.1column %s.het' %(bfile,bfile,)
        cmd += ' | '
        cmd += "awk '{print $1,$5}' > %s_join2.txt" %(bfile,)
        self.execmd(cmd)
        ## join het and imiss
        cmd = "join %s_join1.txt %s_join2.txt > %s.sex.dat" %(
            bfile,bfile,bfile,
            )
        self.execmd(cmd)

        ## clean up
        os.remove('%s_join1.txt' %(bfile))
        os.remove('%s_join2.txt' %(bfile))
        os.remove('%s.sexcheck.samples.1column' %(bfile))

        ## opposite sex
        cmd = 'cat %s.sexcheck' %(bfile,)
        cmd += ''' | awk '{if($5=="PROBLEM" && $4!="0") print $1}' '''
        cmd += '> %s.opposite.txt' %(bfile)
        self.execmd(cmd)

        ## unknown sex
        cmd = 'cat %s.sexcheck' %(bfile,)
        cmd += ''' | awk '{if($5=="PROBLEM" && $4=="0") print $1}' '''
        cmd += '> %s.unknown.txt' %(bfile)
        self.execmd(cmd)

        ## use join instead of grep... safer in case of non-unique IDs...
        for suffix in ['opposite','unknown',]:
            ## join1 (imiss)
            cmd = "grep -f %s.%s.txt %s.imiss" %(bfile,suffix,bfile,)
            cmd += "| awk '{print $1, 1-$6}'"
            cmd += "> %s_join1.txt" %(bfile,)
            self.execmd(cmd)
            ## join2 (het)
            cmd = "grep -f %s.%s.txt %s.het" %(bfile,suffix,bfile,)
            cmd += "| awk '{print $1,$5}'"
            cmd += "> %s_join2.txt" %(bfile,)
            self.execmd(cmd)
            ## rm
            os.remove('%s.%s.txt' %(bfile,suffix,))
            ## join imiss het
            cmd = 'join %s_join1.txt %s_join2.txt > %s.sex.%s.dat' %(
                bfile,bfile,bfile,suffix,
                )
            self.execmd(cmd)
            os.remove('%s_join1.txt' %(bfile))
            os.remove('%s_join2.txt' %(bfile))

        ##
        ##
        ##
        n_samples = int(os.popen('cat %s.het | wc -l' %(bfile)).read())-1
        if bool_with_stddev == True:

            l_cmds,average,stddev,het_min,het_max = self.het2stddev(
                bfile,self.threshold_imiss,bool_execute=True)

        n_SNPs = int(os.popen('cat %s.autosomes.SNPs | wc -l' %(bfile)).read())

        ## horizontal lines
        l_arrows = [
##            'set arrow from graph(0),%f to graph(1),%f nohead lc 0\n' %(
            'set arrow from 0,%f to 1,%f nohead lc 0 lt 2 lw 2\n' %(
                self.threshold_imiss,self.threshold_imiss,
                ),
            'set arrow from 0,%f to 1,%f nohead lc 0 lt 2 lw 2\n' %(
                0.99,0.99,
                ),
            ]
        ## vertical lines
        if bool_with_stddev == True:
            l_arrows += ['set arrow from %f,graph(0) to %f,graph(1) ' %(
                average-3*stddev,average-3*stddev,
                )]
            l_arrows += ['nohead lc 0 lt 1 lw 3\n']
            l_arrows += ['set arrow from %f,graph(0) to %f,graph(1) ' %(
                average+3*stddev,average+3*stddev,
                )]
            l_arrows += ['nohead lc 0 lt 1 lw 3\n']

        if bool_with_stddev == True:
            d_colors = {
                .97:'green',.99:'blue',
                .00:'black',
                #.98:'orange',
                }
            l_thresholds = d_colors.keys()
            l_thresholds.sort()
            l_thresholds.reverse()
            for threshold in l_thresholds:

                (
                    l_cmds,average_w_exclusion,stddev_w_exclusion,het_min,het_max,
                    ) = self.het2stddev(bfile,threshold,bool_execute=True,)
                
                color = d_colors[threshold]
                l_arrows += [
                    'set arrow from %f,graph(0) to %f,graph(1) ' %(
                        average_w_exclusion-3*stddev_w_exclusion,
                        average_w_exclusion-3*stddev_w_exclusion,
                        ),
                    'nohead lc rgb"%s" lt 2 lw 3\n' %(color),
                    'set arrow from %f,graph(0) to %f,graph(1) ' %(
                        average_w_exclusion+3*stddev_w_exclusion,
                        average_w_exclusion+3*stddev_w_exclusion,
                        ),
                    'nohead lc rgb"%s" lt 2 lw 3\n' %(color),
                    ]
    ##            print threshold, len(l_het), average_w_exclusion, 3*stddev_w_exclusion

        if bool_with_stddev == True:
            line_plot = 'plot [%f:%f][:1]' %(
                min(het_min,average-3*stddev,)-0.001,
                max(het_max,average+3*stddev,)+0.001,
                )
        else:
            line_plot = 'plot [:][:1]'
        line_plot += '"%s.imiss.het.joined" ' %(bfile)
        line_plot += 'u ($10):(1-$6) lc 0 ps 2 pt 7 t ""'

        ## sexcheck dots and labels
        l_colors = ['red','orange',]
        for i in xrange(2):
            suffix = ['opposite','unknown',][i]
            color = l_colors[i]
            if not os.path.getsize('%s.sex.%s.dat' %(bfile,suffix,)) > 0:
                continue
            line_plot += ','
            line_plot += '"%s.sex.%s.dat" u 3:2:1 lc rgb"%s" ps 2 pt 7 t ""' %(
                bfile,suffix,color,
                ) ## het and sex points
            line_plot += ','
            line_plot += '"%s.sex.%s.dat" u 3:2:1 w labels ' %(
                bfile,suffix,
                ) ## het and sex labels
            line_plot += 'font "Helvetica,16" noenhanced lc 2 t ""'
            continue
        line_plot += '\n'

        gnuplot.scatter_plot_2d(
            '%s.het.call' %(bfile),
##            s_plot = '"< paste %s.het %s.imiss" u (($5-$3)/$5):(1-$12)' %(bfile,bfile,),
            line_plot = line_plot,
            column1 = '(1-$6)', column2 = 12,
            xlabel = 'heterozygosity',
            ylabel = 'call rate',
            title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                bfile.replace('_','\\\\_'),n_samples,n_SNPs,
                ),
            lines_extra = l_arrows,
            )

        os.remove('%s.sex.opposite.dat' %(bfile))
        os.remove('%s.sex.unknown.dat' %(bfile))
        os.remove('%s.sex.dat' %(bfile))
        os.remove('%s.imiss.het.joined' %(bfile))

        return


    def het2stddev(
        self,bfile,threshold_imiss,bool_execute=False,bool_remove=True,):

        l_cmds = []

        cmd = 'cat %s.imiss' %(bfile)
        cmd += " | awk 'NR>1 {if($6>%f) print $1,$2}'" %(1-threshold_imiss)
        cmd += ' | sort -k2,2'
        cmd += ' > %s.imiss.remove.sorted' %(bfile)
        l_cmds += [cmd]

        cmd = "cat %s.het | awk 'NR>1' | sort -k2,2 > %s.het.sorted" %(bfile,bfile,)
        l_cmds += [cmd]

        cmd = 'join -1 2 -2 2 -v2'
        cmd += ' %s.imiss.remove.sorted %s.het.sorted ' %(bfile,bfile,)
        cmd += ' > %s.het.joined' %(bfile)
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
        cmd += ' cat %s.het.joined' %(bfile)
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

        cmd = 'rm %s.het.sorted %s.imiss.remove.sorted' %(bfile,bfile,)
        l_cmds += [cmd]
        if bool_execute == True:
            self.execmd(cmd)

        if bool_remove == True:
            cmd = 'rm %s.het.joined' %(bfile,)
            l_cmds += [cmd]
            if bool_execute == True:
                self.execmd(cmd)

        return l_cmds, average, stddev, het_min, het_max

    

    def plink_execution(self,bfile,):

        l_plink_cmds = self.l_plink_cmds(bfile,)

        ##
        ## execute plink commands
        ##
        for plink_cmd_full in l_plink_cmds:

            plink_cmd = plink_cmd_full.split()[0].replace('--','')

            ##
            ## check that file input exists (LSF dependency does not work...)
            ##
            bool_continue_in, in_prefix, fn_in = self.check_input_existence(
                bfile,plink_cmd,plink_cmd_full,
                )
            if bool_continue_in == True:
                print 'in does not exist', plink_cmd, fn_in
                continue

            ##
            ## check that file output doesn't already exist
            ## and immediately create it if it doesn't
            ##
            bool_continue_out, out_prefix, fn_out = self.check_output_existence(
                bfile,plink_cmd,plink_cmd_full,
                )
            if bool_continue_out == True:
                print 'out exists', plink_cmd, fn_out
                continue

            ##
            ## initiate command
            ##
            cmd = ''
            if self.bool_verbose == True:
                cmd += self.mail(bfile,plink_cmd,out_prefix,'initiated',)

            ##
            ## before plink
            ##
            cmds_extra_before = self.extra_commands_before(
                bfile,plink_cmd,plink_cmd_full,
                in_prefix, out_prefix,
                )
            if cmds_extra_before != '':
                cmd += '\n%s\n' %(cmds_extra_before,)
            else:
                cmd += '\n'

            ##
            ## plink cmd
            ##
            cmd_plink = self.append_plink(bfile,plink_cmd_full,out_prefix,)
            cmd += cmd_plink

            ##
            ## after plink
            ##
            cmds_extra_after = self.extra_commands_after(
                bfile,plink_cmd,in_prefix,out_prefix,
                )
            if cmds_extra_after != '':
                cmd += '\n%s\n' %(cmds_extra_after)
            else:
                cmd += '\n'

            ##
            ## LSF settings
            ##
            cmd_LSF = self.append_LSF(bfile,plink_cmd,plink_cmd_full,)
            cmd_LSF += './%s_%s.sh' %(out_prefix,plink_cmd,)

            if (
                self.bool_run_all == True
                and
                ## only rerun after completing a step in the "genome loop"
                ## not simply after exiting the "genome loop"
                plink_cmd != 'genome'
                ):
                cmd += self.cmd_rerun(bfile,plink_cmd,)

            ##
            ## terminate command
            ##
            if self.bool_verbose == True:
                cmd += self.mail(bfile,plink_cmd,out_prefix,'finished',)

            ##
            ## write plink command and associated commands
            ## to individual shell scripts
            ##
            fd = open('%s_%s.sh' %(out_prefix,plink_cmd,),'w')
            fd.write(cmd)
            fd.close()
            ## make shell script executable
            os.system('chmod +x %s_%s.sh' %(out_prefix,plink_cmd,))

            if self.bool_verbose == True:
                print
                print plink_cmd
##                print cmd
##                print cmd_LSF

            if not '--remove' in cmd and plink_cmd not in [
                'missing','het','recodeA','check-sex',
                ]:
                print cmd
                stop1
            if '%' in cmd.replace('%02g',''):
                print cmd
                stop2_parenthesis_in_cmd
##            print cmd
##            print cmd_LSF
##            continue
##            stop3

            self.execmd(cmd_LSF)

        return


    def mail(self,bfile,plink_cmd,out_prefix,status,):

        msg = 'job $LSB_JOBID %s' %(status,)
        subject = 'job $LSB_JOBID %s, %s %s %s' %(
            status,bfile,plink_cmd,out_prefix,)

        cmd = ''
        address = '%s@sanger.ac.uk' %(pwd.getpwuid(os.getuid())[0],)
##        address += '%s@sanger.ac.uk\n' %(os.getlogin(),)
##        address = 'tommy.carstensen@gmail.com' ## tmp!!!
        cmd += '\n'
        cmd += 'echo "%s" ' %(msg)
        cmd += '| mail -s "%s" ' %(subject)
        cmd += '%s\n' %(address)

        return cmd


    def cmd_rerun(self,bfile,plink_cmd,nl='\n',):

        cmd = '\n'
##        cmd += 'if [ ! -s %s.SNPQC.bed ]; then\n' %(bfile)
        cmd += 'if [ ! -s %s.SNPQC.bed -a ! -s %s.posthardy.mds ]; then\n' %(
            bfile,bfile,)
        if plink_cmd == 'indep-pairwise':
            cmd += 'sleep 600;'
        else:
            cmd += 'sleep 300;'
        cmd += '/software/bin/python-2.7.3 %s/QC.py' %(
            os.path.dirname(sys.argv[0]),
            )
        ## snould add all variables from .options file if default values are not used...
        cmd += ' --project %s --bfile %s' %(self.project,bfile,)
        ## indep-pairwise
        cmd += ' --indepWindow %i' %(self.indepWindow)
        cmd += ' --indepShift %i' %(self.indepShift)
        cmd += ' --Rsquared %f' %(self.Rsquared)
        ## genome
        cmd += ' --pi_hat_max %f' %(self.pi_hat_max)
        ## het
        cmd += ' --threshold_het_stddev %i' %(self.threshold_het_stddev)
        ## freq
        cmd += ' --threshold_maf %f' %(self.maf_min)
        ## missing
        cmd += ' --threshold_imiss %f' %(self.threshold_imiss)
        cmd += ' --threshold_lmiss %f' %(self.threshold_lmiss)
        ## hardy
        cmd += ' --threshold_hwe_min %.1e' %(self.hwe_min)
        cmd += nl
        cmd += 'fi'
        cmd += '\n'

        return cmd


    def check_output_existence(self,bfile,plink_cmd,plink_cmd_full,):

        bool_continue = False

        if '--out' in plink_cmd_full:
            l = plink_cmd_full.split()
            if plink_cmd == 'genome':
                s = l[l.index('--out')+1]
                out_prefix = s[s.index('/')+1:-6]
            else:
                out_prefix = l[l.index('--out')+1]
        elif '--bfile' in plink_cmd_full:
            l = plink_cmd_full.split()
            out_prefix = l[l.index('--bfile')+1]
        else:
            out_prefix = bfile

        bool_file_out = False
        for out_suffix in self.d_out_suffix[plink_cmd]:
            fn_out = '%s.%s' %(out_prefix,out_suffix,)
            if os.path.isfile(fn_out):
                bool_file_out = True
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
        if bool_file_out == False:
            fn_log = '%s.log' %(out_prefix)
            if os.path.isfile(fn_log):
                if time.time()-os.path.getmtime(fn_log) < 5*60:
                    bool_file_out = True
                    fn_out = fn_log
                    pass
                pass
            ##
            ## create the output so another process doesn't think it's not created
            ##
            else:
                fd = open(fn_log,'w')
                fd.close()
                
        if bool_file_out == True:
            bool_continue = True

        return bool_continue, out_prefix, fn_out


    def check_input_existence(self,bfile,plink_cmd,plink_cmd_full,):

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
            bfile_prefix = bfile
            pass
        l_fp_in += ['%s.bed' %(bfile_prefix)]

        ## append extract, remove, exclude, keep
        for option in [
            'extract','remove','exclude','keep',
            'read-freq', 'read-genome',
            'bmerge',
            ]:
##            if plink_cmd == 'bmerge' and option == 'extract':
##                fp_in = '%s.posthardy.prune.in' %(bfile,)
##                l_fp_in += [fp_in]
##                pass
##            elif '--%s' %(option) in plink_cmd_full:
##                l = plink_cmd_full.split()
##                fp_in = l[l.index('--%s' %(option))+1]
##                l_fp_in += [fp_in]
##                pass
            if '--%s' %(option) in plink_cmd_full:
                l = plink_cmd_full.split()
                if option == 'bmerge':
                    fp_in = l[l.index('--%s' %(option))+2]
                else:
                    fp_in = l[l.index('--%s' %(option))+1]
                l_fp_in += [fp_in]
                pass
            continue

        ## check                
        bool_input = True
        for fp_in in l_fp_in:
            bool_input = self.check_file_in(fp_in)
            if bool_input == False:
                if self.bool_verbose == True:
                    print plink_cmd, fp_in, 'does not exist or was recently modified'
                    pass
                break
            continue
        if bool_input == False:
            bool_continue = True
            pass

        return bool_continue, bfile_prefix, fp_in


    def append_plink(self,bfile,plink_cmd_full,out_prefix,):

        ## http://pngu.mgh.harvard.edu/~purcell/plink/flow.shtml
        
        plink_cmd = plink_cmd_full.split()[0].replace('--','')

        ## initiate plink
        cmd = 'plink \\\n'

        ## --bfile
        if '--bfile' not in plink_cmd_full:
            cmd += '--bfile %s \\\n' %(bfile,)

        ##
        cmd += '%s \\\n' %(plink_cmd_full,)

        ## --extract SNPs
        if '--extract' in plink_cmd_full:
            pass
        elif plink_cmd == 'cluster':
            pass
        else:
            cmd += '--extract %s.autosomes.SNPs \\\n' %(bfile,)

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

        if plink_cmd == 'genome':
            ## backslash causes trouble when using bash -c 'command \ --option'
            cmd = cmd.replace('\\\n','')
            cmd = cmd.replace('$i','$0')
            cmd = cmd.replace('$j','$1')

        return cmd


    def append_LSF(
        self, bfile, plink_cmd, plink_cmd_full=None, JOBID=None, verbose=True,
        ):

        if plink_cmd_full:
            mem = self.assign_memory(bfile,plink_cmd,plink_cmd_full,)
        else:
            mem = 4
        
        cmd = ''
        cmd += 'bsub \\\n'
        cmd += "-M%i000000 -R'select[mem>%i000] rusage[mem=%i000]' \\\n" %(
            mem,mem,mem,
            )
        cmd += '-G %s \\\n' %(self.project)
        cmd += '-q normal \\\n'
        ## redirect stdout and stderr
        if verbose == True and self.bool_verbose == True:
            cmd += '-o stdout/%s.%s.out \\\n' %(bfile,plink_cmd,) ## tmp
            cmd += '-e stderr/%s.%s.err \\\n' %(bfile,plink_cmd,) ## tmp
        else:
            cmd += '-o /dev/null \\\n'
        ## specify JOBID
        if not JOBID:
            JOBID = '%s.%s' %(bfile,plink_cmd,)
##        cmd += '-J"%s" \\\n' %(JOBID,)
        cmd += "-J'%s' \\\n" %(JOBID,)
##            ##
##            ## wait for input to exist and be more than 5 minutes old
##            ## doesn't seem to be supported by LSF7 at the Sanger :(
##            ##
##            if plink_cmd in d_plink_dep_file.keys():
##                cmd += '-w "'
##                l = []
##                for in_suffix in d_plink_dep_file[plink_cmd]:
##                    fp_in = '%s.%s' %(bfile,in_suffix,)
##                    l += ['file(age(%s.%s) > 5M)' %(bfile,in_suffix,)]
##                s = ' && '.join(l)
##                cmd += s
##                cmd += '" \\'
##            ##
##            ## wait for preceding job to finish (done or exit)
##            ## should be changed to check that *all* job IDs finished...
##            ## this only works if all jobs submitted all at once...
##            ##
##            if plink_cmd in d_plink_dep_file.keys():
##                cmd += '-w "ended(%s.%s)" \\' %(bfile,d_plink_dep_done[plink_cmd][0])

##        cmd += '\\\n'

        return cmd


    def assign_memory(self,bfile,plink_cmd,plink_cmd_full,):

        ## approximate memory (MB) required per 1000 samples
        ## depends on number of SNPs as well...
        d_mem = {
            'check-sex':800,
            'cluster':900,
            'het':1100,
            'bmerge':1300,
            'freq':800,
            'indep-pairwise':800,
            'make-bed':800,
            'missing':800,
            'genome':800, ## memory requirements low if EIGENSOFT bed not to be written...
            }

        if '--bfile' not in plink_cmd_full:
            prefix_in = bfile
        else:
            l = plink_cmd_full.split()
            prefix_in = l[l.index('--bfile')+1]

        cmd = 'cat %s.fam | wc -l' %(bfile)
        n_samples = int(os.popen(cmd).read())
        ## additional samples if bmerge
        if plink_cmd == 'bmerge':
            if os.path.isfile(self.fp1000g):
                cmd = 'cat %s.fam | wc -l' %(self.fp1000g)
                if self.bool_verbose == True:
                    print cmd
                n_samples += int(os.popen(cmd).read())

        ## assign appropriate amount of memory (GB)
        if plink_cmd in d_mem.keys():
            mem = int(
                max(3,math.ceil((d_mem[plink_cmd]/1000.)*n_samples/1000.))
                )
        else:
            mem = int(max(3,math.ceil(0.8*n_samples/1000.)))
            pass

        return mem


    def l_plink_cmds(self,bfile,):

        l_plink_cmds = []

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune
        s_indep_pairwise = '--indep-pairwise %i %i %f --maf %f \\\n' %(
            self.indepWindow,
            self.indepShift,
            self.Rsquared,
            self.maf_min,
            )

        ## http://pngu.mgh.harvard.edu/~purcell/plink/strat.shtml#cluster
        s_cluster = '--cluster \\\n'
        s_cluster += '--mds-plot 10 \\\n' ## 10 components
        s_cluster += '--exclude %s.ldregions.SNPs \\\n' %(bfile,)

##        ##
##        ## sequential for Manj
##        ##
##        l_plink_cmds += ['--het --remove %s.imiss.samples --out %s.sequential' %(bfile,bfile,)]
##        l_plink_cmds += ['--check-sex --extract %s.X.SNPs --remove %s.imiss.samples --out %s.sequential' %(bfile,bfile,bfile,)]

        ##
        ## 1) sample QC
        ##

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#missing
        l_plink_cmds += ['--missing']

        ## http://pngu.mgh.harvard.edu/~purcell/plink/ibdibs.shtml#inbreeding
##        l_plink_cmds += ['--het --keep %s.imiss' %(bfile)] ## keep just to make sure it doesn't run until --missing finishes
        l_plink_cmds += ['--recodeA --keep %s.imiss' %(bfile)]
        ## xxx rename het_before/after to recodeA_before/after and include het calculation in recode_after

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#sexcheck
        l_plink_cmds += ['--check-sex --extract %s.X.SNPs' %(bfile)]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed
        s = '--make-bed \\\n'
        s += '--remove %s.sampleQC.samples \\\n' %(bfile)
        s += '--out %s.sampleQC \\\n' %(bfile)
        l_plink_cmds += [s]

        ##
        ## 2) SNP QC, pre Hardy, SNP exclusion (missingness)
        ##

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#missing
        s = '--missing \\\n'
        s +='--out %s.SNPQC \\\n' %(bfile,)
        s += '--remove %s.sampleQC.samples \\\n' %(bfile,)
        l_plink_cmds += [s]

        ##
        ## 3) SNP QC, pre Hardy, sample removal (relatedness)
        ##

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#freq
        s = '--freq \\\n'
        s += '--remove %s.sampleQC.samples \\\n' %(bfile)
        s += '--out %s.sampleQC \\\n' %(bfile)
        s += '--exclude %s.lmiss.SNPs \\\n' %(bfile)
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune
        s = s_indep_pairwise
        s += '--read-freq %s.sampleQC.frq \\\n' %(bfile)
        s += '--out %s.prehardy \\\n' %(bfile,)
        s += '--remove %s.sampleQC.samples \\\n' %(bfile)
        s += '--exclude %s.lmiss.SNPs \\\n' %(bfile)
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/ibdibs.shtml#genome
        s = self.write_genome_cmd(bfile, 'prehardy', 'sampleQC', bfile,)
        l_plink_cmds += [s]
##
##        ## http://pngu.mgh.harvard.edu/~purcell/plink/strat.shtml#cluster
##        s = s_cluster+'--read-genome %s.prehardy.genome \\\n' %(bfile,)
##        s += '--bfile %s \\\n' %(bfile,)
##        s += '--remove %s.sampleQC.samples \\\n' %(bfile,)
##        s += '--extract %s.prehardy.prune.in \\\n' %(bfile,)
##        s += '--out %s.prehardy \\\n' %(bfile,)
##        l_plink_cmds += [s]

        ##
        ## 4) SNP QC, pre Hardy, sample removal (EIGENSOFT)
        ## 5) SNP QC, pre Hardy, SNP exclusion (HWE)
        ##

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#hardy
        s = '--hardy \\\n'
        s += '--exclude %s.lmiss.SNPs \\\n' %(bfile,)
        s += '--remove %s.genome.prehardy.EIGENSOFT.samples \\\n' %(bfile,)
        l_plink_cmds += [s]

        ##
        ## 5) SNP QC, post Hardy, MDS plot without 1000G
        ##

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#freq
        s = '--freq \\\n'
        s += '--out %s.SNPQC \\\n' %(bfile,)
        s += '--remove %s.SNPQC.samples \\\n' %(bfile,)
        s += '--exclude %s.SNPQC.SNPs \\\n' %(bfile,)
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune
        s = s_indep_pairwise
        s += '--read-freq %s.SNPQC.frq \\\n' %(bfile,)
        s += '--out %s.posthardy \\\n' %(bfile,)
        s += '--remove %s.SNPQC.samples \\\n' %(bfile,)
        s += '--exclude %s.SNPQC.SNPs \\\n' %(bfile,)
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/ibdibs.shtml#genome
        s = self.write_genome_cmd(bfile, 'posthardy', 'SNPQC', bfile,)
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/strat.shtml#cluster
        s = s_cluster+'--read-genome %s.posthardy.genome \\\n' %(bfile,)
        s += '--bfile %s \\\n' %(bfile,)
        s += '--remove %s.SNPQC.samples \\\n' %(bfile,)
        s += '--extract %s.posthardy.prune.in \\\n' %(bfile,)
        s += '--out %s.posthardy \\\n' %(bfile,)
        l_plink_cmds += [s]

        ##
        ## 7) post SNP QC, make-bed
        ##

        ## http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed
        s = '--make-bed \\\n'
        s += '--remove %s.SNPQC.samples \\\n' %(bfile)
        s += '--out %s.SNPQC \\\n' %(bfile)
        s += '--exclude %s.SNPQC.SNPs \\\n' %(bfile,)
        l_plink_cmds += [s]

        ##
        ## SNP QC, chromosome X
        ##

        for sex in ['males','females']:

            ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#freq
##            s = '--freq \\\n'
##            s += '--out %s.X.%s \\\n' %(bfile,sex,)
##            s += '--remove %s.SNPQC.samples \\\n' %(bfile,)
##            s += '--extract %s.X.SNPs \\\n' %(bfile)
##            s += '--filter-%s \\\n' %(sex)
##            l_plink_cmds += [s]

            ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#missing
            s = '--missing \\\n'
            s += '--out %s.X.%s \\\n' %(bfile,sex,)
            s += '--remove %s.SNPQC.samples \\\n' %(bfile,)
            s += '--extract %s.X.SNPs \\\n' %(bfile)
            s += '--filter-%s \\\n' %(sex)
            l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#hardy
        s = '--hardy \\\n'
        s += '--out %s.X.females \\\n' %(bfile,)
        s += '--remove %s.genome.prehardy.EIGENSOFT.samples \\\n' %(bfile,)
        s += '--extract %s.X.SNPs \\\n' %(bfile)
        s += '--exclude %s.X.females.lmiss.SNPs \\\n' %(bfile,)
        s += '--filter-females \\\n'
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed
        s = '--make-bed \\\n'
        s += '--remove %s.SNPQC.samples \\\n' %(bfile)
        s += '--extract %s.X.SNPs \\\n' %(bfile)
        s += '--out %s.X \\\n' %(bfile,)
        s += '--exclude %s.X.lmiss.hwe.SNPs \\\n' %(bfile,)
        l_plink_cmds += [s]

        ##
        ## MDS with 1000G
        ##

        suffix = self.fn1000g

        ## http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#merge
        s = '--bmerge \\\n'
        s += '%s.bed \\\n' %(self.fp1000g,)
        s += '%s.bim \\\n' %(self.fp1000g,)
        s += '%s.fam \\\n' %(self.fp1000g,)
        s += '--make-bed \\\n'
        s += '--out %s.%s \\\n' %(bfile,suffix,)
        s += '--bfile %s.SNPQC \\\n' %(bfile,)
        s += '--remove %s.SNPQC.samples \\\n' %(bfile)
        s += '--extract %s.%s.comm.SNPs \\\n' %(bfile,self.fn1000g,)
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#freq
        s = '--freq \\\n'
        s += '--bfile %s.%s \\\n' %(bfile,suffix,)
        s += '--out %s.%s \\\n' %(bfile,suffix,)
        s += '--remove %s.SNPQC.samples \\\n' %(bfile,)
        s += '--extract %s.%s.comm.SNPs \\\n' %(bfile,suffix,)
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune
        s = s_indep_pairwise
        s += '--bfile %s.%s \\\n' %(bfile,suffix,)
        s += '--read-freq %s.%s.frq \\\n' %(bfile,suffix,)
        s += '--out %s.%s \\\n' %(bfile,suffix,)
        s += '--remove %s.SNPQC.samples \\\n' %(bfile,)
        s += '--extract %s.%s.comm.SNPs \\\n' %(bfile,suffix,)
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/ibdibs.shtml#genome
        s = self.write_genome_cmd(
            bfile, suffix, 'SNPQC', '%s.%s' %(bfile,suffix,),
            )
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/strat.shtml#cluster
        s = s_cluster+'--read-genome %s.%s.genome \\\n' %(bfile,suffix,)
        s += '--bfile %s.%s \\\n' %(bfile,suffix,)
        s += '--remove %s.SNPQC.samples \\\n' %(bfile,)
        s += '--extract %s.%s.prune.in \\\n' %(bfile,suffix,)
        s += '--out %s.%s \\\n' %(bfile,suffix,)
        l_plink_cmds += [s]

        return l_plink_cmds


    def write_genome_cmd(
        self, bfile, suffix, suffix_remove, bfile_in,
        ):

        out = 'genome/%s.%s.$i.$j' %(bfile,suffix,),
        remove = '%s.%s.samples' %(bfile,suffix_remove,),
        extract = '%s.%s.prune.in' %(bfile,suffix,),
        genome_lists = 'fam/%s.%s.fam.$i fam/%s.%s.fam.$j' %(
            bfile,suffix,bfile,suffix,
            )

        cmd_genome = '--genome \\\n'
        cmd_genome += '--genome-lists %s \\\n' %(genome_lists)
        cmd_genome += '--bfile %s \\\n' %(bfile_in)
        cmd_genome += '--out %s \\\n' %(out)
        cmd_genome += '--remove %s \\\n' %(remove)
        cmd_genome += '--extract %s' %(extract)

        return cmd_genome


    def check_file_in(self,fp_in,):

        minutes = 5
        if '.prune.in' in fp_in:
            minutes = 10

        bool_check = True
        if not os.path.isfile(fp_in):
            bool_check = False
##        elif not os.path.getsize(fp_in) > 0:
##            bool_check = False
        elif not time.time()-os.path.getmtime(fp_in) > minutes*60:
            bool_check = False

        return bool_check


    def extra_commands_before(
        self,bfile,plink_cmd,plink_cmd_full,in_prefix,out_prefix,
        ):

        l_cmds = []

        cmd = ''
        for out_suffix in self.d_out_suffix[plink_cmd]:
            fn_out = '%s.%s' %(out_prefix,out_suffix,)
            cmd += 'if [ -f %s ]; then\nexit\nfi\n' %(fn_out)
        l_cmds += [cmd]

        if plink_cmd == None:
            pass

        elif plink_cmd == 'genome':
            l_cmds = self.genome_before(bfile,in_prefix,out_prefix,)

##        elif plink_cmd == 'hardy':
##            l_cmds = self.hardy_before(bfile,out_prefix,)

        cmds = '\n\n'.join(l_cmds)
        
        return cmds


    def extra_commands_after(self,bfile,plink_cmd,in_prefix,out_prefix,):

        l_cmds = []

        if plink_cmd == 'missing':
            l_cmds = self.missing_after(bfile)

        elif plink_cmd == 'recodeA':
            l_cmds = self.recodeA_after(bfile)
            cmd = self.concatenate_sampleQC_remove_lists(bfile,)
            l_cmds += [cmd]

        elif plink_cmd == 'check-sex':
            cmd = '''awk '{if($5=="PROBLEM") print $1,$2}' '''
            cmd += ' %s.sexcheck > %s.sexcheck.samples' %(bfile, bfile,)
            l_cmds += [cmd]
            cmd = self.concatenate_sampleQC_remove_lists(bfile,)
            l_cmds += [cmd]

        elif plink_cmd == 'genome':
            l_cmds = self.genome_after(bfile,in_prefix,out_prefix,)

        elif plink_cmd == 'hardy':
            l_cmds = self.hardy_after(bfile,out_prefix,)

##        elif plink_cmd == 'cluster':
##            l_cmds += ['fi\n']

        elif plink_cmd == 'bmerge':
            l_cmds = self.bmerge_after(bfile,)

        cmds = '\n\n'.join(l_cmds)

        return cmds


    def hardy_after(self,bfile,out_prefix,):

        l_cmds = []

        ##
        ## autosome
        ##
        if out_prefix == bfile:

            ##
            ## if output exists
            ##
            cmd = 'if [ -s %s.hwe ]; then\n' %(bfile,)

            cmd += "awk '{if ($9 < %.1e) print $2}' %s.hwe > %s.hwe.SNPs\n" %(
                self.hwe_min, bfile,bfile,)
            cmd += 'cat %s.lmiss.SNPs %s.hwe.SNPs > %s.SNPQC.SNPs\n' %(
                bfile,bfile,bfile,
                )

            ##
            ## only merge SNPs found in 1000g and current dataset
            ##
            for fn_in,fn_out in [
                [self.fp1000g,self.fn1000g,],
                [bfile,bfile,],
                ]:
                ## 1) sort
                cmd += "\ncat %s.bim | awk '{print $2}' | sort > %s.SNPs.sorted" %(
                    fn_in,fn_out,)

            cmd += '\ncomm -12 %s.SNPs.sorted %s.SNPs.sorted > %s.%s.comm.SNPs' %(
                bfile,self.fn1000g,bfile,self.fn1000g,)

            cmd += '\nrm %s.SNPs.sorted %s.SNPs.sorted' %(bfile,self.fn1000g,)

########            cmd += '\nsort %s.posthardy.prune.in' %(bfile)
########            cmd += ' > %s.posthardy.prune.in.sorted' %(bfile)
########
########            cmd += '\ncomm -12 %s.posthardy.prune.in.sorted %s.%s.comm.SNPs > %s.%s.prune.in' %(
########                bfile,bfile,self.fn1000g,bfile,self.fn1000g,)
########
########            cmd += '\nrm %s.%s.comm.SNPs %s.posthardy.prune.in.sorted' %(
########                bfile,self.fn1000g,bfile,)

            ##
            ## fi
            ##
            cmd += '\nfi\n'
            l_cmds += [cmd]

        ##
        ## X
        ##
        if out_prefix == '%s.X.females' %(bfile):

            cmd = 'if [ -s %s.X.females.hwe -a -s %s.X.males.lmiss.SNPs ]; then\n' %(bfile,bfile,)

            cmd += 'cat %s.X.females.hwe | ' %(bfile)
            cmd += "awk '{if ($9 < %.1e) print $2}' > %s.X.females.hwe.SNPs\n" %(
                self.hwe_min, bfile,)
            cmd += 'cat '
            cmd += '%s.X.males.lmiss.SNPs ' %(bfile,)
            cmd += '%s.X.females.lmiss.SNPs ' %(bfile,)
            cmd += '%s.X.females.hwe.SNPs ' %(bfile,)
            cmd += '| sort -u '
            cmd += '> %s.X.lmiss.hwe.SNPs\n' %(bfile,)

            cmd += 'fi\n'
            l_cmds += [cmd]

        return l_cmds


    def EIGENSOFT(self,bfile,):

        ## http://computing.bio.cam.ac.uk/local/doc/eigenstrat.txt
        ## http://helix.nih.gov/Applications/README.eigenstrat

        ##Estimated running time of the smartpca program is 
        ##  2.5e-12 * nSNP * NSAMPLES^2 hours            if not removing outliers.
        ##  2.5e-12 * nSNP * NSAMPLES^2 hours * (1+m)    if m outlier removal iterations.
        ##Thus, under the default of up to 5 outlier removal iterations, running time is 
        ##  up to 1.5e-11 * nSNP * NSAMPLES^2 hours.

        cmd = ''

        if self.project == 'uganda_gwas':
            bool_eigensoft = True
        elif self.project == 'agv':
            bool_eigensoft = False ## cp8 17jan2013
            
        if bool_eigensoft == True:

            ##
            ## commands executed now
            ##
            
            cmd = ''
            cmd += 'cat %s.fam' %(bfile)
    ##        cmd += "awk '{print $1,$2,substr($1,12,10),substr($2,12,10)}'"
            ## init awk
            cmd += " | awk '{"
            cmd += 'sub(/-RECLUSTER/,"",$1);sub(/-RECLUSTER/,"",$2)'
            cmd += ';print $1,$2,'
            cmd += 'substr($1,length($1)-9,10),substr($2,length($2)-9,10)'
            ## term awk
            cmd += "}'"
            cmd += ' > %s.recoded.txt' %(bfile)
            cmd += '\n\n'
            self.execmd(cmd)

            cmd = ''
            cmd += 'cat %s.recoded.txt' %(bfile)
            cmd += " | awk '{print $4}'"
            cmd += ' | sort'
            cmd += ' | uniq -d'
            if os.popen(cmd).read().strip() != '' and os.popen(cmd).read().strip() != 'APP5212239':
                print cmd
                stop

            self.write_EIGENSOFT_parameter_file(bfile)

            ##
            ## command to be appended to shell script
            ##

            cmd = ''
            cmd += 'echo EIGENSOFT\n\n'
            
            cmd += 'cat %s.lmiss.SNPs %s.ldregions.SNPs > %s.EIGENSOFT.exclude.SNPs\n\n' %(
                bfile,bfile,bfile,)

            ## format IIDs of sample removal list
            cmd += 'cat %s.genome.samples' %(bfile)
            cmd += " | awk '{"
            cmd += 'sub(/-RECLUSTER/,"",$2);sub(/-RECLUSTER/,"",$2); print '
            cmd += 'substr($1,length($1)-9,10),substr($2,length($2)-9,10)'
            cmd += "}'"
            cmd += ' > %s.genome.EIGENSOFT.samples\n\n' %(bfile)

            cmd += 'if [ ! -f %s.EIGENSOFT.bed ]; then\n\n' %(bfile)
            cmd += 'plink \\\n'
            cmd += '--bfile %s \\\n' %(bfile,)
            cmd += '--exclude %s.EIGENSOFT.exclude.SNPs \\\n' %(bfile,)
            cmd += '--extract %s.prehardy.prune.in \\\n' %(bfile,)
            cmd += '--update-ids %s.recoded.txt \\\n' %(bfile,)
            cmd += '--remove %s.genome.EIGENSOFT.samples \\\n' %(bfile,)
            cmd += '--make-bed \\\n'
            cmd += '--out %s.EIGENSOFT \\\n' %(bfile)
            cmd += '\nfi'
            cmd += '\n\n'

            ##
            cmd += 'cat %s.genome.prehardy.samples' %(bfile)
            cmd += " | awk '{"
            cmd += 'sub(/-RECLUSTER/,"",$2);sub(/-RECLUSTER/,"",$2); print '
            cmd += 'substr($1,length($1)-9,10),substr($2,length($2)-9,10)'
            cmd += "}'"
            cmd += ' | sort -k2,2'
            cmd += ' > %s.genome.prehardy.samples.sorted\n\n' %(bfile)
            ##
            cmd += 'cat %s.EIGENSOFT.fam' %(bfile)
            cmd += ' | sort -k2,2'
            cmd += ' > %s.EIGENSOFT.fam.sorted\n\n' %(bfile)
            ##
            for s1,s2,s3 in [
                [' -v2','fou','>',],
                ['','nof','>>',],
                ]:
                cmd += 'join -1 2 -2 2 %s -o 2.1,2.2,2.3,2.4,2.5,2.6' %(s1)
                cmd += ' %s.genome.prehardy.samples.sorted %s.EIGENSOFT.fam.sorted' %(bfile,bfile,)
                cmd += " | awk '{"
                cmd += 'sub(/-RECLUSTER/,"",$2);sub(/-RECLUSTER/,"",$2); print '
                cmd += 'substr($1,length($2)-9,10),substr($2,length($2)-9,10),$3,$4,$5,"%s"' %(s2)
                cmd += "}'"
                cmd += ' %s %s.EIGENSOFT.fam\n\n' %(s3,bfile)
  
            cmd += 'echo "fou" > %s.EIGENSOFT.poplist\n' %(bfile)

            ##
            cmd += 'rm %s.EIGENSOFT.fam.sorted %s.genome.prehardy.samples.sorted\n\n' %(bfile,bfile,)

            ## run with parameter file options
            cmd += '%s \\\n' %(self.eigensoft)
            cmd += '-p %s.par\n\n' %(bfile)

##            ## run with command line options
##            cmd += '/software/varinf/bin/eigensoft/bin/smartpca.perl \\\n'
##            cmd += ' -i %s.EIGENSOFT.bed' %(bfile)
##            cmd += ' -a %s.EIGENSOFT.bim' %(bfile)
##            cmd += ' -b %s.EIGENSOFT.fam' %(bfile)
##            cmd += ' -o %s.EIGENSOFT.evec' %(bfile)
##            cmd += ' -p %s.EIGENSOFT.plot' %(bfile)
##            cmd += ' -e %s.EIGENSOFT.eval' %(bfile)
##            cmd += ' -l %s.EIGENSOFT.log' %(bfile)
##            cmd += ' -w %s.EIGENSOFT.poplist' %(bfile)
##            cmd += '\n\n'

##            cat omni2.5-8_20120809_gwa_uganda_gtu_flipped.outlier | awk '{print $3}' | awk -F ":" '{print $1}' > tmp ; fgrep -f tmp omni2.5-8_20120809_gwa_uganda_gtu_flipped.prehardy.genome | awk '{if($10>0.05) print}' | head

        else:

            fd = open('%s.EIGENSOFT.outlier' %(bfile),'w')
            fd.close()

        ##
        ## parse outliers and append
        ##
        ## print EIGENSOFT sample names to file
        cmd += 'cat %s.EIGENSOFT.outlier' %(bfile)
        cmd += " | awk '{print $3}'"
        cmd += ''' | awk 'BEGIN{FS=":"}{print $2}' '''
        cmd += ' > %s.outlier.EIGENSOFTsamples' %(bfile)
        ## convert EIGENSOFT sample names to original sample names
        cmd += ' ; fgrep -w -f %s.outlier.EIGENSOFTsamples %s.recoded.txt' %(bfile,bfile)
        cmd += " | awk '{print $1,$2}'"
        cmd += ' > %s.EIGENSOFT.samples' %(bfile)
        ## before hardy
        cmd += ' ; cat %s.genome.prehardy.samples %s.EIGENSOFT.samples' %(bfile,bfile)
        cmd += ' > %s.genome.prehardy.EIGENSOFT.samples' %(bfile)
        ## after hardy
        cmd += ' ; cat %s.genome.samples %s.EIGENSOFT.samples' %(bfile,bfile)
        cmd += ' > %s.SNPQC.samples\n\n' %(bfile)
        ## clean-up
        cmd += 'rm %s.outlier.EIGENSOFTsamples\n\n' %(bfile)

        return cmd


##    def cluster_after(self,bfile,in_prefix,out_prefix,):
##
##        l_cmds = []
##        cmd = '\n;\n'
##
##        cmd += '"'
##        l_cmds += [cmd]
##
##        return l_cmds


    def write_EIGENSOFT_parameter_file(self,bfile):

        ## http://helix.nih.gov/Applications/README.popgen
        ## http://computing.bio.cam.ac.uk/local/doc/popgen.txt
        ##
        ## input
        par = 'genotypename: %s.EIGENSOFT.bed\n' %(bfile)
        par += 'snpname: %s.EIGENSOFT.bim\n' %(bfile)
        par += 'indivname: %s.EIGENSOFT.fam\n' %(bfile)
        ##
        ## output
        par += 'evecoutname: %s.EIGENSOFT.evec\n' %(bfile)
        par += 'evaloutname: %s.EIGENSOFT.eval\n' %(bfile)
        ##
        ## optional parameters
        ##
        ## numoutlieriter: maximum number of outlier removal iterations.
        ## Default is 5.  To turn off outlier removal, set this parameter to 0.
##        par += 'numoutlieriter: 5\n'
        ##
        ## outliersigmathresh: number of standard deviations which an individual must 
        ## exceed, along one of the top (numoutlierevec) principal components, in
        ## order for that individual to be removed as an outlier.  Default is 6.0.
##        par += 'outliersigmathresh: 4\n'
        ##
        ## outlieroutname: output logfile of outlier individuals removed. If not specified,
        ## smartpca will print this information to stdout, which is the default.
        par += 'outlieroutname: %s.EIGENSOFT.outlier\n' %(bfile)
        ##
        ## numoutlierevec: number of principal components along which to
        ## remove outliers during each outlier removal iteration.  Default is 10.
##        par += 'numoutlierevec: 2\n'
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

        par += 'poplistname: %s.EIGENSOFT.poplist\n' %(bfile)

        fd = open('%s.par' %(bfile),'w')
        fd.write(par)
        fd.close()

        return


    def genome_before(self,bfile,in_prefix,out_prefix,):

        l_cmds = []

        if out_prefix == '%s.prehardy' %(bfile):
            fn_samples_remove = '%s.sampleQC.samples' %(bfile)
        elif out_prefix == '%s.posthardy' %(bfile):
            fn_samples_remove = '%s.SNPQC.samples' %(bfile)
        elif out_prefix == '%s.%s' %(bfile,self.fn1000g,):
            fn_samples_remove = '%s.SNPQC.samples' %(bfile)
        else:
            print in_prefix
            print out_prefix
            print self.fn1000g
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
        if out_prefix != '%s.%s' %(bfile,self.fn1000g,):
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


    def genome_after(self,bfile,in_prefix,out_prefix,):

        l_cmds = []

        ##
        ## continuation of nested loop and command from genome_before
        ##
        cmd = ''

        ## append to command
        if self.bool_run_all == True:
            cmd += ';'
            cmd += self.cmd_rerun(bfile,'genome',nl=';\n')
        ## terminate command
        cmd += "\n'\n\n"
        ## evaluate command
        cmd_LSF = self.append_LSF(bfile,'genome',JOBID='genome.${i}.${j}',)
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

        if out_prefix == '%s.prehardy' %(bfile):

            ##
            ## run IBD prune
            ##

            ## disallow duplicates (i.e. IBD>0.90)
            if self.project == 'uganda_gwas':
                cmd += 'python %s/QC_IBD_prune.py ' %(
                    os.path.dirname(sys.argv[0]),
                    )
                cmd += '--pi_hat_max %.2f --genome %s.prehardy --imiss %s --out %s\n\n' %(
                    .9,bfile,bfile,bfile,
                    )

            ## disallow related samples for HWE
            ## write minimal sample exclusion list
            cmd += 'python %s/QC_IBD_prune.py ' %(
                os.path.dirname(sys.argv[0]),
                )
            cmd += '--pi_hat_max %.2f --genome %s.prehardy --imiss %s --out %s\n\n' %(
                self.pi_hat_max,bfile,bfile,bfile,
                )

            ##
            ## concatenate sample removal lists
            ##

            ## after HWE
            if self.project == 'uganda_gwas':
                ## concatenate sample exclusion lists (i.e. IBD>0.90)
                cmd += 'cat %s.sampleQC.samples %s.genome.%.2f.samples' %(
                    bfile,bfile,.9,
                    )
                cmd += ' > %s.genome.samples\n\n' %(bfile,)
            elif self.project == 'agv':
                ## concatenate sample exclusion lists (i.e. IBD>low threshold)
                cmd += 'cat %s.sampleQC.samples %s.genome.%.2f.samples' %(
                    bfile,bfile,self.pi_hat_max,
                    )
                cmd += ' > %s.genome.samples\n\n' %(bfile,)
            else:
                stop_unknown_project

            ## before HWE
            cmd += 'cat %s.sampleQC.samples %s.genome.%.2f.samples' %(
                bfile,bfile,self.pi_hat_max,
                )
            cmd += ' > %s.genome.prehardy.samples\n' %(bfile,)

            ##
            ## EIGENSOFT before HWE
            ##
            cmd_EIGENSOFT = self.EIGENSOFT(bfile,)
            cmd += '\n##\n## EIGENSOFT\n##\n'+cmd_EIGENSOFT+'\n\n'

        if self.bool_run_all == True:
            cmd += self.cmd_rerun(bfile,'genome',)

##        cmd += self.mail(bfile,'genome',out_prefix,'finished',)

        l_cmds += [cmd]

        return l_cmds


    def recodeA_after(self,bfile,):

        l_cmds = []

        cmd = 'cat %s.raw' %(bfile)
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
        cmd += ' > %s.het' %(bfile)
        l_cmds += [cmd]

        ## append header
        cmd = "sed -i '1i FID IID N(HET) N(NM) HET' %s.het" %(bfile)
##        cmd = 'echo "FID IID N(HET) N(NM) HET" > %s.het' %(bfile)
        l_cmds += [cmd]
        
        ## a) calculate heterozygosity mean and stddev
        l, mean,average, het_min, het_max = self.het2stddev(
                bfile,self.threshold_imiss,bool_execute=False,)
        l_cmds += l

        l_cmds += ['mean=$(echo $s | cut -d " " -f1)']
        l_cmds += ['stddev=$(echo $s | cut -d " " -f2)']

        l_cmds += ['echo mean $mean stddev $stddev']

        ## b) write list of samples outside mean+-3stddev
        cmd = 'awk -v mean=$mean -v stddev=$stddev '
        cmd += " 'NR>1 {if( "
        cmd += '$5<mean-%i*stddev || $5>mean+%i*stddev ' %(
            self.threshold_het_stddev, self.threshold_het_stddev,
            )
        cmd += ") print $1,$2}' "
        cmd += ' %s.het > %s.het.samples' %(bfile, bfile,)
        l_cmds += [cmd]

        return l_cmds


    def missing_after(self,bfile,):

        l_cmds = []

        ##
        ## imiss.samples
        ##
        cmd = 'if [ -s %s.imiss -a ! -f %s.imiss.samples ]\n' %(bfile,bfile,)
        cmd += 'then\n'
        ## exclusion list
        cmd += "awk 'NR>1 {if($6>%f) print $1,$2}' %s.imiss " %(
            1-self.threshold_imiss,
            bfile,
            )
        cmd += ' > %s.imiss.samples\n' %(bfile)
##        cmd += 'rm %s.lmiss\n' %(bfile)
        cmd += self.concatenate_sampleQC_remove_lists(bfile,)
        cmd += 'fi\n'
        l_cmds += [cmd]

        ##
        ## lmiss.SNPs
        ##
        cmd = 'if [ -s %s.SNPQC.lmiss -a ! -f %s.lmiss.SNPs ]\n' %(bfile,bfile,)
        cmd += 'then\n'
        cmd += 'cat %s.SNPQC.lmiss' %(bfile)
        cmd += ' | '
        cmd += "awk 'NR>1 {if ((1-$5)<%.2f) print $2}' " %(
            self.threshold_lmiss,
            )
        cmd += ' > '
        cmd += '%s.lmiss.SNPs\n' %(bfile)
        cmd += 'rm %s.SNPQC.imiss\n' %(bfile)
        cmd += 'fi\n'
        l_cmds += [cmd]

        ##
        ## lmiss.X.SNPs
        ##
        for sex in ['males','females',]:
            cmd = 'if [ -s %s.X.%s.lmiss -a ! -f %s.X.%s.lmiss.SNPs ]\n' %(
                bfile,sex,bfile,sex,)
            cmd += 'then\n'
            cmd += 'cat %s.X.%s.lmiss' %(bfile,sex,)
            cmd += ' | '
            cmd += "awk 'NR>1 {if ((1-$5)<%.2f) print $2}' " %(
                self.threshold_lmiss,
                )
            cmd += ' > '
            cmd += '%s.X.%s.lmiss.SNPs\n' %(bfile,sex,)
            cmd += 'rm %s.X.%s.imiss\n' %(bfile,sex,)
            cmd += 'fi\n'
            l_cmds += [cmd]

        return l_cmds


    def bmerge_after(self,bfile,):

        l_cmds = []

        cmd = 'rm %s.%s.prune.in' %(bfile,self.fn1000g,)
        l_cmds += [cmd]

        return l_cmds


    def concatenate_sampleQC_remove_lists(self,bfile,):

        cmd = 'if [ '
        ## check file out
        cmd += '! -s %s.sampleQC.samples ' %(bfile)
        ## check file in
        cmd += '-a -s %s.imiss ' %(bfile)
        cmd += '-a -s %s.het ' %(bfile)
        cmd += '-a -s %s.sexcheck ' %(bfile)
        cmd += '-a -f %s.imiss.samples ' %(bfile)
        cmd += '-a -f %s.het.samples ' %(bfile)
        cmd += '-a -f %s.sexcheck.samples ' %(bfile)
        cmd += ']\n'
        cmd += 'then\n'
        ## concatenate sample files
        cmd += 'cat '
        cmd += '%s.imiss.samples ' %(bfile,)
        cmd += '%s.het.samples ' %(bfile,)
        cmd += '%s.sexcheck.samples ' %(bfile,)
        cmd += '| sort -u > %s.sampleQC.samples\n' %(bfile,)
        cmd += 'fi\n'

        return cmd
    

    def histogram_het(self,bfile,bool_with_stddev=True,):

        if os.path.isfile('%s.het.png' %(bfile)):
            return

        print 'histogram het', bfile

        n_samples = int(os.popen('cat %s.het | wc -l' %(bfile)).read())-1
        cmd = 'cat %s.imiss.samples | wc -l' %(bfile)
        n_samples -= int(os.popen(cmd).read())
        cmd = 'cat %s.autosomes.SNPs | wc -l' %(bfile)
        n_SNPs = int(os.popen(cmd).read())

        x_step=0.0005

        if bool_with_stddev == True:

            fn_plot = '%s.het.joined' %(bfile)

            l_cmds,average,stddev,het_min,het_max = self.het2stddev(
                bfile,self.threshold_imiss,bool_execute=True,bool_remove=False,)

            l_arrows = [
                'set arrow from %f,0 to %f,graph(1) nohead lc 0\n' %(
                    average-3*stddev,average-3*stddev,
                    ),
                'set arrow from %f,0 to %f,graph(1) nohead lc 0\n' %(
                    average+3*stddev,average+3*stddev,
                    ),
                ]

            s_plot = 'plot [%f:%f]"%s" ' %(
                min(het_min,average-3*stddev,)-2*x_step,
                max(het_max,average+3*stddev,)+2*x_step,
                fn_plot,
                )
            s_plot += 'u (hist($5,width)):(1.0) smooth freq w boxes '
            s_plot += 'lc rgb"red" notitle\n'

            s_average = '%.4f' %(average)
            s_stddev = '%.4f' %(stddev)

            pass

        else:

            fn_plot = '%s.het' %(bfile)

            s_plot = 'plot "%s" ' %(fn_plot)
            s_plot += 'u (hist($5,width)):(1.0) smooth freq w boxes '
            s_plot += 'lc rgb"red" notitle\n'
            l_arrows = None

            s_average = 'N/A'
            s_stddev = 'N/A'

        gnuplot.histogram2(
            '%s.het' %(bfile),
##            x_min=het_min,x_max=het_max,#tic_step=0.05,
            x_step = x_step,
            s_plot = s_plot,
            xlabel='heterozygosity',
            title='%s (n_{samples}=%i, n_{SNPs}=%i, mean=%s, SD=%s)' %(
                bfile.replace('_','\\\\_'),n_samples,n_SNPs,
                s_average, s_stddev,
                ),
            lines_extra = l_arrows,
            color = 'red',
            )

        return


    def histogram_genome(self,bfile,):

##awk 'NR>1{if($10!="PI_HAT") print $0;}' $chip.genome > $chip.genome2

        if os.path.isfile('%s.prehardy.genome.png' %(bfile)):
            return

        print 'histogram PI_HAT', bfile

        ## no data to plot
        if int(os.popen('head %s.prehardy.genome | wc -l' %(bfile,)).read()) == 1:
            return

##        cmd = 'cat %s.prehardy.genome | wc -l' %(bfile)
##        if self.verbose == True:
##            print cmd
##        n = int(os.popen(cmd).read())-1
##        ## http://simple.wikipedia.org/wiki/Quadratic_equation
##        n_samples1 = int(((--.5)+math.sqrt((-.5)**2-4*.5*-n))/(2*.5))

        n_samples2 = int(os.popen('cat %s.fam | wc -l' %(bfile)).read())
        n_samples2 -= int(os.popen('cat %s.sampleQC.samples | wc -l' %(bfile)).read())

##        if n_samples1 != n_samples2:
##            print n_samples1
##            print n_samples2
##            stop
##        else:
##            n_samples = n_samples1 = n_samples2

        n_samples = n_samples2

        n_SNPs = int(os.popen('cat %s.prehardy.prune.in | wc -l' %(bfile)).read())

        ##
        ## find highest pairwise IBD per sample
        ##
        d_IBD = {}
        ## populate dictionary with keys
        fd = open('%s.fam' %(bfile))
        lines = fd.readlines()
        fd.close()
        for line in lines:
            IID = line.split()[0]
            d_IBD[IID] = -1

        fd = open('%s.prehardy.genome' %(bfile),'r')
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

        l = ['%s\n' %(v) for v in d_IBD.values()]
        fd = open('%s.genome.max' %(bfile),'w')
        fd.writelines(l)
        fd.close()

        for prefix,column,xlabel,x_max,x_step in [
            [
                '%s.prehardy.genome' %(bfile),
                '$10','PI CIRCUMFLEX',self.pi_hat_max,0.002,
                ],
            ['%s.genome.max' %(bfile),'$1','PI CIRCUMFLEX MAX',1,0.01,],
            ]:

            gnuplot.histogram2(
                prefix,
                x_min=0,
                x_step=x_step,
                x_max=x_max,
                column=column,
                xlabel=xlabel,
                title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                    bfile.replace('_','\\\\_'),n_samples,n_SNPs,
                    ),
                color = 'green',
                )

        os.remove('%s.genome.max' %(bfile))

        return


    def histogram_frq(self,bfile,):

        if os.path.isfile('%s.sampleQC.frq.png' %(bfile)):
            return

        print 'histogram MAF', bfile

        cmd = 'cat %s.fam | wc -l' %(bfile)
        n_samples = int(os.popen(cmd).read())
        cmd = 'cat %s.sampleQC.samples | wc -l' %(bfile,)
        n_samples -= int(os.popen(cmd).read())
        cmd = 'cat %s.sampleQC.frq | wc -l' %(bfile)
        n_SNPs = int(os.popen(cmd).read())-1

        if n_samples < 200:
            x_step=.01
        else:
            x_step=.005
    
        gnuplot.histogram2(
            '%s.sampleQC.frq' %(bfile),
            x_step=x_step,
            x_max=0.5,
            column='$5',
##            s_plot = 'plot "%s.frq" u (hist($5,width)):(1.0) smooth freq w boxes lc rgb"blue" notitle\n' %(
##                bfile,
##                ),
            xlabel='MAF',
            title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                bfile.replace('_','\\\\_'),n_samples,n_SNPs,
                ),
            color = 'blue',
            )

        return


    def histogram_lmiss(self,bfile,):

        if os.path.isfile('%s.lmiss.png' %(bfile)):
            return

        print 'histogram call lmiss', bfile

        n_samples = int(os.popen('cat %s.fam | wc -l' %(bfile)).read())
        n_samples -= int(os.popen('cat %s.sampleQC.samples | wc -l' %(bfile,)).read())
        n_SNPs = int(os.popen('cat %s.SNPQC.lmiss | wc -l' %(bfile)).read())
        n_SNPs -= 1

        if n_samples > 200:
            x_step = .0002
        elif n_samples < 20:
            x_step = .1
        else:
            x_step = .01
    
        gnuplot.histogram2(
            '%s.SNPQC.lmiss' %(bfile),
            prefix_out='%s.lmiss' %(bfile),
            column='(1-$5)',
            x_step=x_step,
##            x_min=self.threshold_lmiss,
            x_min=0.95,
            x_max=1,
            xlabel='call rate, SNPs',
            title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                bfile.replace('_','\\\\_'),n_samples,n_SNPs,
                ),
            color = 'purple',
            )

        return


    def histogram_imiss(self,bfile,):

        if os.path.isfile('%s.imiss.png' %(bfile)):
            return

        print 'histogram call imiss', bfile

        n_samples = int(os.popen('cat %s.imiss | wc -l' %(bfile)).read())-1
        if n_samples < 200:
            x_step=.005
        else:
            x_step=.002

        n_SNPs = int(os.popen('cat %s.autosomes.SNPs | wc -l' %(bfile)).read())

        cmd = "awk 'NR>1{print 1-$6}' %s.imiss" %(bfile)
        x_min = min(
            0.95,
            self.threshold_imiss,
            min([float(s) for s in os.popen(cmd).readlines()]),
            )
   
        gnuplot.histogram2(
            '%s.imiss' %(bfile),
            column='(1-$6)',
##            x_step=.0001,
            x_step=x_step,
            x_max=1,
            x_min=x_min,
##            x_min=self.threshold_imiss,
            xlabel='call rate, samples',
            title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                bfile.replace('_','\\\\_'),n_samples,n_SNPs,
                ),
            color = 'orange',
            )

        return


    def init(self,bfile,):

        ##
        ## create dirs
        ##
        for dn in ['stdout','stderr','fam','genome',]:
            if not os.path.isdir(dn):
                os.mkdir(dn)
                pass
            continue

        ##
        ##
        ##
        if not os.path.isfile('%s.bed' %(self.bfile)):
            print 'bed not found:', self.bfile
            sys.exit(0)

        ##
        ## write options in use to file
        ##
        d_options = vars(self.opts)
        l = []
        for k,v in d_options.items():
            l += [[k,v,]]
            continue
        l.sort()
        s = ''
        for k,v in l:
            s += '%-20s\t%s\n' %(k,v,)
            continue
        fd = open('%s.options' %(self.bfile),'w')
        fd.write(s)
        fd.close()

        ##
        ## write lists of autosomal SNPs and X chromosome SNPs
        ## the latter is used for --check-sex
        ## subtract concatenated SNP exclusion list
        ##
        self.write_SNPs_by_chromosome(bfile,)

        ##
        ## convert range exclusion list to SNP exclusion list
        ## to allow --exclude and --extract to be combined
        ## which they cannot be, when --range is submitted for one of them
        ##
        self.convert_range_exclusion_to_SNP_exclusion(bfile,)

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


    def convert_range_exclusion_to_SNP_exclusion(self,bfile,):

        fn_out = '%s.ldregions.SNPs' %(bfile)

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
            if not chrom in d_ranges.keys():
                d_ranges[chrom] = []
            d_ranges[chrom] += [[i1,i2,]]

        fd = open('%s.bim' %(bfile),'r')
        lines = []
        for line in fd:
            ## split line
            l = line.split()
            ## parse chromosome
            chrom = l[0]
            if not chrom in d_ranges.keys():
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

        if '%s' in cmd:
            print cmd
            stop
        if self.bool_verbose == True:
            print
            print inspect.stack()[1][3]
            print cmd
        os.system(cmd)

        return


    def write_SNPs_by_chromosome(self,bfile,):

        '''somehow I can never remember what this function does.
maybe I should rename it.'''

        for out_prefix, condition in [
            ['X','$1==23',],
            ['autosomes','$1>=1 && $1<=22',],
            ]:

            fn_out = '%s.%s.SNPs' %(bfile,out_prefix,)

            ## continue if file already exists
            if os.path.isfile(fn_out):
                continue

            ## generate inclusion list
            cmd = "cat %s.bim | awk '{if(%s) print $2}'  > %s" %(
                bfile,condition,fn_out,
                )
            self.execmd(cmd)

        return


    def parse_options(self,):

        parser = optparse.OptionParser()

        parser.add_option(
            '-p', '--project', dest='project',
            help='The functionality of the script depends on the proejct you choose.',
            metavar='STR',
            )

        parser.add_option(
            '--pi_hat_max', '--pi_circumflex_max', '--threshold_ibd', '--threshold_pi_hat',
            dest='threshold_pi_hat_max',
            help='The HWE step is dependent on the chose PI CIRCUMFLEX MAX',
            metavar='FLOAT',default=0.05,
            )

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune
        parser.add_option(
            '--indepWindow', dest='indepWindow',
            help='The SNP pruning step (indep-pairwise) is dependent on the chosen value.',
            metavar='INT',default=50,
            )

        parser.add_option(
            '--indepShift', dest='indepShift',
            help='The SNP pruning step (indep-pairwise) is dependent on the chosen value.',
            metavar='INT',default=5,
            )

        parser.add_option(
            '--Rsquared', dest='indepRsquared',
            help='The SNP pruning step (indep-pairwise) is dependent on the chosen value.',
            metavar='FLOAT',default=.2,
            )

        parser.add_option(
            '--threshold_het_stddev',
            dest='threshold_het_stddev',
            help='Used to determine, which samples fail the heterozygosity check.',
            metavar='INT',default=3,
            )

        ## http://pngu.mgh.harvard.edu/~purcell/plink/thresh.shtml
        parser.add_option(
            '-i', '--mind', '--imiss', '--threshold_imiss', '--threshold_callrate_samples', '--threshold_callrate_sample', '--threshold_missing_sample',
            dest='threshold_imiss_min',
            help='The sample call rate calculation (--missing) is dependent on this value.',
            metavar='FLOAT',default=.97,
            )

        ## http://pngu.mgh.harvard.edu/~purcell/plink/thresh.shtml
        parser.add_option(
            '-l', '--geno', '--lmiss', '--threshold_lmiss', '--threshold_callrate_SNPs', '--threshold_callrate_SNP', '--threshold_missing_SNP',
            dest='threshold_lmiss_min',
            help='The SNP call rate calculation (--missing) is dependent on this value.',
            metavar='FLOAT',default=.97,
            )

        ## http://pngu.mgh.harvard.edu/~purcell/plink/thresh.shtml
        parser.add_option(
            '--maf','--maf_min','--threshold_maf',
            dest='threshold_maf_min',
            help='Allele frequency threshold',
            metavar='FLOAT',default=.05,
            )

        ## http://pngu.mgh.harvard.edu/~purcell/plink/thresh.shtml
        parser.add_option(
            '--hwe','--hwe_min', '--threshold_hwe_min', '--threshold_hwe',
            dest='threshold_hwe_min',
            help='HWe threshold used to determine failing SNPs after the HWE step.',
            metavar='FLOAT',default=10**-8,
            )

        parser.add_option(
            '-b', '--bfile',
            dest='bfile',
            help='Specify an input file.',
            metavar='FILE',
            )

        parser.add_option("-v", '--verbose', action="store_true", dest="verbose", default=True)
        parser.add_option("-q", '--silent', '--quiet', action="store_false", dest="verbose")

        ## parse option arguments
        (opts, args) = parser.parse_args()

        if opts.bfile is None:
            parser.error('--bfile not specified')
        if opts.project is None:
            parser.error('--project not specified')

        return opts


    def __init__(self,):

        ##
        ## options
        ##

        opts = self.parse_options()

        self.bfile = opts.bfile
        self.project = str(opts.project)
        self.bool_verbose = self.verbose = bool(opts.verbose)
        ## indep-pairwise
        self.indepWindow = int(opts.indepWindow)
        self.indepShift = int(opts.indepShift)
        self.Rsquared = float(opts.indepRsquared)
        ## het
        self.threshold_het_stddev = int(opts.threshold_het_stddev)
        ## freq
        self.maf_min = float(opts.threshold_maf_min)
        ## hwe
        self.hwe_min = float(opts.threshold_hwe_min)
        ## genome
        self.pi_hat_max = float(opts.threshold_pi_hat_max)
        ## missing
        self.threshold_imiss = float(opts.threshold_imiss_min)
        self.threshold_lmiss = float(opts.threshold_lmiss_min)

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

        ## perform various sanity checks of file lengths etc.
        self.bool_check_file_io = True

        ##
        ## paths
        ##
        
        self.dn1000g = '/lustre/scratch107/projects/uganda_gwas/users/tc9/QC'
        self.fn1000g = '1000G_True' ## quad, build37, excl "chip effect" SNPs
##        self.dn1000g = '/lustre/scratch107/projects/agv/users/tc9/QC/'
##        self.fn1000g = 'omni25_b37_1kg'
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

        self.opts = opts

        return

if __name__ == '__main__':
    instance_main = main()
    instance_main.main(instance_main.bfile)
