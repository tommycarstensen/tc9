#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

## built-ins
import os,sys,time,math,inspect,optparse
## add-ons
import numpy,pwd
## http://docs.scipy.org/doc/scipy/reference/stats.html
from scipy import stats
## other
sys.path.append('/nfs/users/nfs_t/tc9/github/tc9/misc')
import gnuplot

class main:

    ## Essential document describing the order of operations in PLINK:
    ## http://pngu.mgh.harvard.edu/~purcell/plink/flow.shtml

    def main(self,bfile,):

        self.check_logs() ## tmp!!!
        sys.exit(0) ## tmp!!!

        self.init(bfile,)

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
            if os.path.getsize('stderr/%s' %(fn)) > 0:
                print 'stderr/%s' %(fn)
                stop
            continue
        l = os.listdir('%s/stdout' %(os.getcwd()))
        for fn in l:
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
        cmd = 'cat %s.SNPQC.samples | wc -l' %(bfile)
        n_samples -= int(os.popen(cmd).read())

        n_SNPs = int(os.popen('cat %s.hwe.dat | wc -l' %(bfile)).read())

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
            [['flipped.bed'],self.table_flip,],
            [['sampleQC.frq'],self.table_frq,],
            [['SNPQC.lmiss'],self.table_lmiss,],
            [['prehardy.genome'],self.table_genome,],
            [['hwe'],self.table_hwe,],
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

        cmd = '''awk 'NR>1{print $5}' %s.sampleQC.frq''' %(bfile,)
        if self.verbose == True:
            print cmd
        l_frq = []
        for s in os.popen(cmd).readlines():
            if s.strip() == 'NA': ## what is this caused by???
                pass
            else:
                l_frq += [float(s)]

##        l_scores = [.0,.1,.2,.3,.4,.5,1.,2.,3.,4.,5.,10.,20.,30.,40.,50.,]
        l_scores = [.0,.1,.2,.5,1.,2.,5.,10.,20.,50.,]

        s = '%s' %(bfile)
        for score in l_scores:
            score /= 100
            n = int(len(l_frq)*stats.percentileofscore(l_frq, score, 'weak',))/100
            if self.verbose == True:
                print score, n
            s += '\t%i' %(n)
        s += '\n'

        if not os.path.isfile('frq.table'):
            l = [str(score) for score in l_scores]
            s_header = '#\t'+'\t'.join(l)+'\n'
            fd = open('frq.table','w')
            fd.write(s_header)
            fd.close()

        fd = open('frq.table','a')
        fd.write(s)
        fd.close()

        os.system('sort frq.table -k1,1 -u -o frq.table')

        return


    def table_sexcheck(self,bfile,):

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

        fn_out = 'sexcheck.table'

        if not os.path.isfile(fn_out):
            fd = open(fn_out,'w')
            fd.write('# IID PEDSEX SNPSEX F\n')
            fd.close()

        fd = open(fn_out,'a')
        fd.writelines(lines)
        fd.close()

        os.system('sort -u %s | sort -k1,1 -o %s' %(fn_out,fn_out,))


        return l_sexcheck


    def table_hwe(self,bfile,):

        ## do loop over thresholds even though it's not optimal...
        s = '%s' %(bfile)
        for i in xrange(-1,9+1):
            f = 10**-i
            cmd = 'grep ALL %s.hwe' %(bfile,)
            cmd += " | awk '{if ($9 < %.9f) print $2}' " %(f,)
            cmd += ' | wc -l'
            n = int(os.popen(cmd).read())
            print i, n
            s += '\t%i' %(n)
        s += '\n'

        if not os.path.isfile('hwe.table'):
            l = ['']
            for i in xrange(3,9+1):
                l += ['10^-%i' %(i)]
            s_header = '\t'.join(l)+'\n'
            fd = open('hwe.table','w')
            fd.write(s_header)
            fd.close()

        fd = open('hwe.table','a')
        fd.write(s)
        fd.close()

        os.system('sort hwe.table -k1,1 -u -o hwe.table')

        return


    def table_lmiss(self,bfile,):

        s = ''
        cmd = '''awk 'NR>1{print 1-$5}' %s.SNPQC.lmiss''' %(
            bfile,
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

        self.table_header_cumfreq('lmiss',lowerreallimit,binsize,)

        s = self.cumfreq2string(l_cumfreq,)
        line = '%s\t%s\n' %(bfile,s,)
        fd = open('lmiss.table','a')
        fd.write(line)
        fd.close()

        os.system('sort lmiss.table -k1,1 -u -o lmiss.table')

        return


    def table_header_cumfreq(self,prefix,lowerreallimit,binsize,):

        if not os.path.isfile('%s.table' %(prefix)):
            l = ['#']
            for i in xrange(int((1-lowerreallimit)/binsize)+1):
                l += ['%.3f' %(lowerreallimit+i*binsize)]
            s = '\t'.join(l)+'\n'
            fd = open('%s.table' %(prefix),'w')
            fd.write(s)
            fd.close()

        return
    

    def table_genome(self,bfile,):

        if (
            ## no need to append if project is uganda with only one population
            self.project == 'uganda'
            and
            os.path.isfile('genome.table')
            and
            os.path.isfile('genome_outliers.table')
            ):
            return

        print '\n\ngenome'

        ##
        ## related samples at different thresholds
        ##

        n_samples = int(os.popen('cat %s.fam | wc -l' %(bfile)).read())
        n_samples -= int(os.popen('cat %s.sampleQC.samples | wc -l' %(bfile)).read())

        l = [bfile]
        l_header = ['#']
        numbins = 4
        l_pi_hat_max = []
        for i in xrange(numbins):
            pi_hat_max = (i+1)*((.20-0.)/numbins)
            l_pi_hat_max += [pi_hat_max]
        l_pi_hat_max += [0.5,0.9,]
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
        s = '\t'.join(l)+'\n'
        s_header = '\t'.join(l_header)+'\n'

        if not os.path.isfile('genome.table'):
            fd = open('genome.table','w')
            fd.write(s_header)
            fd.close()

        fd = open('genome.table','a')
        fd.write(s)
        fd.close()

        os.system('sort genome.table -k1,1 -u -o genome.table')

        ##
        ## highly related samples
        ##
        cmd = '''awk 'NR>1{if($10>0.65) print $1,$3,$10}' %s.prehardy.genome''' %(bfile)
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
        cmd = '''awk 'NR>1{if($10>0.90) print $1,$3,$10}' %s.prehardy.genome > %s.genome.high''' %(bfile,bfile,)
        self.execmd(cmd)
        cmd = 'join -1 2 -2 1 -o 0,2.3,1.6 %s.imiss %s.genome.high > genome.high.table' %(bfile,bfile,)
        self.execmd(cmd)
        cmd = 'join -1 2 -2 2 -o 0,2.3,1.6 %s.imiss %s.genome.high > genome.high.table' %(bfile,bfile,)
        self.execmd(cmd)
        os.remove('%s.genome.high' %(bfile))

        return


    def table_imiss(self,bfile,):

        print
        print 'imiss'
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

        self.table_header_cumfreq('imiss',lowerreallimit,binsize,)

        s = self.cumfreq2string(l_cumfreq,)
        line = '%s\t%s\n' %(bfile,s,)
        fd = open('imiss.table','a')
        fd.write(line)
        fd.close()

        os.system('sort imiss.table -k1,1 -u -o imiss.table')

        return l_imiss


    def cumfreq2string(self,l_cumfreq,):

        s = '\t'.join(
            [str(l_cumfreq[3])]+[str(int(v)) for v in l_cumfreq[0]+l_cumfreq[3]]
            )

        return s


    def table_het(self,bfile,):

        cmd = 'cat %s.imiss.samples ' %(bfile)
        cmd += "| awk '{print $1}' > %s.imiss.1column.samples" %(bfile)
        self.execmd(cmd)

        cmd = 'fgrep -w -v -f %s.imiss.1column.samples %s.het' %(bfile,bfile,)
        cmd += '''| awk 'NR>1{print ($5-$3)/$5}' '''
        if self.verbose == True:
            print cmd
        l_het = os.popen(cmd).readlines()
        os.remove('%s.imiss.1column.samples' %(bfile))

        ## convert to float        
        l_het = [float(het) for het in l_het]
        ## calculate statistics
        mean = average = numpy.mean(l_het)
        stddev = numpy.std(l_het)
        print 'average', average, 'stddev', stddev

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
        fd = open('het.table','a')
        fd.write(line)
        fd.close()

        os.system('sort het.table -k1,1 -u -o het.table')

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
##            [['SNPQC.lmiss','frq',],self.scatter_lmiss_frq,],
            [['mds',],self.scatter_mds,],
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


    def scatter_mds_excl1000g(self,bfile,):

        '''this function needs to be rewritten.
it's ugly and I will not understand it 1 year form now.'''

        if not os.path.isdir('mds'):
            os.mkdir('mds')

        if self.project == 'uganda':

            ##
            ## tmp scatter plot of mds
            ##
            d = {}
            for ethnicity in ['Murundi','Muganda','Munyarwanda',]:
                fd = open('samples/%s.samples' %(ethnicity),'r')
                s = fd.read()
                fd.close()
                l = s.split('\n')
                for s in l:
                    d[s] = ethnicity
                fd = open('%s.mds' %(ethnicity),'w')
                fd.close()
            fd = open('%s.mds' %(bfile),'r')
            lines = fd.readlines()
            fd.close()
            for line in lines[1:]:
                FID = line.split()[0][-10:]
                if FID in d.keys():
                    prefix = d[FID]
                else:
                    prefix = 'other'
                fd = open('%s.mds' %(prefix),'a')
                fd.write(line)
                fd.close()
    ##            print FID, os.popen('grep %s samples/*' %(FID)).read()
    ##        stop
            
            line_plot = 'plot '
    ##        line_plot += '"%s.mds" u 4:5 lc 0 ps 2 pt 7 t ""\n' %(bfile)
            line_plot += '"Murundi.mds" u 4:5 lc 0 ps 2 pt 7 t "",' %()
            line_plot += '"Muganda.mds" u 4:5 lc 1 ps 2 pt 7 t "",' %()
            line_plot += '"Munyarwanda.mds" u 4:5 lc 2 ps 2 pt 7 t "",' %()
            line_plot += '"other.mds" u 4:5 lc 3 ps 2 pt 7 t ""\n' %()
            cmd = 'cat %s.mds | wc -l' %(bfile)
            n_samples = int(os.popen(cmd).read().strip())-1
            n_SNPs = int(os.popen('cat %s.posthardy.prune.in | wc -l' %(bfile)).read())
            gnuplot.scatter_plot_2d(
                '%s.mds' %(bfile),
    ##            s_plot = '"< paste %s.het %s.imiss" u (($5-$3)/$5):(1-$12)' %(bfile,bfile,),
    ##            s_plot = s_plot,
                line_plot = line_plot,
                column1 = 4, column2 = 5,
                xlabel = 'C1',
                ylabel = 'C2',
                title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                    bfile.replace('_','\\\\_'),n_samples,n_SNPs,
                    ),
                )

        return


    def scatter_mds(self,bfile,):

        '''this function needs to be rewritten.
it's ugly and I will not understand it 1 year form now.'''

        if not os.path.isdir('mds'):
            os.mkdir('mds')

        if self.project == 'uganda':

            self.scatter_mds_excl1000g(bfile,)

        elif self.project == 'agv':

            ##
            ## scatter plot of mds
            ##
            d = {}
            d_pops = {}
            for fn in os.listdir('samples'):
                ethnicity = fn[:-8]
    ##            if len(ethnicity) != 3 and ethnicity not in ['Murundi','Muganda','Munyarwanda',]:
    ##                continue
                fd = open('samples/%s' %(fn),'r')
                s = fd.read()
                fd.close()
                l = s.split('\n')
                d_pops[ethnicity] = []
                for s in l:
                    if s == 'NA':
                        continue
                    if ethnicity == '':
                        continue
                    if s in d.keys():
                        if ethnicity == 'Unknown':
                            continue
                        print s, d[s], ethnicity
                        stop
                    d[s] = ethnicity
                    d_pops[ethnicity] += [s]
                fd = open('mds/%s_%s.mds' %(bfile,ethnicity),'w')
                fd.close()

            fd = open('%s.mds' %(bfile),'r')
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
                prefix = d[FID]
    ##            print prefix
                fd = open('mds/%s_%s.mds' %(bfile,prefix),'a')
                fd.write(line)
                fd.close()
    ##            print FID, os.popen('grep %s samples/*' %(FID)).read()
    ##        stop

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

            if self.project == 'uganda':
                l_tribes = [
                    'Muganda',
                    'Munyarwanda',
                    'RwandeseUgandan',
                    'Murundi',
                    'Munyankole',
                    'Mukiga',
                    'Other',
                    'Mutanzania',
                    'Musoga',
                    'Mutooro',
                    'Mufumbira',
                    'Nyanjiro(Tanzania)',
                    'Unknown',
                    ]
            else:
                l_tribes = [bfile,]
            
            i = 0
            d_colors = {}
            for tribe in l_tribes:
                d_colors[tribe] = [i*20,i*20,i*20,]
                i += 1

            l_pops = list(set(d_pops.keys())-set(l_tribes))
            l_pops.sort()
            l_pops += l_tribes
            
            line_plot = 'plot [:0.25]'
    ##        line_plot = 'plot [:0][-0.01:0.02]'
            i_color = 0
            d_plot = {'uganda':'','1000g':'',}
            for i_pop in xrange(len(l_pops)):
                pop = l_pops[i_pop]
                fn_mds = 'mds/%s_%s.mds' %(bfile,pop)
                if not os.path.isfile(fn_mds):
                    print 'skipping', pop
                    continue
                if os.path.getsize(fn_mds) == 0:
                    print 'skipping', pop
                    continue
                if pop in d_colors.keys():
                    color = d_colors[pop]
                    ps = 2
                    pt = 7
                else:
                    print pop
##                    color = l_colors[i_color%len(l_colors)]
                    color = l_colors[i_color]
                    i_color += 1
                    ps = 2
                    pt = 7
                s = '"%s" u 4:5 lc rgb "#%s" ps %i pt %s t "%s",' %(
                    fn_mds,
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
            n_samples = int(os.popen('cat %s.mds | wc -l' %(bfile)).read())
            n_samples -= 1
            cmd = 'cat %s.%s.prune.in | wc -l' %(bfile,self.fn1000g,)
            n_SNPs = int(os.popen(cmd).read())
            gnuplot.scatter_plot_2d(
                '%s.mds' %(bfile),
    ##            s_plot = '"< paste %s.het %s.imiss" u (($5-$3)/$5):(1-$12)' %(bfile,bfile,),
    ##            s_plot = s_plot,
                line_plot = line_plot,
                column1 = 4, column2 = 5,
                xlabel = 'C1',
                ylabel = 'C2',
                title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                    bfile.replace('_','\\\\_'),n_samples,n_SNPs,
                    ),
                )
##            if self.project == 'agv':
##                stop_color_by_chip_cf_dg11_bp7

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
        cmd += "awk '{print $1, ($5-$3)/$5}' > %s_join2.txt" %(bfile,)
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
            cmd += "| awk '{print $1, ($5-$3)/$5}'"
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
            cmd = '''awk 'NR>1{print ($5-$3)/$5}' %s.het''' %(bfile)
            l_het = os.popen(cmd).readlines()
            l_het = [float(het) for het in l_het]
            n_samples = n = len(l_het)
            average = numpy.mean(l_het)
            stddev = numpy.std(l_het)

        n_SNPs = int(os.popen('cat %s.autosomes.SNPs | wc -l' %(bfile)).read())

##        print average, 3*stddev

        l_arrows = [
##            'set arrow from graph(0),%f to graph(1),%f nohead lc 0\n' %(
            'set arrow from 0,%f to 1,%f nohead lc 0 lt 2 lw 2\n' %(
                self.threshold_imiss,self.threshold_imiss,
                ),
            'set arrow from 0,%f to 1,%f nohead lc 0 lt 2 lw 2\n' %(
                0.99,0.99,
                ),
            ]
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
                #.98:'orange',
                }
            for threshold in d_colors.keys():
                cmd = '''
    awk 'NR>1 {if(1-$6<%f) print $1}' %s.imiss > %s.imiss.filtered
    grep -v -f %s.imiss.filtered %s.het | awk 'NR>1{print ($5-$3)/$5}'
    ''' %(
        threshold,bfile,bfile,bfile,bfile,
        )

                l_het = os.popen(cmd).readlines()
                l_het = [float(het) for het in l_het]
                average_w_exclusion = numpy.mean(l_het)
                stddev_w_exclusion = numpy.std(l_het)
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

            os.remove('%s.imiss.filtered' %(bfile))

        if bool_with_stddev == True:
            line_plot = 'plot [%f:%f][:1]' %(
                min(min(l_het),average-3*stddev,)-0.001,
                max(max(l_het),average+3*stddev,)+0.001,
                )
        else:
            line_plot = 'plot [:][:1]'
        line_plot += '"%s.imiss.het.joined" ' %(bfile)
        line_plot += 'u (($10-$8)/$10):(1-$6) lc 0 ps 2 pt 7 t ""'

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
                bfile.replace('_','\_'),n_samples,n_SNPs,
                ),
            lines_extra = l_arrows,
            )

        os.remove('%s.sex.opposite.dat' %(bfile))
        os.remove('%s.sex.unknown.dat' %(bfile))
        os.remove('%s.sex.dat' %(bfile))
        os.remove('%s.imiss.het.joined' %(bfile))

        return

    

    def plink_execution(self,bfile,):

        l_plink_cmds = self.l_plink_cmds(bfile,)

        ##
        ## execute plink commands
        ##
        for plink_cmd_full in l_plink_cmds:

            plink_cmd = plink_cmd_full.split()[0].replace('--','')

            ##
            ## check that file output doesn't already exist
            ##
            bool_continue_out, out_prefix = self.check_output_existence(
                bfile,plink_cmd,plink_cmd_full,
                )

            ##
            ## check that file input exists (LSF dependency does not work...)
            ##
            bool_continue_in, in_prefix = self.check_input_existence(
                bfile,plink_cmd,plink_cmd_full,
                )

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

            ## run all in one go (not finished...)
            ## .. will not finish this as results should be checked during the process
            if self.bool_run_all == True and plink_cmd != 'genome':
                cmd += self.cmd_rerun(bfile,)

            ##
            ## terminate command
            ##
            if self.bool_verbose == True:
                cmd += self.mail(bfile,plink_cmd,out_prefix,'finished',)

            if bool_continue_out == True:
                print 'out exists', plink_cmd
                continue
            if bool_continue_in == True:
                print 'in does not exist', plink_cmd
                continue

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
                'missing','check-sex','het','flip',
                ]:
                print cmd
                stop1
            if '%' in cmd.replace('%02g',''):
                print cmd
                stop2
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


    def cmd_rerun(self,bfile,):

        cmd = '\n'
        cmd += 'if [ ! -s %s.SNPQC.bed ]\nthen\n' %(bfile)
        cmd += 'sleep 300\n\n'
        cmd += 'python %s/QC.py ' %(
            os.path.dirname(sys.argv[0]),
            )
        ## add all variables if default values are not used...
        cmd += '--project %s --bfile %s' %(self.project,bfile,)
        cmd += '\nfi'
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
        if bool_file_out == True:
            bool_continue = True

        return bool_continue, out_prefix


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
            ## 'flip',
            ]:
            if plink_cmd == 'bmerge' and option == 'extract':
                fp_in = '%s.posthardy.prune.in' %(bfile,)
                l_fp_in += [fp_in]
                pass
            elif '--%s' %(option) in plink_cmd_full:
                l = plink_cmd_full.split()
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
                    print plink_cmd, fp_in, 'does not exist'
                    pass
                break
            continue
        if bool_input == False:
            bool_continue = True
            pass

        return bool_continue, bfile_prefix


    def append_plink(self,bfile,plink_cmd_full,out_prefix,):

##        ## http://pngu.mgh.harvard.edu/~purcell/plink/flow.shtml
####        cmd += ' --keep-before-remove --keep %s.non1000g.samples' %(bfile)
##        cmd += '--keep %s.non1000g.samples \\\n' %(bfile)
        
        plink_cmd = plink_cmd_full.split()[0].replace('--','')

        ## initiate plink
        cmd = 'plink \\\n'

        ## --bfile
        if '--bfile' not in plink_cmd_full:
            cmd += '--bfile %s \\\n' %(bfile,)
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
            cmd += '-o stdout/plink.%s.%s.out \\\n' %(bfile,plink_cmd,) ## tmp
            cmd += '-e stderr/plink.%s.%s.err \\\n' %(bfile,plink_cmd,) ## tmp
        else:
            cmd += '-o /dev/null \\\n'
        ## specify JOBID
        if not JOBID:
            JOBID = '%s.%s' %(bfile,plink_cmd,)
        cmd += '-J"%s" \\\n' %(JOBID,)
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
            'check-sex':900,
            'cluster':900,
            'het':1000,
            'bmerge':1700,
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
            cmd = 'cat %s.fam | wc -l' %(self.fp1000g)
            if self.bool_verbose == True:
                print cmd
            n_samples += int(os.popen(cmd).read())

        ## assign appropriate amount of memory (GB)
        if plink_cmd in d_mem.keys():
            mem = int(
                max(2,math.ceil((d_mem[plink_cmd]/1000.)*n_samples/1000.))
                )
        else:
            mem = int(max(2,math.ceil(0.8*n_samples/1000.)))
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
        s_cluster += '--mds-plot 4 \\\n'
        s_cluster += '--exclude %s.ldregions.SNPs \\\n' %(bfile,)
        s_cluster += '--extract %s.posthardy.prune.in \\\n' %(bfile,)
        s_cluster += '--remove %s.SNPQC.samples \\\n' %(bfile,)

        ##
        ## sequential for Manj
        ##
        l_plink_cmds += ['--het --remove %s.imiss.samples --out %s.sequential' %(bfile,bfile,)]
        l_plink_cmds += ['--check-sex --extract %s.X.SNPs --remove %s.imiss.samples --out %s.sequential' %(bfile,bfile,bfile,)]

        ##
        ##
        ##
        l_plink_cmds += [
            ## http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#flip
            '--flip %s.flip \\\n--make-bed --out %s.flipped' %(bfile,bfile,),
            ]

        ##
        ## sample QC
        ##

        l_plink_cmds += [
            ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#missing
            '--missing',
            ## http://pngu.mgh.harvard.edu/~purcell/plink/ibdibs.shtml#inbreeding
            '--het',
            ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#sexcheck
            '--check-sex --extract %s.X.SNPs' %(bfile),
            ]
        ## http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed
        s = '--make-bed \\\n'
        s += '--flip %s.flip \\\n' %(bfile)
        s += '--remove %s.sampleQC.samples \\\n' %(bfile)
        s += '--out %s.sampleQC \\\n' %(bfile)
        l_plink_cmds += [s]

        ##
        ## SNP QC, pre Hardy
        ##

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#freq
        s = '--freq \\\n'
        s += '--remove %s.sampleQC.samples \\\n' %(bfile)
        s += '--out %s.sampleQC \\\n' %(bfile)
        s += '--exclude %s.lmiss.SNPs \\\n' %(bfile)
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#missing
        s = '--missing \\\n'
        s +='--out %s.SNPQC \\\n' %(bfile,)
        s += '--remove %s.sampleQC.samples \\\n' %(bfile,)
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

        ## http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#hardy
        s = '--hardy \\\n'
        s += '--exclude %s.lmiss.SNPs \\\n' %(bfile,)
        s += '--remove %s.genome.prehardy.samples \\\n' %(bfile,)
        l_plink_cmds += [s]

        ##
        ## SNP QC, post Hardy
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
        cmd_cluster = s_cluster+'--read-genome %s.posthardy.genome \\\n' %(
            bfile,
            )
        cmd_cluster += '--bfile %s \\\n' %(bfile,)
        l_plink_cmds += [cmd_cluster]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed
        s = '--make-bed \\\n'
        s += '--flip %s.flip \\\n' %(bfile)
        s += '--remove %s.SNPQC.samples \\\n' %(bfile)
        s += '--out %s.SNPQC \\\n' %(bfile)
        s += '--exclude %s.SNPQC.SNPs \\\n' %(bfile,)
        l_plink_cmds += [s]

        ##
        ## MDS with 1000G
        ##

        ## http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#merge
        s = '--bmerge \\\n'
        s += '%s.%s.flipped.bed \\\n' %(bfile,self.fn1000g,)
        s += '%s.%s.flipped.bim \\\n' %(bfile,self.fn1000g,)
        s += '%s.%s.flipped.fam \\\n' %(bfile,self.fn1000g,)
        s += '--make-bed \\\n'
        s += '--out %s.%s \\\n' %(bfile,self.fn1000g,)
        s += '--bfile %s.flipped \\\n' %(bfile,)
        s += '--remove %s.SNPQC.samples \\\n' %(bfile)
        s += '--extract %s.%s.prune.in \\\n' %(bfile,self.fn1000g,)
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/ibdibs.shtml#genome
        s = self.write_genome_cmd(
            bfile, self.fn1000g, 'SNPQC', '%s.%s' %(bfile,self.fn1000g,),
            )
        l_plink_cmds += [s]

        ## http://pngu.mgh.harvard.edu/~purcell/plink/strat.shtml#cluster
        s = s_cluster+'--read-genome %s.%s.genome \\\n' %(bfile,self.fn1000g,)
        s += '--bfile %s.%s \\\n' %(bfile,self.fn1000g,)
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

        ##
        ## exclude samples from previous steps
        ##
        if plink_cmd == 'bmerge':
            l_cmds = self.bmerge_before(bfile,)

        elif plink_cmd == 'genome':
            l_cmds = self.genome_before(bfile,in_prefix,out_prefix,)

        elif plink_cmd == 'flip':
            cmd = '''awk '{if($5=="-") print $1}' %s | sort > %s.flip''' %(
                self.strand,bfile,
                )
            l_cmds += [cmd]

        cmds = '\n\n'.join(l_cmds)
        
        return cmds


    def extra_commands_after(self,bfile,plink_cmd,in_prefix,out_prefix,):

        l_cmds = []

        if plink_cmd == 'missing':
            l_cmds = self.missing_after(bfile)

        elif plink_cmd == 'het':
            l_cmds = self.het_after(bfile)
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
            cmd = "awk '{if ($9 < %.1e) print $2}' %s.hwe > %s.hwe.SNPs" %(
                self.hwe_min,
                bfile,bfile,
                )
            l_cmds += [cmd]
            cmd = 'cat %s.lmiss.SNPs %s.hwe.SNPs > %s.SNPQC.SNPs' %(
                bfile,bfile,bfile,
                )
            l_cmds += [cmd]

##        elif plink_cmd == 'cluster':
##            l_cmds += ['fi\n']

        elif plink_cmd == 'bmerge':
            l_cmds += ['rm %s.%s.flipped.*' %(bfile, self.fn1000g,)] ## --bfile

##        elif plink_cmd == 'flip':
##            l_cmds += ['rm %s.flip' %(bfile,)]

        elif plink_cmd == 'make-bed' and out_prefix == '%s.SNPQC' %(bfile):
            l_cmds += self.makebed_after(bfile,in_prefix,out_prefix,)

        cmds = '\n\n'.join(l_cmds)

        return cmds


    def makebed_after(self,bfile,in_prefix,out_prefix,):

        l_cmds = []
        cmd = self.append_LSF(bfile,'eigensoft',verbose=False,)
        cmd += '%s \\\n' %(self.eigensoft)
        ## genotype file in any format (see ../CONVERTF/README)
        cmd += '-i %s.bed \\\n' %(os.path.join(os.getcwd(),out_prefix))
        ## snp file in any format (see ../CONVERTF/README)
        cmd += '-a %s.bim \\\n' %(os.path.join(os.getcwd(),out_prefix))
        ## indiv file in any format (see ../CONVERTF/README)
        cmd += '-b %s.fam \\\n' %(os.path.join(os.getcwd(),out_prefix))
        ## (Default is 10) number of principal components to output
        cmd += '-k 10 \\\n'
        ## output file of principal components.  Individuals removed
        ## as outliers will have all values set to 0.0 in this file.
        cmd += '-o %s.pca \\\n' %(out_prefix.replace('.','_'))
        ## prefix of output plot files of top 2 principal components.
        ## (labeling individuals according to labels in indiv file)
        cmd += '-p %s.plot \\\n' %(out_prefix.replace('.','_'))
        ## output file of all eigenvalues
        cmd += '-e %s.eval \\\n' %(out_prefix.replace('.','_'))
        ## output logfile        
        cmd += '-l %s.eigensoft.log \\\n' %(out_prefix.replace('.','_'))
        ## (Default is 5) maximum number of outlier removal iterations.
        ## To turn off outlier removal, set -m 0.
        cmd += '-m 5 \\\n'
        ## (Default is 10) number of principal components along which
        ## to remove outliers during each outlier removal iteration.
        cmd += '-t 10 \\\n'
        ## (Default is 6.0) number of standard deviations which an
        ## individual must exceed, along one of topk top principal
        ## components, in order to be removed as an outlier.
        cmd += '-s 6 \\\n'

##        l_cmds += [cmd]
        fd = open('%s.pca.sh' %(bfile),'w')
        fd.write(cmd)
        fd.close()

        cmd = 'echo "'
        cmd += 'Hi Deepti. This is an automated message sent to you, because my script ran to completion. '
        cmd += 'Would you mind running this EIGENSOFT script for me: '
        cmd += '%s/%s.pca.sh ' %(os.getcwd(),bfile,)
        cmd += 'Thank you! Tommy'
        cmd += '"'
        cmd += '| mail -s "EIGENSOFT" '
        cmd += 'dg11@sanger.ac.uk\n'
        l_cmds += [cmd]

        cmd = 'echo "'
        cmd += 'Hi Deepti. This is an automated message sent to you, because my script ran to completion. '
        cmd += 'Would you mind running this EIGENSOFT script for me: '
        cmd += '%s/%s.pca.sh ' %(os.getcwd(),bfile,)
        cmd += 'Thank you! Tommy'
        cmd += '"'
        cmd += '| mail -s "EIGENSOFT" '
        cmd += 'tommy.carstensen@gmail.com\n'
        l_cmds += [cmd]

        return l_cmds


##    def cluster_after(self,bfile,in_prefix,out_prefix,):
##
##        l_cmds = []
##        cmd = '\n;\n'
##
##        cmd += '"'
##        l_cmds += [cmd]
##
##        return l_cmds


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

        cmd = ''

        ## sort fam file
        cmd = "cat %s.fam | sort -k1,1 -k2,2" %(in_prefix,)
        cmd += " | awk '{print $1,$2}' > %s.fam.sorted" %(out_prefix,)
        l_cmds += [cmd]

        ## sort sample removal file
        cmd = 'cat %s | sort -k1,1 -k2,2 ' %(fn_samples_remove,)
        cmd += '> %s.sorted' %(fn_samples_remove,)
        l_cmds += [cmd]
        
        ## exclude samples
        cmd = "comm -23 %s.fam.sorted %s.sorted" %(out_prefix,fn_samples_remove,)
        cmd += " > %s.filtered.fam" %(out_prefix,)
        l_cmds += [cmd]

        ## parse FID and IID columns
        cmd = 'cat %s.filtered.fam' %(out_prefix,)
        cmd += ' | '
        cmd += " awk '{print $1,$2}' "
        ## and pipe them to individual fam files
        cmd += ' | split -d -a 2 -l %i - fam/%s.fam.' %(
            self.len_lists,out_prefix,
            )
        l_cmds += [cmd]

        ## clean up
        cmd = 'rm %s.sorted\n' %(fn_samples_remove)
        cmd += 'rm %s.fam.sorted\n' %(out_prefix)
        cmd += 'rm %s.filtered.fam\n' %(out_prefix)
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
        cmd = '\n'
        cmd += 'for i in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
        cmd += 'for j in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
        cmd += 'if [ $j -gt $i ]\nthen continue\nfi\n\n'
        cmd += 'if [ -s genome/%s.$i.$j.genome ]\nthen continue\nfi\n\n' %(
            out_prefix,
            )
        ## initiate command
        cmd += 'cmd="\n'
        l_cmds += [cmd]

        return l_cmds


    def genome_after(self,bfile,in_prefix,out_prefix,):

        l_cmds = []

        ## terminate command
        cmd = '\n"\n\n'
        ## evaluate command
        cmd += '%s$cmd\n' %(
            self.append_LSF(bfile,'genome',JOBID='genome.${i}.${j}',),
            )
        ## end loops
        cmd += '\ndone\ndone\n'
        l_cmds += [cmd]

        ## count lines
        cmd = 'nlines=0\n'
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
        cmd += 'if [ $nlines -eq $(($n*($n-1)/2)) ]\nthen\n'

        ## initiate .genome file with header
        s = '                   FID1                   IID1                   FID2                   IID2 RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO'
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
            ## disallow duplicates (i.e. IBD>0.90)
            if self.project == 'uganda':
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
            if self.project == 'uganda':
                ## concatenate sample exclusion lists (i.e. IBD>0.90)
                cmd += 'cat %s.sampleQC.samples %s.genome.%.2f.samples' %(
                    bfile,bfile,.9,
                    )
                cmd += ' > %s.SNPQC.samples\n\n' %(bfile,)
            elif self.project == 'agv':
                ## concatenate sample exclusion lists (i.e. IBD>low threshold)
                cmd += 'cat %s.sampleQC.samples %s.genome.%.2f.samples' %(
                    bfile,bfile,self.pi_hat_max,
                    )
                cmd += ' > %s.SNPQC.samples\n\n' %(bfile,)
            else:
                stop_unknown_project
            ## concatenate sample exclusion lists (related)
            cmd += 'cat %s.sampleQC.samples %s.genome.%.2f.samples' %(
                bfile,bfile,self.pi_hat_max,
                )
            cmd += ' > %s.genome.prehardy.samples\n' %(bfile,)

        if self.bool_run_all == True:
            cmd += self.cmd_rerun(bfile,)

        cmd += self.mail(bfile,'genome',out_prefix,'finished',)

        cmd += '\nfi\n'
        l_cmds += [cmd]

        return l_cmds


    def het_after(self,bfile,):

        l_cmds = []
        
        ## a) calculate heterozygosity mean and stddev and use for exclusion of samples
        ## make list of samples going to be excluded by sample call rate threshold
        cmd = "awk '{print $1}' %s.imiss.samples > %s.awk.tmp" %(bfile,bfile,)
        l_cmds += [cmd]
        ## use var=$(cmd) if bash...
        ## use set var=`cmd` if tcsh...
        ## set variable s equal to stdout of command in parenthesis
        cmd = 's=$( '
        ## exclude samples below call rate threshold before calc of stddev
        ## I probably should use sort and join instead of grep... if IDs not unique...
        cmd += 'fgrep -w -v -f %s.awk.tmp %s.het' %(bfile,bfile,)
        cmd += ' | '
        ## calculate mean and stddev for remaining samples
        cmd += "awk 'NR>1 "
        cmd += " {het=($5-$3)/$5; sum+=het; sumsq+=het**2} "
        cmd += " END {print sum/(NR-1), sqrt(sumsq/(NR-1) - (sum/(NR-1))**2)}' "
        cmd += ')'
        cmd += ';mean=$(echo $s | cut -d " " -f1)'
        cmd += ';stddev=$(echo $s | cut -d " " -f2)'
        cmd += ';rm %s.awk.tmp' %(bfile)
        l_cmds += [cmd]
        ## b) write list of samples outside mean+-3stddev
        cmd = 'awk -v mean=$mean -v stddev=$stddev '
        cmd += " 'NR>1 {if( "
        cmd += '(($5-$3)/$5)<mean-%i*stddev || (($5-$3)/$5)>mean+%i*stddev ' %(
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
        cmd += "awk 'NR>1 {if($6>%f) print $1,$2}' %s.imiss " %(
            1-self.threshold_imiss,
            bfile,
            )
        cmd += ' > %s.imiss.samples\n' %(bfile)
        cmd += 'rm %s.lmiss\n' %(bfile)
        cmd += self.concatenate_sampleQC_remove_lists(bfile,)
        cmd += 'fi\n'
        l_cmds += [cmd]

        ##
        ## lmiss.SNPs
        ##
        cmd += '\n'
        cmd = 'if [ -s %s.SNPQC.lmiss -a ! -s %s.lmiss.SNPs ]\n' %(bfile,bfile,)
        cmd += 'then\n'
        cmd += 'cat %s.SNPQC.lmiss' %(bfile)
        cmd += ' | '
        cmd += "awk 'NR>1 {if ((1-$5)<%.2f) print $2}' " %(
            self.threshold_lmiss,
            )
        cmd += ' > '
        cmd += '%s.lmiss.SNPs\n' %(bfile)
        cmd += 'fi\n'
        l_cmds += [cmd]

        return l_cmds


    def bmerge_before(self,bfile,):

        l_cmds = []

        ## agreed with Manj at meeting 15oct2012
        ## to find SNPs common between the chip and 1000g
        ## by position instead of rsID
        l_cmds += ['##\n## 1a) rename 1000g\n##']
        l_cmds += self.rename_1000g(bfile,)

        l_cmds += ['##\n## 1b) flip 1000g\n##']

        ## 2a) 1000g.flip
        ## 2a1) sort chip bim
##        cmd = 'cat %s.bim | sort -k2,2 > %s.bim.sorted' %(bfile,bfile,)
        cmd = 'cat %s.flipped.bim | sort -k2,2 > %s.bim.sorted' %(bfile,bfile,)
        l_cmds += [cmd]

        ## 2a2) sort renamed 1000g bim
        cmd = 'cat %s.%s_renamed.bim | sort -k2,2 > %s.%s_renamed.bim.sorted' %(
            bfile,self.fn1000g,
            bfile,self.fn1000g,
            )
        l_cmds += [cmd]

        ## 2a3) join bims by name (join should be by chr+pos!!!)
        cmd = 'join -o 0 1.5 1.6 2.5 2.6 '
        cmd += '-1 2 -2 2 %s.bim.sorted %s.%s_renamed.bim.sorted '%(
            bfile, bfile, self.fn1000g,
            )
        cmd += ' > %s.%s.bim.joined' %(bfile,self.fn1000g,)
        l_cmds += [cmd]
        l_cmds += ['rm %s.bim.sorted' %(bfile)]

        ## 2a4) only allow alleles to be flipped,
        ## if they are flipped differently
        ## *before* flipping
        cmd = "cat %s.%s.bim.joined" %(bfile,self.fn1000g,)
        cmd += " | awk '{if (!("
        cmd += '(($2==$4 || $2==0 || $4==0) && ($3==$5 || $3==0 || $5==0))'
        cmd += ' || '
        cmd += '(($2==$5 || $2==0 || $5==0) && ($3==$4 || $3==0 || $4==0))'
        cmd += ")) print $1}'"
        cmd += " > %s.%s.different.flip" %(bfile,self.fn1000g,)
        l_cmds += [cmd]

        ## clean up
##        l_cmds += ['rm %s.%s.bim.joined' %(bfile,self.fn1000g,)]

        ## 2b) flip renamed 1000g
##        cmd = 'if [ ! -s %s.flipped.bed ]\n' %(self.fn1000g)
##        cmd += 'then\n\n'
        cmd = 'echo "---flip 1000g---"\n'
        cmd += 'plink \\\n'
        cmd += '--bfile %s.%s_renamed \\\n' %(bfile,self.fn1000g,)
        cmd += '--flip %s.%s.different.flip \\\n' %(bfile,self.fn1000g,)
        cmd += '--make-bed \\\n'
        cmd += '--out %s.%s.flipped \\\n' %(
            bfile,self.fn1000g,
            )
        cmd += '--extract %s.%s.prune.in \\\n' %(bfile,self.fn1000g,)
        cmd += '--noweb \\\n'
        cmd += '--allow-no-sex \\\n'
        cmd += '--nonfounders \\\n'
##        cmd += '\nfi\n'
        l_cmds += [cmd]

        ## clean up
        l_cmds += ['rm %s.%s_renamed*' %(bfile, self.fn1000g,)] ## --bfile
        l_cmds += ['rm %s.%s.different.flip' %(bfile,self.fn1000g,)] ## --flip

        l_cmds += ['##\n## 2) merge flipped files\n##\n']

        return l_cmds


    def rename_1000g(self,bfile,):

        l_cmds = []

        col_chr = 1
        col_name = 2
        col_pos = 4

        ##
        ## 1) create prune.in.chrpos.sorted
        ##

        ## 1a) sort .prune.in
        cmd = 'cat %s.posthardy.prune.in | sort > %s.prune.in.sorted' %(bfile,bfile,)
        l_cmds += [cmd]

        ## 1b) sort .bim
        cmd = 'cat %s.bim | sort -k%i,%i > %s.bim.sorted' %(
            bfile,col_name,col_name,bfile,
            )
        l_cmds += [cmd]

        ## 1c) join .prune.in and .bim by name (i.e. find common names)
        cmd = 'join -1 1 -2 %i' %(col_name)
        cmd += ' -o 2.1,2.2,2.3,2.4,2.5,2.6'
        cmd += ' %s.prune.in.sorted %s.bim.sorted' %(bfile,bfile,)
        ## sort by chr and pos
        cmd += ' | sort -k%i,%i -k%i,%i' %(col_chr,col_chr,col_pos,col_pos,)
        ## print chr+pos column
        cmd += ''' | awk '{print $0,$%i":"$%i}' ''' %(col_chr,col_pos,)
        ## generate join file 2
        cmd += ' > %s.prune.in.chrpos.sorted' %(bfile,)
        l_cmds += [cmd]
        l_cmds += ['rm %s.prune.in.sorted %s.bim.sorted' %(bfile,bfile,)]

        ##
        ## 2) create 1000g.bim.chrpos.sorted
        ##

        ## 2a) sort 1000g.bim by chr and pos
        cmd = 'cat %s.bim | ' %(self.fp1000g)
        cmd += 'sort -k%i,%i -k%i,%i' %(col_chr,col_chr,col_pos,col_pos,)
        ## print chr+pos and name column
        cmd += ''' | awk '{print $0,$%i":"$%i}' ''' %(col_chr,col_pos,)
        ## generate join file 1
        cmd += ' > %s.%s.bim.chrpos.sorted' %(
            bfile,self.fn1000g,
            )
        l_cmds += [cmd]

        ##
        ## 3) join prune.in.chrpos.sorted and 1000g.bim.chrpos.sorted by chrpos
        ##
        
        ## join bims by position
        cmd = 'join -1 7 -2 7 -o 1.2,2.2,1.5,1.6,2.5,2.6 '
        cmd += '%s.%s.bim.chrpos.sorted %s.prune.in.chrpos.sorted' %(
            bfile,self.fn1000g,
            bfile,
            )
        cmd += ' > %s.%s.bim.prune.joined' %(
            bfile, self.fn1000g,
            )
        l_cmds += [cmd]

        ##
        ## 4)
        ##
        
        ## 4a) write 2 column dic
        cmd = "cat %s.%s.bim.prune.joined" %(bfile,self.fn1000g,)
        cmd += " | awk '{print $1,$2}'"
        cmd += ' > %s.%s.prune.dic' %(
            bfile, self.fn1000g,
            )
        l_cmds += [cmd]
        l_cmds += ['rm %s.%s.bim.chrpos.sorted %s.prune.in.chrpos.sorted' %(
            bfile,self.fn1000g,
            bfile,
            )]

        ## 4b) do not extract mismatch SNPs (e.g. kgp11141322 AG/GA/TC/CT vs GT/TG/CA/AC)
        cmd = "cat %s.%s.bim.prune.joined" %(bfile,self.fn1000g,)
        cmd += " | awk '{if (!("
        cmd += '(($3==0 || $5==0 || $3==$5) && ($4==0 || $6==0 || $4==$6))'
        cmd += ' || '
        cmd += '(($3==0 || $6==0 || $3==$6) && ($4==0 || $5==0 || $4==$5))'
        cmd += ")) print $1,$2}'"
        cmd += ' > %s.%s.mismatch.SNPs.2columns' %(bfile,self.fn1000g,)
        l_cmds += [cmd]
        l_cmds += ['rm %s.%s.bim.prune.joined' %(bfile,self.fn1000g,)]        
        ##
        cmd = "cat %s.%s.mismatch.SNPs.2columns" %(bfile,self.fn1000g,)
        cmd += " | awk '{print $1}' "
        cmd += "> %s.%s.mismatch.SNPs.1column" %(bfile,self.fn1000g,)
        l_cmds += [cmd]
        ##
        cmd = 'cat %s.%s.mismatch.SNPs.2columns' %(bfile,self.fn1000g,)
        cmd += " | awk '{print $2}' "
        cmd += '>> %s.%s.mismatch.SNPs.1column' %(bfile,self.fn1000g,)
        l_cmds += [cmd]
        ##
        l_cmds += ['rm %s.%s.mismatch.SNPs.2columns' %(bfile, self.fn1000g,),]

        ## 4c) only extract SNPs common to both beds
        cmd = 'cat %s.%s.prune.dic | ' %(bfile, self.fn1000g,)
        cmd += "awk '{print $2}' | "
        cmd += 'sort > %s.%s.prune.in.sorted' %(bfile, self.fn1000g,)
        l_cmds += [cmd]

        ## 4d) exclude mismatch SNPs from extraction list
        cmd = 'sort %s.%s.mismatch.SNPs.1column ' %(bfile,self.fn1000g,)
        cmd += '> %s.%s.mismatch.SNPs.sorted' %(bfile,self.fn1000g,)
        l_cmds += [cmd]
        l_cmds += ['rm %s.%s.mismatch.SNPs.1column' %(bfile,self.fn1000g,),] ## sort
        cmd = 'comm -23 '
        cmd += '%s.%s.prune.in.sorted ' %(bfile,self.fn1000g,)
        cmd += '%s.%s.mismatch.SNPs.sorted ' %(bfile,self.fn1000g,)
        cmd += '> %s.%s.prune.in' %(bfile,self.fn1000g,)
        l_cmds += [cmd]
        cmd = 'rm %s.%s.prune.in.sorted' %(bfile,self.fn1000g,)
        l_cmds += [cmd] ## comm1
        cmd = 'rm %s.%s.mismatch.SNPs.sorted' %(bfile,self.fn1000g,)
        l_cmds += [cmd] ## comm2

        ##
        ## 5) rename the 1000g SNPs
        ##
        l_cmds += ['## --update-name']
        l_cmd = ['plink']
        l_cmd += ['--bfile %s' %(self.fp1000g,)]
        l_cmd += ['--make-bed --out %s.%s_renamed' %(
            bfile, self.fn1000g,
            )]
        ## the file quad2octo.dic is generated by the function extract_and_translate
        l_cmd += ['--update-map %s.%s.prune.dic --update-name' %(
            bfile, self.fn1000g,
            )]
        l_cmd += ['--extract %s.%s.prune.in' %(bfile,self.fn1000g,)]
        l_cmd += ['--noweb --nonfounders']
        cmd = ' \\\n'.join(l_cmd)
        l_cmds += [cmd]

        ##
        ## clean up (need --extract file for --bmerge)
        ##
        ## --update-map
        l_cmds += ['rm %s.%s.prune.dic' %(bfile, self.fn1000g,),]

        return l_cmds


    def concatenate_sampleQC_remove_lists(self,bfile,):

        cmd = 'if [ '
        ## check file out
        cmd += '! -s %s.sampleQC.samples ' %(bfile)
        ## check file in
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
        cmd = 'cat %s.X.SNPs | wc -l' %(bfile)
        n_SNPs = int(os.popen(cmd).read())

        x_step=0.001

        if bool_with_stddev == True:

            cmd = "cat %s.imiss.samples | awk '{print $1}' > %s.imiss.1column.samples" %(
                bfile,bfile,
                )
            self.execmd(cmd)

            cmd = 'fgrep -w -v -f %s.imiss.1column.samples %s.het' %(bfile,bfile,)
            cmd += '''| awk 'NR>1{print ($5-$3)/$5}' '''
            l_het = os.popen(cmd).readlines()
            os.remove('%s.imiss.1column.samples' %(bfile))
            ## exclude imiss outliers before stddev calc...!!!
            l_het = [float(het) for het in l_het]
            n_samples = len(l_het)
            average = numpy.mean(l_het)
            stddev = numpy.std(l_het)
            print average,stddev

            l_arrows = [
                'set arrow from %f,0 to %f,graph(1) nohead lc 0\n' %(
                    average-3*stddev,average-3*stddev,
                    ),
                'set arrow from %f,0 to %f,graph(1) nohead lc 0\n' %(
                    average+3*stddev,average+3*stddev,
                    ),
                ]

            s_plot = 'plot [%f:%f]"%s.het" ' %(
                min(min(l_het),average-3*stddev,)-x_step,
                max(max(l_het),average+3*stddev,)+x_step,
                bfile,
                )
            s_plot += 'u (hist(($5-$3)/$5,width)):(1.0) smooth freq w boxes '
            s_plot += 'lc rgb"red" notitle\n'

            pass

        else:

            s_plot = 'plot "%s.het" ' %(bfile)
            s_plot += 'u (hist(($5-$3)/$5,width)):(1.0) smooth freq w boxes '
            s_plot += 'lc rgb"red" notitle\n'
            l_arrows = None

        gnuplot.histogram2(
            '%s.het' %(bfile),
##            x_min=het_min,x_max=het_max,#tic_step=0.05,
            x_step = 0.001,
            s_plot = s_plot,
            xlabel='heterozygosity',
            title='%s (n_{samples}=%i, n_{SNPs}=%i)' %(
                bfile.replace('_','\\\\_'),n_samples,n_SNPs,
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

        cmd = 'cat %s.prehardy.genome | wc -l' %(bfile)
        if self.verbose == True:
            print cmd
        n = int(os.popen(cmd).read())-1
        ## http://simple.wikipedia.org/wiki/Quadratic_equation
        n_samples1 = int(((--.5)+math.sqrt((-.5)**2-4*.5*-n))/(2*.5))

        n_samples2 = int(os.popen('cat %s.fam | wc -l' %(bfile)).read())
        n_samples2 -= int(os.popen('cat %s.sampleQC.samples | wc -l' %(bfile)).read())

        if n_samples1 != n_samples2:
            print n_samples1
            print n_samples2
            stop
        else:
            n_samples = n_samples1 = n_samples2

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
        else:
            x_step = .002
    
        gnuplot.histogram2(
            '%s.SNPQC.lmiss' %(bfile),
            prefix_out='%s.lmiss' %(bfile),
            column='(1-$5)',
            x_step=x_step,
            x_min=self.threshold_lmiss,
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
        ## download and unzip strand file
        ##
        self.download_and_unzip_strand_file(bfile,)

        ##
        ## write SNPs with mismatched chromosome and/or position
        ## between strand/multiple files and bim file
        ## N.B. this step will also exclude SNPs unique to bim file
        ## N.B. this step will also exclude unplaced SNPs
        ##
        self.strand_bim_mismatch_position(bfile,)
        ##
        ## write SNPs in strand.miss file
        ##
        self.strand_miss(bfile,)
        ##
        ## write SNPs with duplicate chromosomal positions
        ##
        self.bim_duplicates(bfile,)

##        gnuplot.venn3(
##            f1='%s.position.SNPs' %(bfile),
##            f2='%s.miss.SNPs' %(bfile),
##            f3='%s.duplicates.SNPs' %(bfile),
##            suffix='%s.SNPs' %(bfile),
##            text1='chromosome+position mismatch',
##            text2='in the miss file',
##            text3='duplicate chromosome+position',
##            )

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
            'flip':['bed',],
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


    def bim_duplicates(self,bfile,):

        ## this function assumes that the only multiplets are duplets

        if (
            os.path.isfile('%s.duplicates.SNPs' %(bfile,))
            and
            os.path.isfile('%s.autosomes.SNPs' %(bfile,))
            and
            os.path.isfile('%s.X.SNPs' %(bfile,))
##            and
##            os.path.getsize('%s.duplicates.SNPs' %(bfile,)) > 0
            ):
            return

        cmd = '''cat %s.bim | sort -k1,1 -k4,4 | awk '{print $2,$1,$4}' ''' %(bfile)
        cmd += '> %s.bim.sorted' %(bfile,)
        self.execmd(cmd)

        fd_in = open('%s.bim.sorted' %(bfile),'r')
        fd_out = open('%s.duplicates.SNPs' %(bfile,),'w')
        fd_out_inclusion = open('%s.duplicates.inclusion.SNPs' %(bfile,),'w')
        line_prev = 'x x x'
        while True:
            line = fd_in.readline()
            if not line:
                break
            l = line.split()
            l_prev = line_prev.split()
            if (
                ## chr
                l_prev[1] == l[1]
                and
                ## pos
                l_prev[2] == l[2]
                and
                ## not unplaced
                l[1] != '0'
                ):
                if l[0][:2] == 'rs':
                    fd_out.write('%s\n' %(l_prev[0]))
                    fd_out_inclusion.write('%s\n' %(l[0]))
                elif l_prev[0][:2] == 'rs':
                    fd_out.write('%s\n' %(l[0]))
                    fd_out_inclusion.write('%s\n' %(l_prev[0]))
                else:
                    fd_out.write('%s\n' %(l[0]))
                    fd_out_inclusion.write('%s\n' %(l_prev[0]))
            line_prev = line
        fd_in.close()
        fd_out.close()
        fd_out_inclusion.close()

        os.remove('%s.bim.sorted' %(bfile))

        return


    def strand_miss(self,bfile,):

        if (
            os.path.isfile('%s.miss.SNPs' %(bfile,))
            and
            os.path.isfile('%s.autosomes.SNPs' %(bfile,))
            and
            os.path.isfile('%s.X.SNPs' %(bfile,))
            ):
            return

        ## cat HumanOmni2.5-8v1_A-b37.miss | awk 'NR>1{print $3}' | wc -l
        ## =
        ## fgrep -w -f %s.miss.SNPs omni2.5-8_20120809_gwa_uganda_gtu.bim | wc -l
        cmd = "cat %s | awk 'NR>1{print $3}' > %s.miss.SNPs" %(self.miss,bfile,)
        self.execmd(cmd)

        return


    def strand_bim_mismatch_position(self,bfile,):

        if (
            os.path.isfile('%s.position.SNPs' %(bfile,))
            and
            os.path.isfile('%s.autosomes.SNPs' %(bfile,))
            and
            os.path.isfile('%s.X.SNPs' %(bfile,))
            ):
            return

        cmd = '''cat %s | sort -k1,1 | ''' %(self.strand,)
        cmd += '''awk '{sub(/X/,23,$2);sub(/Y/,24,$2)'''
        cmd += ''';sub(/XY/,25,$2);sub(/MT/,26,$2)'''
        cmd += ''';print $1,$2,$3}' > %s.%s.sorted''' %(bfile,self.strand,)
        self.execmd(cmd)
        cmd = '''cat %s.bim | sort -k2,2 | awk '{print $2,$1,$4}' > %s.bim.sorted''' %(
            bfile,bfile,
            )
        self.execmd(cmd)
        cmd = 'comm -23 %s.bim.sorted %s.%s.sorted | ' %(
            bfile,bfile,self.strand,
            )
        cmd += "awk '{print $1}' > %s.position.SNPs" %(bfile,)
        self.execmd(cmd)

        os.remove('%s.%s.sorted' %(bfile,self.strand,))
        os.remove('%s.bim.sorted' %(self.bfile,))

        return


    def download_and_unzip_strand_file(self,bfile,):

        n_SNPs = int(os.popen('cat %s.bim | wc -l' %(bfile)).read())
        ## it's really stupid to choose the strand based on the number of SNPs
        ## this is a temporary solution...
        if n_SNPs == 2449626: ## quad
            fn_zip = 'HumanOmni2.5M-b37-strand-v2'
            strand = 'HumanOmni2.5M-b37-v2.strand'
            miss = 'Strand-HumanOmni2.5-b37.miss'
        elif n_SNPs == 2379855:
            fn_zip = 'HumanOmni2.5-8v1_A-b37-strand.zip'
            strand = 'HumanOmni2.5-8v1_A-b37.strand'
            miss = 'HumanOmni2.5-8v1_A-b37.miss'
        elif n_SNPs == 2379514: ## octo
            fn_zip = 'HumanOmni2.5-8v1_A-b37-strand.zip'
            strand = 'HumanOmni2.5-8v1_A-b37.strand'
            miss = 'HumanOmni2.5-8v1_A-b37.miss'
        elif n_SNPs == 2362478: ## merged quad+octo
            fn_zip = 'HumanOmni2.5-8v1_A-b37-strand.zip'
            strand = 'HumanOmni2.5-8v1_A-b37.strand'
            miss = 'HumanOmni2.5-8v1_A-b37.miss'
        ## assume octo (add command line option...)
        else:
            print n_SNPs
            stop
        self.strand = strand
        self.miss = miss

        self.strand1000g = 'HumanOmni2.5M-b37-v2.strand'
        self.miss1000g = 'Strand-HumanOmni2.5-b37.miss'
        fn_zip_1000g = 'HumanOmni2.5M-b37-strand-v2.zip' ## quad

        for strand,miss,fn in [
            [self.strand,self.miss,fn_zip,],
            [self.strand1000g,self.miss1000g,fn_zip_1000g,],
            ]:
            if not os.path.isfile('%s' %(strand)):
                cmd = 'wget http://www.well.ox.ac.uk/~wrayner/strand/%s' %(fn)
                os.system(cmd)
                os.system('unzip %s' %(fn))
                os.remove('%s' %(fn))

                os.system('dos2unix %s' %(strand))
                os.system('dos2unix %s' %(miss))

            if not os.path.isfile(strand):
                stop1
        if not os.path.isfile(miss):
            stop2

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

        for out_prefix, condition in [
            ['X','$1==23',],
            ['autosomes','$1>=1 && $1<=22',],
            ]:

            fn_out = '%s.%s.SNPs' %(bfile,out_prefix,)

            ## continue if file already exists
            if os.path.isfile(fn_out):
                continue

            ## concatenate exclusion lists
            if not os.path.isfile('%s.preQC.SNPs.sorted' %(bfile)):
                cmd = 'cat %s.position.SNPs %s.miss.SNPs %s.duplicates.SNPs' %(
                    bfile,bfile,bfile,
                    )
                cmd += ' | sort -u > %s.preQC.SNPs.sorted' %(bfile,)
                self.execmd(cmd)
                

            ## generate inclusion list
            cmd = "cat %s.bim | awk '{if(%s) print $2}' | sort > %s.sorted" %(
                bfile,condition,fn_out,
                )
            self.execmd(cmd)

            ## subtract exclusion list from inclusion list
            cmd = 'comm -23 %s.sorted %s.preQC.SNPs.sorted > %s' %(
                fn_out,bfile,fn_out,
                )
            self.execmd(cmd)

            os.remove('%s.sorted' %(fn_out))
            os.remove('%s.preQC.SNPs.sorted' %(bfile))

        ## clean up after ourselves; i.e. delete exclusion lists
##        for suffix in ['position','miss','duplicates','preQC',]:
##            if os.path.isfile('%s.%s.SNPs' %(bfile,suffix,)):
##                os.remove('%s.%s.SNPs' %(bfile,suffix,))

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
            metavar='FLOAT',default=0.10,
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
            '--hwe','--hwe_min', '--threshold_hwe_min',
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
        ## write options in use to file
        ##
        d_options = vars(opts)
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
        ## other options
        ##

        ## length of lists used for pairwise IBD estimation run in parallel
        self.len_lists = 400

        self.bool_run_all = False
        if '--run-all' in sys.argv:
            self.bool_run_all = True
        self.bool_run_all = True

        self.bool_sequential = False

        ## perform various sanity checks of file lengths etc.
        self.bool_check_file_io = True

        ##
        ## paths
        ##
        
        self.dn1000g = '/nfs/african_diversity/data/1KG/genotypes/release20110527/working_data/omni25_b37_bed_autosomes'
        self.fn1000g = '1kg_r20110527_omni25_b37_autosomes'
        self.fp1000g = os.path.join(self.dn1000g,self.fn1000g,)

        ## dg11 2012oct09 06:32
        self.cluster_exclude = '/nfs/t149_influenza_exomes/working/analysis2/ldregions.txt'

        self.eigensoft = '/nfs/team149/Software/usr/share/eigensoft/EIG4.2/bin/smartpca.perl'

        ##
        ## create dirs
        ##
        for dn in ['stdout','stderr','fam','genome',]:
            if not os.path.isdir(dn):
                os.mkdir(dn)
                pass
            continue

        return

if __name__ == '__main__':
    instance_main = main()
    instance_main.main(instance_main.bfile)
