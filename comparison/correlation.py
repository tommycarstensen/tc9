#!/software/bin/python3

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012-2013

## built-ins
import os, math, sys, argparse, re
sys.path.append('/nfs/users/nfs_t/tc9/github/sandbox')
import gnuplot

l_chroms = [str(chrom) for chrom in xrange(1,22+1)]+['X']

class main():


    def main(self,):

        if self.extension1 in ['bed','ped',]:
            self.PLINK_transpose(self.fp1)
            self.fn1 = os.path.basename(self.fp1)
            self.extension1 = 'tped'
            self.fp1 = self.fn1

        if self.extension1 == 'tped' and self.extension2 == 'bgl':
            self.compare_tped_with_bgl_or_gen_or_gprobs(
                self.l_fp1,self.l_fp2,self.affix,header2=True,)
        elif self.extension1 == 'tped' and self.extension2 == 'gprobs':
            self.compare_tped_with_bgl_or_gen_or_gprobs(
                self.l_fp1,self.l_fp2,self.affix,header2=True,)
        elif self.extension1 == 'tped' and self.extension2 == 'gen':
            self.compare_tped_with_bgl_or_gen_or_gprobs(
                self.l_fp1,self.l_fp2,self.affix,header2=False,)
        elif self.extension1 == 'gprobs' and self.extension2 == 'imputed':
            self.compare_tped_with_bgl_or_gen_or_gprobs(
                self.l_fp1,self.l_fp2,self.affix,header2=False,)
        elif self.format1 == 'IMPUTE2' and self.format2 == 'IMPUTE2':
            self.compare_tped_with_bgl_or_gen_or_gprobs(
                self.l_fp1,self.l_fp2,self.affix,header2=False,)
        else:
            print(self.extension1, self.extension2)
            print('did I write code for handling this???')
            sys.exit(0)

    ####    ped = 'omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped'
    ####    sepjoin = 'sep'
    ####
    ####    PLINK_transpose(ped,'ped',)
    ##
    ##    instance = GATK_pipeline.main()
    ##    d_chromosome_lengths = instance.parse_chromosome_ranges()
    ##    d_centromere_ranges = parse_centromere_ranges()
    ##
    ##    IMPUTE2_tped(ped,d_chromosome_lengths,sepjoin,d_centromere_ranges,)
    ##
    ##    bgl_ped(ped,d_chromosome_lengths,sepjoin,d_centromere_ranges,)

        return


    def parse_line_bgl(self,line_bgl):

        l_bgl = line_bgl.strip().split()
        marker_bgl = l_bgl[0]
        l = marker_bgl.split(':')
        chrom_bgl = l[0]
        pos_bgl = int(l[1])

        return l_bgl, chrom_bgl, pos_bgl


    def compare_tped_with_bgl_or_gen_or_gprobs(
        self,l_fp1,l_fp2,affix,header2=True,):

        ## todo: replace all instances of bgl in this function with something generic...

        ##
        ## open fp2 file stream
        ##
        index_fp2 = 0
        chrom2 = str(index_fp2+1)
        fp2 = l_fp2[index_fp2]
        ##
        ## open file2
        ##
        fd2 = open(fp2,'r')

        if self.extension1 in ['tped']:
            header1 = False
        elif self.extension1 == 'gprobs':
            header1 = True
        else:
            print(self.extension1)
            stop

        ##
        ## parse file2 header/samples
        ##
        if self.extension2 in ['bgl','gprobs',] and self.fp_samples2 == None:
            line2 = fd2.readline()
            l = line2.split()
            l_samples2 = [l[i] for i in xrange(3,len(l),3)]
            parse = 'gprobs'
            i_alleleA = 1
            i_alleleB = 2
        elif self.extension2 in ['gen','imputed',] and self.fp_samples2 != None:
            fd = open(self.fp_samples2,'r')
            s = fd.read()
            fd.close()
            l_samples2 = s.strip().split('\n')
            parse = 'gen'
            i_alleleA = 3
            i_alleleB = 4
        else:
            print(self.extension2)
            print(self.fp_samples2)
            stop

        ##
        ## find relevant line numbers in ped
        ## instead of doing lookup in samples from gprobs all the time
        ##
        if self.extension1 == 'tped':
            l_samples1 = self.parse_samples(self.fp_samples1)
            func1 = self.parse_line_tped
        elif self.fp_samples1:
            print(self.fp_samples1)
            stoptmp
        elif self.extension1 == 'gprobs':
            print(l_fp1[0])
            fd = open(l_fp1[0],'r')
            s = fd.readline()
            fd.close()
            l_samples1 = s.strip().split()[3:-1:3]
            func1 = self.parse_line_gprobssss
            stop
            ## parse_line_gprobssss and parse_line_tped must rreturn the same output
            ## i.e. same number of variables and a list with genotype probabilities
        else:
            print(self.extension1)
            stop

        d_samples = {}
        if self.fp_sampledic:
            fd = open(self.fp_sampledic)
            for line in fd.read().split('\n'):
                if line.strip() == '': continue
                sample1 = line.split()[1]
                sample2 = line.split()[0]
                d_samples[sample1] = sample2
            fd.close()
        else:
            for sample in l_samples1:
                d_samples[sample] = sample

        ##
        ## match sequence of file2 and tfam sample IDs
        ##
        l_indexes1 = []
        l_indexes2 = []
        if parse == 'gprobs':
            for sample1 in l_samples1:
                if not sample1 in d_samples.keys():
                    continue
                sample2 = d_samples[sample1]
                if not sample2 in l_samples2:
                    continue
                l_indexes2 += [
                    ## marker alleleA alleleB
                    3+3*l_samples2.index(sample2)]
##            l_indexes2 = [
##                ## marker alleleA alleleB
##                3+3*l_samples2.index(sample1) for sample1 in l_samples1]
        else:
            for index1 in xrange(len(l_samples1)):
                sample1 = l_samples1[index1]
                if not sample1 in l_samples2:
                    continue
                ## --- rsID position alleleA alleleB
                index2 = 5+3*l_samples2.index(sample1)
                l_indexes2 += [index2]
        for index1 in xrange(len(l_samples1)):
            sample1 = l_samples1[index1]
            if not sample1 in d_samples.keys():
                continue
            sample2 = d_samples[sample1]
            if not sample2 in l_samples2:
                continue
            l_indexes1 += [4+2*index1]

        n_samples = len(l_indexes1)# = len(l_indexes2)

        if n_samples == 0:
            print(l_samples1)
            print(l_samples2)
            print('n_samples=0')
            print(d_samples)
            sys.exit()

        ##
        ## set counts and lists
        ##
        ## sets
        cnt1 = 0
        cnt2 = 0
        ## complement set
        cnt_1_not_2 = 0
        cnt_2_not_1 = 0
        ## intersection set
        cnt_1_and_2 = 0
        ## maybe use unix utility join instead of memory hungry python lists...
        ## list of rsIDs sets; complement and intersection
        l_1_not_2 = []
        l_1_and_2 = []
        l_corr = []
        l_conc = []
        ## discordant count
        cnt_discordant_all_SNPs = 0
        cnt_missing_all_SNPs = 0
        index_fp1 = 0
        fp1 = l_fp1[index_fp1]
        ##
        ## open file1
        ##
        fd1 = open(fp1,'r')
        fd_discordant = open('%s.discordant' %(affix),'w')

        if header1 == True:
            fd1.readline()

        ##
        ## parse first line
        ##
        bool_EOF1 = False
        bool_EOF2 = False
        bool_read1 = True
        bool_read2 = True
        while True:

            ##
            ## parse lines
            ##
            if bool_EOF1 == False:
                if bool_read1 == True:
                    line1 = fd1.readline()
                    if line1=='':
                        if index_fp1+1 == len(l_fp1):
                            bool_EOF1 = True
                        else:
                            index_fp1 += 1
                            fp1 = l_fp1[index_fp1]
                            fd1 = open(fp1,'r')
                            if header1 == True:
                                line1 = fd1.readline()
                            line1 = fd1.readline()
                            l1, chrom1, pos1 = func1(line1)
                            cnt1 += 1
                    else:
                        l1, chrom1, pos1 = func1(line1)
                        cnt1 += 1
                        if cnt1 % 10000 == 0: print('line/variant/SNP', cnt1, chrom1, pos1)
            if bool_EOF2 == False:
                if bool_read2 == True:
                    line2 = fd2.readline()
                    if line2=='':
                        ## last file
                        if index_fp2+1 == len(l_fp2):
                            bool_EOF2 = True
                        ## open next file
                        else:
                            index_fp2 += 1
                            fp2 = l_fp2[index_fp2]
                            fd2 = open(fp2,'r')
                            if header2 == True:
                                line2 = fd2.readline()
                            line2 = fd2.readline()
                            if parse == 'gprobs':
                                l2, chrom2, pos2 = self.parse_line_bgl(line2)
                            else:
                                chrom2 = str(index_fp2+1)
                                l2, pos2 = self.parse_line_gen(line2)
                            cnt2 += 1
                    else:
                        if parse == 'gprobs':
                            l2, chrom2, pos2 = self.parse_line_bgl(line2)
                        else:
                            chrom2 = str(index_fp2+1)
                            l2, pos2 = self.parse_line_gen(line2)
                        cnt2 += 1

            ##
            ## continue or break or pass?
            ##
            if bool_EOF1 == True:
                bool_read2 = True
                if bool_EOF2 == True:
                    break
                else:
                    cnt_2_not_1 += 1
                    continue
            elif bool_EOF2 == True:
                bool_read1 = True
                cnt_1_not_2 += 1
                l_1_not_2 += [l1[1]]
                continue

            (
                bool_continue, bool_read1, bool_read2,
                cnt_1_not_2, cnt_2_not_1, cnt_1_and_2,
                l_1_not_2, l_1_and_2,) = self.check_positions(
                    chrom1,chrom2,pos1,pos2,
                    cnt_1_not_2, cnt_2_not_1, cnt_1_and_2,
                    l_1_not_2,l_1_and_2,
                    l1,
                    )
            if bool_continue == True:
                continue

            alleleA2 = l2[i_alleleA]
            alleleB2 = l2[i_alleleB]

            ## deletion
            if len(alleleA2) > 1:
                continue
            ## insertion
            if len(alleleB2) > 1:
                continue

            cnt_discordant, cnt_missing, r2 = self.count_discordant(
                l1, l2, l_indexes1, l_indexes2,
                fd_discordant,
                n_samples,
                alleleA2, alleleB2,
                chrom1, pos1,
                )
            ## maybe move slow file I/O outside of loop...
            l_corr += ['%s %s\n' %(l1[1],r2,)]
            l_conc += ['%s %s\n' %(l1[1],cnt_discordant,)]
            ## append count for current SNP to cumulative count across all SNPs
            cnt_discordant_all_SNPs += cnt_discordant
            cnt_missing_all_SNPs += cnt_missing

            continue
        
        ##
        ## close files
        ##
        fd_discordant.close()
        fd1.close()
        fd2.close()

        ##
        ## correlation / concordance
        ##
        fd = open('%s.correlation' %(affix),'w')
        fd.writelines(l_corr)
        fd.close()
        
        fd = open('%s.concordance' %(affix),'w')
        fd.writelines(l_conc)
        fd.close()

        cmd = 'cat %s.correlation' %(affix)
        cmd += ''' | awk '{if($2!="None") {sum+=$2;n++}} END{print sum/n}' '''
        r2_avg = float(os.popen(cmd).read())

        ##
        ##
        ##
        s = ''
        s += 'fp1: %s\n' %(fp1)
        s += 'fp2: %s\n' %(fp2)
        s += 'fp1: %i variants/SNPs\n' %(cnt_1_not_2+cnt_1_and_2)
        s += 'fp2: %i variants/SNPs\n' %(cnt_2_not_1+cnt_1_and_2)
        s += 'intersection: %i variants/SNPs\n' %(cnt_1_and_2)
        s += 'n_samples: %i\n' %(n_samples)
        if cnt1 != cnt_1_not_2+cnt_1_and_2:
            print(cnt1, cnt_1_not_2+cnt_1_and_2)
            print(cnt_1_not_2, cnt_1_and_2)
            stoptmp1
        if cnt2 != cnt_2_not_1+cnt_1_and_2:
            print(cnt2, cnt_2_not_1+cnt_1_and_2)
            print(cnt_2_not_1, cnt_1_and_2)
            stoptmp2
        if chrom1 != '22' or chrom2 not in['22','Y',]:
            print(chrom1, chrom2, pos1, pos2)
            stoptmp3
        ## discordant
        s += 'cnt_discordant_all_SNPs: %i\n' %(cnt_discordant_all_SNPs)
        ## missing
        s += 'cnt_missing_all_SNPs: %i\n' %(cnt_missing_all_SNPs)
        s += 'concordance: %.4f\n' %(
            1-(float(cnt_discordant_all_SNPs)/float(n_samples*cnt_1_and_2)))
        s += 'correlation: %.4f\n' %(r2_avg)            
        print(s)
        fd = open('%s.stats' %(affix),'a')
        fd.write(s)
        fd.close()
        for fn_SNPs,l_SNPs in [
            ['%s.1_not_2.SNPs' %(affix),l_1_not_2,],
            ['%s.1_and_2.SNPs' %(affix),l_1_and_2,],
            ]:
            s_SNPs = '\n'.join(l_SNPs)+'\n'
            fd = open(fn_SNPs,'w')
            fd.write(s_SNPs)
            fd.close()

        return


    def parse_samples(self,fp_samples,):
        
        if fp_samples:
            if fp_samples[fp_samples.rindex('.'):] in ['.sample','.tfam',]:
                cmd = "cat %s | awk '{print $2}'" %(fp_samples)
                l_samples1 = os.popen(cmd).read().strip().split('\n')
        else:
            print(fp_samples)
            stop

        ## convert plate_well_sample format to sample format
        ## remove info about plate and well
##        keyword = re.compile(r'(\d\d\d\d\d\d_[A-H]\d\d_)(APP\d\d\d\d\d\d\d)')
        keyword = re.compile(r'(\d\d\d\d\d\d_[A-H]\d\d_)(.+\d\d\d\d\d\d\d)')
        l = []
        for sampleID in l_samples1:
            match = result = keyword.search(sampleID)
            if match:
                l += [match.group(2)]
            else:
                l += [sampleID]
        l_samples1 = l

        return l_samples1


    def parse_line_gen(self,line):

        ## slow split
        l = line.strip().split()
##        marker = l[0]
        pos = int(l[2])

        return l, pos


    def check_positions(
        self,
        chrom1,chrom2,pos1,pos2,
        cnt_1_not_2, cnt_2_not_1, cnt_1_and_2,
        l_1_not_2, l_1_and_2,
        l1,
        ):

        bool_read1 = False
        bool_read2 = False
        if chrom1 != chrom2:
            if l_chroms.index(chrom1) > l_chroms.index(chrom2):
                bool_read2 = True
                cnt_2_not_1 += 1
                bool_continue = True
            else:
                bool_read1 = True
                cnt_1_not_2 += 1
                l_1_not_2 += [l1[1]]
                bool_continue = True
        elif pos1 > pos2:
            bool_read2 = True
            cnt_2_not_1 += 1
            bool_continue = True
        elif pos2 > pos1:
            bool_read1 = True
            cnt_1_not_2 += 1
            l_1_not_2 += [l1[1]]
            bool_continue = True
        else:
            ## continue loop over both lines
            bool_read1 = True
            bool_read2 = True
            cnt_1_and_2 += 1
            l_1_and_2 += [l1[1]]
            bool_continue = False

        return (
            bool_continue, bool_read1, bool_read2,
            cnt_1_not_2, cnt_2_not_1, cnt_1_and_2,
            l_1_not_2, l_1_and_2,
            )


    def count_discordant(
        self, l_tped, l_bgl,
        l_indexes1, l_indexes2,
        fd_discordant,
        n_samples, ## derived from l_indexes_bgl/l_indexes2
        alleleA_bgl, alleleB_bgl, ## derived from l_bgl
        chrom_tped, pos_tped, ## derived from l_tped
        ):

        sum_xy = 0
        sum_xx = 0
        sum_yy = 0
        sum_x = 0
        sum_y = 0
        n = 0

        cnt_discordant = 0
        cnt_missing = 0
        bool_return = False
        for i in xrange(n_samples):
##            index_tfam = 4+2*i
            index_tfam = l_indexes1[i]
            alleleA_tped = l_tped[index_tfam]
            alleleB_tped = l_tped[index_tfam+1]
            index_bgl = l_indexes2[i]
            ## l_probs_bgl is the BEAGLE genotype triplet
            l_probs_bgl = l_bgl[index_bgl:index_bgl+3]
            x = None
            y = None
            if alleleA_tped == '0' and alleleB_tped == '0':
                cnt_missing += 1
                continue
            ## insertion/deletion
            elif alleleA_bgl == '0' and alleleB_bgl == '1':
                continue
            elif alleleA_tped == alleleA_bgl:
                if alleleB_tped == alleleA_bgl:
                    x = 0.
                elif alleleB_tped == alleleB_bgl:
                    x = 1
                else:
                    bool_return = True
                    break
                    print(alleleA_tped, alleleB_tped, alleleA_bgl, alleleB_bgl)
                    stop
            elif alleleA_tped == alleleB_bgl:
                if alleleB_tped == alleleA_bgl:
                    x = 1
                elif alleleB_tped == alleleB_bgl:
                    x = 2
                else:
                    bool_return = True
                    break
                    print(alleleA_ped, alleleB_ped, alleleA_gprobs, alleleB_gprobs)
                    stop
            else:
                if (
                    alleleA_tped in [alleleA_bgl, alleleB_bgl]
                    and
                    alleleB_tped in [alleleA_bgl, alleleB_bgl]
                    ):
                    print(alleleA_tped, alleleB_tped, alleleA_bgl, alleleB_bgl, pos_tped)
                    stop
                bool_return = True
                break
##            if True:
##                pass ## tmp!!!
            if False: pass
            ## homozygous, AA
            elif float(l_probs_bgl[0]) > 0.9:
##                stop1
##                y = 0
                if (
                    alleleA_bgl == alleleA_tped
                    and
                    alleleA_bgl == alleleB_tped
                    ):
                    if x != 0: print(x, stop)
                    x = 0
                    pass
                else:
                    ## homozygous
                    if alleleA_tped == alleleB_tped:
                        if x != 2: print(x, stop)
                        x = 2
                    ## heterozygous
                    else:
                        if x != 1: print(x, stop)
                        x = 1
                    cnt_discordant += 1
            ## heterozygous, Aa or aA
            elif float(l_probs_bgl[1]) > 0.9:
##                stop2
##                y = 1
                if any([
                    all([
                        alleleA_bgl == alleleA_tped,
                        alleleB_bgl == alleleB_tped,
                        ]),
                    all([
                        alleleB_bgl == alleleA_tped,
                        alleleA_bgl == alleleB_tped,
                        ]),
                    ]):
                    if x != 1: print(x, stop)
                    x = 1
                    pass
                else:
                    ## homozygous 1
                    if (
                        alleleA_tped == alleleA_bgl
                        and
                        alleleB_tped == alleleA_bgl
                        ):
                        if x != 0: print(x, stop)
                        x = 0
                    ## homozygous 2
                    else:
                        if x != 2: print(x, stop)
                        x = 2
                    cnt_discordant += 1
            ## homozygous, aa
            elif float(l_probs_bgl[2]) > 0.9:
##                stop3
##                y = 2
                if (
                    alleleB_bgl == alleleA_tped
                    and
                    alleleB_bgl == alleleB_tped
                    ):
                    if x != 2: print(x, stop)
                    x = 2
                    pass
                else:
                    ## homozygous
                    if alleleA_tped == alleleB_tped:
                        if x != 0: print(x, stop)
                        x = 0
                    ## heterozygous
                    else:
                        if x != 1: print(x, stop)
                        x = 1
                    cnt_discordant += 1
##            ## pmax<.9
##            else:
##                stop4
##                cnt_missing += 1

            if x != None:
                ## x is chip "dosage", y is sequence/imputation dosage
                ## use dosage instead
                y = 0
##                y += 0*float(l_probs_bgl[0])
                y += 1*float(l_probs_bgl[1])
                y += 2*float(l_probs_bgl[2])
                sum_xy += x*y
                sum_xx += x*x
                sum_yy += y*y
                sum_x += x
                sum_y += y
                n += 1.

            ## continue loop over samples
            continue

        if bool_return == True:
            cnt_discordant = 0
            cnt_missing = 0
            r2 = None
            return cnt_discordant, cnt_missing, r2

        fd_discordant.write('%s %s:%s %s %s\n' %(
            l_tped[1],chrom_tped,pos_tped,cnt_discordant,cnt_missing,))

        ## apply a 90% sample call rate threshold
        if n > .9*n_samples:
##        if n > 0: ## tmp!!!
            nom = (sum_xy-sum_x*sum_y/n)
            den_sq = (sum_xx-sum_x**2/n)*(sum_yy-sum_y**2/n)
            ## problem if n=94, xx=2 x=2 yy=0.00274104 y=0.5076
            ## den_sq = -1.21430643318e-17
##            den = math.sqrt(round(den_sq,16))
            den = math.sqrt(abs(den_sq))
            ## likely monomorphic
            if den == 0:
##                ## monomorphic
##                if sum_xx == 2*sum_x or sum_yy == 2*sum_y:
##                    r2 = None
##                else:
##                    print sum_xx, sum_x, sum_yy, sum_y
##                    stop
                r2 = None
            else:
                r = nom/den
                r2 = r**2
        else:
            r2 = None

##        if r2 != None:
####            for i in xrange(n_samples):
####                index_tfam = l_indexes1[i]
####                alleleA_tped = l_tped[index_tfam]
####                alleleB_tped = l_tped[index_tfam+1]
####                index_bgl = l_indexes2[i]
####                l_probs_bgl = l_bgl[index_bgl:index_bgl+3]
####                print alleleA_tped, alleleB_tped, l_probs_bgl
##            print r2

##        if l_tped[1] == 'rs9958437' or r2 == 1/3. and l_tped[1] not in ['kgp15191706']:
##            ii = 0
##            for i in xrange(n_samples):
##                index_tfam = l_indexes1[i]
##                alleleA_tped = l_tped[index_tfam]
##                alleleB_tped = l_tped[index_tfam+1]
##                index_bgl = l_indexes2[i]
##                l_probs_bgl = l_bgl[index_bgl:index_bgl+3]
##                if alleleA_tped == '0' and alleleB_tped == '0':
##                    continue
##                if max(
##                    float(l_probs_bgl[0]),
##                    float(l_probs_bgl[1]),
##                    float(l_probs_bgl[2]),
##                    ) <= 0.9:
##                    print 'xxx', i
##                    continue
##                if lx[ii] != ly[ii]:
##                    print alleleA_tped, alleleB_tped, l_probs_bgl, lx[i], ly[i]
##                ii += 1
##            import scipy
##            from scipy import stats
##            r = scipy.stats.pearsonr(lx,ly)
##            print lx
##            print ly
##            print l_tped[1]
##            print r2
##            print nom, den
##            print sum_x, sum_xx, sum_y, sum_yy, n
##            print r[0]**2
##            print r
##            print r2
##            print n_samples, n
##            print len(lx)
##            stop

##        if n > 0.9*n_samples and r2 == None:
##            import scipy
##            from scipy import stats
##            import numpy
##            r = scipy.stats.pearsonr(lx,ly)
##            if r[0] == numpy.float('nan'):
##                print lx
##                print ly
##                print r
##                print sum_x, sum_xx, sum_y, sum_yy
##                print nom, den
##                print type(r[0])
##                print r[0]
##                stop2

        return cnt_discordant, cnt_missing, r2


    def parse_line_tped(self,line_tped,):

        ## slow split
        l_tped = line_tped.strip().split()
        chrom_tped = l_tped[0]
        pos_tped = int(l_tped[3])

        return l_tped, chrom_tped, pos_tped


    def execmd(self,cmd):

        print(cmd)
        os.system(cmd)

        return


    def parse_l_tped(self,l_tped):

##        chrom_tped = l_tped[0]
##        rsID_tped = l_tped[1]
        pos_tped = int(l_tped[3])

        return


    def parse_centromere_ranges(self,):

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


    def parse_line_gen_old(self,line):

        ## slow split
        l = line.strip().split()
        marker = l[0]
        pos = int(l[2])
        alleleA = l[3]
        alleleB = l[4]

        return l, marker, pos, alleleA, alleleB


    def parse_line_gprobs(self,line_gprobs):

        ## slow split
        l_gprobs = line_gprobs.strip().split()
        marker_gprobs = l_gprobs[0]
        pos_gprobs = int(marker_gprobs.split(':')[1])
        alleleA_gprobs = l_gprobs[1]
        alleleB_gprobs = l_gprobs[2]

        return l_gprobs, marker_gprobs, pos_gprobs, alleleA_gprobs, alleleB_gprobs


    def BEAGLE_ped(self,ped):

        stop_check_fd_write_at_end_is_generic_parse_some_prefix_from_cmd_line

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
            ## slow split
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
                        r2 = allelic_R2_ped(ped,row_map,l_indexes_samples,l_gprobs,alleleA_gprobs,alleleB_gprobs,)
    ##                    t2 = time.time()
    ##                    print t2-t1
    ##                    print (t2-t1)*1622321/3600.
    ##                    print r2, marker_map, marker_gprobs
    ##                    stop
                        l_r2_allelic += ['%s %s\n' %(marker_map,r2,)]
            ## identical marker found in map file
            else:
                r2 = allelic_R2_ped(ped,row_map,l_indexes_samples,l_gprobs,alleleA_gprobs,alleleB_gprobs,)
                print(r2, marker_map, marker_gprobs)
                l_r2_allelic += ['%s %s\n' %(marker_map,r2,)]
        fd_gprobs.close()

        fd = open('allelicR2_BEAGLEped.txt','w')
        fd.writelines(l_r2_allelic)
        fd.close()

        return


    def allelic_R2_tped(
        self,
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
                print(l_tped)
                stop
            if alleleA_ped == alleleA_gprobs:
                if alleleB_ped == alleleA_gprobs:
                    y = dosage_genotype = 0.
                    count_0 += 1
                elif alleleB_ped == alleleB_gprobs:
                    y = dosage_genotype = 1.
                    count_1 += 1
                else:
                    print(alleleA_ped, alleleB_ped, alleleA_gprobs, alleleB_gprobs)
                    stop
            elif alleleA_ped == alleleB_gprobs:
                if alleleB_ped == alleleA_gprobs:
                    y = dosage_genotype = 1.
                    count_1 += 1
                elif alleleB_ped == alleleB_gprobs:
                    y = dosage_genotype = 2.
                    count_2 += 1
                else:
                    print(alleleA_ped, alleleB_ped, alleleA_gprobs, alleleB_gprobs)
                    stop
            else:
                print('discordance',)
                print(alleleA_ped, alleleB_ped,)
                print(alleleA_gprobs, alleleB_gprobs, l_gprobs[0])
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


    def allelic_R2_ped(
        self,ped,row_map,l_indexes_samples,l_gprobs,
        alleleA_gprobs,alleleB_gprobs,):

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
                    print(alleleA_ped, alleleB_ped, alleleA_gprobs, alleleB_gprobs)
                    stop
            elif alleleA_ped == alleleB_gprobs:
                if alleleB_ped == alleleA_gprobs:
                    y = dosage_genotype = 1.
                elif alleleB_ped == alleleB_gprobs:
                    y = dosage_genotype = 2.
                else:
                    print(alleleA_ped, alleleB_ped, alleleA_gprobs, alleleB_gprobs)
                    stop
            else:
                print(alleleA_ped, alleleB_ped)
                print(alleleA_gprobs, alleleB_gprobs)
                print(sample_ped)
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
            print(n,)
        else:
            print(alleleA_ped)
            print(alleleB_ped)
            print(n)
            print(marker_map)
            print(marker_gprobs)
            print(pos_map)
            print(pos_gprobs)
            stop

        return r2


    def PLINK_transpose(self,fp,):

        dn = os.path.dirname(fp)
        fn = os.path.basename(fp)
        extension = fn[fn.rindex('.')+1:]
        basename = fn[:fn.rindex('.')]

        ## concantenate command arguments
        cmd = 'plink \\\n'
        cmd += '--noweb \\\n'
        if extension == 'ped':
            cmd += '--file %s \\\n' %(os.path.join(dn,basename))
        elif extension == 'bed':
            cmd += '--bfile %s \\\n' %(os.path.join(dn,basename))
        else:
            print(extension)
            print(fn)
            stop
        cmd += '--recode \\\n'
        cmd += '--transpose \\\n'
        cmd += '--out %s \\\n' %(basename)

        ## execute command
        if not os.path.isfile('%s.tped' %(basename)):
            execmd(cmd)

        ## check that output was generated
        if not os.path.isfile('%s.tped' %(basename)):
            stop

        return


    def IMPUTE2_tped(self,ped,d_chromosome_lengths,sepjoin,d_centromere_ranges,):

        stop_check_fd_write_at_end_is_generic_parse_some_prefix_from_cmd_line

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
            ) = parse_line_gen_old(line_gen)

        ##
        ## loop over markers in tped file
        ##
        l_r2_allelic = []
        pos_add = 0
        print('chromosome', chromosome)
        fd_tped = open('%s.tped' %(ped),'r')
        for line_tped in fd_tped:
            l_tped = line_tped.strip().split()
            chromosome_tped = l_tped[0]
            marker_tped = l_tped[1]
            pos_tped = int(l_tped[3])
            if chromosome != chromosome_tped:
                print('chromosome', chromosome_tped)
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
                    ) = parse_line_gen_old(line_gen)
            if pos_tped < pos_gen:
                continue
            elif pos_tped > pos_gen:
                for line_gen in fd_gen:
                    (
                        l_gen, marker_gen, pos_gen,
                        alleleA_gen, alleleB_gen
                        ) = parse_line_gen_old(line_gen)
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


    def IMPUTE2_low_mem(self,n_individuals):

        sum_xy = 0
        sum_xx = 0
        sum_yy = 0
        sum_x = 0
        sum_y = 0

        fd = fd_map_omni = open('omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.map','r')
        fd_ped_omni = open('omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.ped','r')
        print(len(fd_map_omni.readlines()))
        print(len(fd_ped_omni.readlines()))
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
                print(pos_wgs, pos_omni)
                if int(pos_wgs) > int(pos_omni):
                    bool_continue = False
                elif pos_wgs == pos_omni:
                    print(line)
                    print(line_map_omni)
                    stop0
                    bool_continue = False
                ## wgs still lagging
                else:
                    continue

            for line_map_omni in fd_map_omni:
                pos_omni = line_map_omni.split()[3]
                if line_map_omni[:2] != 'xx':
                    continue
                print(pos_wgs, pos_omni)
                if int(pos_wgs) < int(pos_omni):
                    bool_continue = True
                    break
                stop
                if pos_wgs == pos_omni:
                    print(line)
                    print(line_map_omni)
                    stop1

            if bool_continue == True:
                continue
            if d.has_key(pos_wgs):
                print(d[pos_wgs])
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


    def omni2dic(self,):

        print('read map')

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


    def BEAGLE_count(self,):

        ##
        d_map, l_map = omni2dic()

        ##
        ## loop chromosomes
        ##
        print('read gprobs')
        l_wgs_markers_all = []
        wgs_not_omni = 0
        wgs = 0
        wgs_omni_intersect = 0
        for chromosome in xrange(1,22+1):
            print(chromosome)
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
        print('wgs', wgs)
        print('wgs not omni', wgs_not_omni)
        print('wgs and omni overlap', wgs_omni_intersect)

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
        print('omni', len(set_omni_markers))
        print('omni not wgs', len(set_omni_markers-set_wgs_markers))
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
                print(d_map[marker])
                print(l[:6])
                for i_sample in xrange(n_samples):
    ##                x = dosage_imputation = 0*l[6+i_sample+0]+l[6+i_sample+1]+2*l[6+i_sample+1]
                    x = dosage_imputation = l[6+i_sample+1]+2*l[6+i_sample+1]
                    sample = l_samples[i_sample]
                    alleles = d_ped[sample][marker]
                    print(alleleA_imputation)
                    print(alleleB_imputation)
                    alleleA_omni = alleles[0]
                    alleleB_omni = alleles[1]
                    print(alleleA_omni)
                    print(alleleB_omni)
                    ## homozygous
                    if alleleA_omni == alleleA_imputation and alleleB_omni == alleleA_imputation:
                        y = dosage_genotype = 0
                    ## homozygous
                    elif alleleA_omni == alleleB_imputation and alleleB_omni == alleleB_imputation:
                        y = dosage_genotype = 2
                    ## homozygous
                    else:
                        print(alleleA_omni, alleleA_imputation)
                        print(alleleB_omni, alleleB_imputation)
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


    def bgl_gprobs(self,ped,d_chromosome_lengths,sepjoin,d_centromere_ranges,):

        stop_check_fd_write_at_end_is_generic_parse_some_prefix_from_cmd_line

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
        print('chromosome', chromosome)
        fd_tped = open('%s.tped' %(ped),'r')
        for line_tped in fd_tped:
            l_tped = line_tped.strip().split()
            chromosome_tped = l_tped[0]
            marker_tped = l_tped[1]
            pos_tped = int(l_tped[3])
            if chromosome != chromosome_tped:
                print('chromosome', chromosome_tped)
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


    def gnuplot(self,prefix,d_chromosome_lengths,d_centromere_ranges,):

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

        execmd('gnuplot %s.gnuplotsettings' %(prefix))

    ##    cmd = 'convert -rotate 90 %s.ps %s.png' %(prefix,prefix,)
        cmd = 'convert %s.ps %s.png' %(prefix,prefix,)
        print(cmd)
        execmd(cmd)
        if os.path.isfile('%s.png' %(prefix)):
            os.remove('%s.ps' %(prefix))
            os.remove('%s.gnuplotsettings' %(prefix))

        return


    def __init__(self,):

        parser = argparse.ArgumentParser()

        ## /lustre/scratch107/projects/uganda_gwas/QC/omni2.5-8_20120809_gwa_uganda_gtu_flipped.postQC.autosomes.bed
        ## /lustre/scratch107/projects/agv/users/tc9/WGS/uganda_20130113/out_ProduceBeagleInput/ProduceBeagleInput.bgl

        parser.add_argument(
            '--fp1',
            dest='l_fp1',nargs='+',
            help='List of file paths 1 in sequential order',
            metavar='FILE',default=None,
            required = True,
            )

        parser.add_argument(
            '--fp2',
            dest='l_fp2',nargs='+',
            help='List of file paths 2 in sequential order (e.g. chrom1 chrom2...).',
            metavar='FILE',default=None,
            required = True,
            )

        parser.add_argument(
            '--affix',
            dest='affix',
            help='File prefix/suffix',
            metavar='FILE',default=None,
            required = True,
            )

        parser.add_argument(
            '--fp_samples1',
            dest='fp_samples1',
            help='File with sample IDs in sequence corresponding to sequence of genotype probabilities',
            metavar='FILE',default=None,
            required = False,
            )

        parser.add_argument(
            '--fp_samples2',
            dest='fp_samples2',
            help='File with sample IDs in sequence corresponding to sequence of genotype probabilities',
            metavar='FILE',default=None,
            required = False,
            )

        parser.add_argument(
            '--fp_sampledic',
            dest='fp_sampledic',
            help='Translation of sample IDs from file1 to file2',
            metavar='FILE',default=None,
            required = False,
            )

        parser.add_argument(
            '--format1',
            dest='format1',
            metavar='STRING',default=None,
            required = False,
            )

        parser.add_argument(
            '--format2',
            dest='format2',
            metavar='STRING',default=None,
            required = False,
            )

        ## http://docs.python.org/2/library/functions.html#vars
        for k,v in vars(parser.parse_args()).items():
            setattr(self,k,v)
            continue

        for l_fp in [self.l_fp1,self.l_fp2,]:
            for fp in l_fp:
                if not os.path.isfile(fp):
                    print('missing', fp)
                    sys.exit(0)

        self.extension1 = self.l_fp1[0][self.l_fp1[0].rindex('.')+1:]
        self.extension2 = self.l_fp2[0][self.l_fp2[0].rindex('.')+1:]

        if self.fp_samples1 == None:
            if self.extension1 == 'tped':
                self.fp_samples1 = '%s.tfam' %(self.l_fp1[0][:-5])

        return


if __name__ == '__main__':
    self = main()
    self.main()
