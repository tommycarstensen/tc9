#!/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

## http://www.broadinstitute.org/gatk/about#typical-workflows
## http://www.broadinstitute.org/gatk/guide/topic?name=best-practices

## built-ins
import math, os, sys, time, re
## New in version 2.7
import argparse

##
## todo20120724: tc9: split chromosomes into smaller parts before BEAGLE step
## to avoid requesting 16gb of memory as those nodes are limited
##

##
## todo20120816: tc9: automate split bgl file step
## between producebeagleinput and BEAGLE
##

##
## todo20120809: tc9/dg11: add a ReduceReads step before UnifiedGenotyper for speed purposes
## and for better variant calling? cf. slide 14 of https://www.dropbox.com/sh/e31kvbg5v63s51t/ajQmlTL6YH/ReduceReads.pdf
## "VQSR Filters are highly empowered by calling all samples together"
## "Reduced BAMs provides better results for large scale analysis projects ( > 100 samples) because it doesn't require batching."
## http://www.broadinstitute.org/gatk/guide/topic?name=best-practices
## "Even for single samples ReduceReads cuts the memory requirements, IO burden, and CPU costs of downstream tools significantly (10x or more) and so we recommend you preprocess analysis-ready BAM files with ReducedReads."
##

## todo20130204: tc9: make memory sample size dependent... only tested on 3 datasets with 100 samples each...

## todo20130206: tc9: mention problem with splitting and centromeres/telomers to ms23 and dg11 and come up with a better solution if necessary... for now opt for the quick/easy solution...

## todo20130207: tc9: get rid of orphant functions

## todo20130210: tc9: download markers from http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes.phase1_release_v3/ instead of asking user for their location by default... Make those arguments non-required... --- alternatively download as specified here http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes.phase1_release_v3/READ_ME_beagle_phase1_v3

class main():

    def main(self):

        ## define list of chromosomes
        l_chroms = [str(i) for i in xrange(1,22+1,)]+['X','Y',]

        self.init(l_chroms,)

        self.check_logs()

        ## parse chromsome lengths from reference sequence
        d_chrom_lens = self.parse_chrom_lens()

        ##
        ## write shell scripts
        ##
        self.UnifiedGenotyper(l_chroms,d_chrom_lens,)

        self.VariantRecalibrator(l_chroms,d_chrom_lens,)

        self.ApplyRecalibration(l_chroms,d_chrom_lens,)

        self.ProduceBeagleInput()

        self.BEAGLE(l_chroms,)

        self.IMPUTE2(l_chroms,d_chrom_lens,)

        return


    def check_logs(self,):

        bool_error = False
        l_fn = os.listdir('stderr')
        for fn in l_fn:
            if os.path.isdir('stderr/%s' %(fn)): continue
            if os.path.getsize('stderr/%s' %(fn)) == 0: continue
            fd = open('stderr/%s' %(fn),'r')
            s = fd.read()
            fd.close()
            print s
            print 'stderr/%s' %(fn)
            bool_error = True
            break

        if bool_error == True:
            sys.exit(0)

        return


    def init(self,l_chroms,):

        self.mkdirs(l_chroms,)

        return


    def BEAGLE_unite(self,chrom):

        fd = open('BEAGLE_divide_indexes.txt','r')
        lines = fd.readlines()
        fd.close()
        d_index2pos = {}
        for line in lines:
            l = line.strip().split()
            if l[0] != chrom: continue
            index = int(l[1])
##            pos = int(l[2])
            pos1 = int(l[2])
            pos2 = int(l[3])
##            d_index2pos[index] = pos
            d_index2pos[index] = pos1
            d_index2pos[index] = pos2
        fp_out = 'out_BEAGLE/%s/' %(chrom,)
        fp_out += 'BeagleOutput.%s.bgl.gprobs' %(chrom,)
        for index in xrange(1,max(d_index2pos.keys())+1,):
            print 'BEAGLE_unite', chrom, index
            fp_in = 'out_BEAGLE/%s/%s.%i.gprobs' %(chrom,chrom,index,)
            if index == 1:
                cmd = 'head -n1 %s > %s' %(fp_in,fp_out,)
                os.system(cmd)
##            pos = d_index2pos[index]
            pos1 = d_index2pos[index][0]
            pos2 = d_index2pos[index][1]
            cmd = "awk 'NR>1{pos=substr($1,%i);" %(len(chrom)+2)
            cmd += "if(pos>=%i&&pos<%i) print $1}'" %(
##                pos,pos+1000000*self.i_BEAGLE_size)
                pos1,pos2,)
            cmd += ' %s >> %s' %(fp_in,fp_out,)
            os.system(cmd)

        return


    def gprobs2gen():

        ##
        ## IMPUTE2 takes gen files...
        ##
        ##
        ## gprobs > gen
        ##
        ## skip header
        s_more = 'more +2 %s' %(fp_gprobs,)
        ## print columns 1 and 2 when using field separator : (i.e. a replacement of field operator)
        s_awk1 = "awk -F : '{print $1, $2}'"
        ## print columns 1 and 2 with default field separator (i.e. append columns)
        s_awk2 = '''awk '{print $1":"$2, $1":"$2, $0}' '''
        ## cut out column 3 and do *not* print this column (i.e. remove "marker" from header and chromosomeID from subsequent rows)
        s_cut = 'cut -d " " -f 3 --complement > out_BEAGLE/BeagleOutput.$CHROMOSOME.gen'
        s_convert = ' | '.join([s_more,s_awk1,s_awk2,s_cut])
        cmd += '\n'+s_convert+'\n\n'
##        print s_convert
##        os.system(s_convert)
##        stop

        return


    def IMPUTE2_parse_input_files(self,):

        ## open file
        fd = open('BEAGLE_divide_indexes.txt','r')
        ## read all lines into memory
        lines = fd.readlines()
        ## close file
        fd.close()

        l_fp_in = []
        for line in lines:
            l = line.strip().split()
            chrom = l[0]
            index = int(l[1])
            fp_in = 'out_BEAGLE/%s/' %(chrom,)
            fp_in += 'BeagleOutput.%s.bgl.%i.' %(chrom,index,)
            fp_in += 'ProduceBeagleInput.%s.bgl.%i.gprobs' %(chrom,index,)
            l_fp_in += [fp_in]

        return l_fp_in


    def IMPUTE2(self,l_chroms,d_chrom_lens,):

        ##
        ## 1) check input existence
        ##
        l_fp_in = self.IMPUTE2_parse_input_files()
        bool_exit = self.check_in('BEAGLE',l_fp_in,)

        ##
        ## 2) touch
        ##
        bool_return = self.touch('IMPUTE2')
        if bool_return == True: return

        for chrom in l_chroms:
##            if chrom != '22': continue ## tmp!!!
            self.BEAGLE_unite(chrom)

        self.gprobs2gen()

        ##
        ## gprobs > sample
        ##
        fp_out = 'out_BEAGLE/chromosome$CHROMOSOME.sample'
        cmd += '''echo "ID_1 ID_2 missing\n0 0 0" > %s ; \
        head -1 %s \
        | sed 's/ /\\n/g' \
        | fgrep -v "marker" | fgrep -v "allele" \
        | uniq \
        | awk '{print $1, $1, "NA"}' \
        >> %s''' %(
            fp_out, fp_gprobs, fp_out,
            )
        cmd += '\n\n'

##        ## gen2ped (gen and sample to ped and map)
##        ## http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html
##        out_prefix = 'out_GTOOL/gtool.%s' %(chrom)
##        cmd += 'gtool=/nfs/team149/Software/usr/share/gtool/gtool'
##        cmd += '$gtool -G '
##        ## in
##        cmd += '--g out_BEAGLE/BeagleOutput.%s.gen ' %(chrom)
##        cmd += '--s out_BEAGLE/chromosome%s.sample ' %(chrom)
##        ## out
##        cmd += '--ped %s.ped ' %(out_prefix)
##        cmd += '--map %s.map ' %(out_prefix)
##        ## parameters
##        cmd += '--phenotype phenotype_1 --threshold 0.9 --snp'
##        cmd += '\n'


        fp_in = 'out_BEAGLE/BeagleOutput.$CHROMOSOME.gen'

##        intervals_IMPUTE2 = int(math.ceil(
##            d_chrom_lens[chrom]
##            /self.bps_per_interval_IMPUTE2
##            ))

        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html
        ## /lustre/scratch107/projects/agv/imputation/code_FINAL/imputechunks.sh
        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2_instructions.pdf
        ## http://mathgen.stats.ox.ac.uk/impute/example_one_phased_panel.html
        ## http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html

        lines = []

        lines += ['#!/bin/bash']

        lines += ['CHROMOSOME=$1']
        lines += ['LENCHROMOSOME=$2']

        ##
        ## define files
        ##
        ## Input file options
        ## http://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1integrated_v3_impute.tgz
        s_haps = self.fp_impute2_hap
        s_legend = self.fp_impute2_legend

        print self.fp_impute2_map
        stop
        lines += ['if [ $CHROMOSOME=="X" ]\nthen']
        lines += ['map=%s' %(self.fp_impute2_map.replace('$CHROMOSOME','$CHROMOSOME_PAR1',),)]
        lines += ['else']
        lines += ['map=%s' %(self.fp_impute2_map)]
        lines += ['fi']
        
        ## http://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-g
        fn_genotype    = '%s' %(fp_in)
##        ## Strand alignment options
##        s_strand       = '/lustre/scratch107/projects/agv/imputation/trial3/mock_strand_file/strand_file_chr%s.txt' %(chrom)

        ## IMPUTE2 requires that you specify an analysis interval
        ## in order to prevent accidental whole-chromosome analyses.
        ## If you want to impute a region larger than 7 Mb
        ## (which is not generally recommended),
        ## you must activate the -allow_large_regions flag.
        lines += ['posmax=$((($LSB_JOBINDEX+0)*5000000))']
        lines += ['if test $posmax -ge $LENCHROMOSOME']
        lines += ['then posmax=$LENCHROMOSOME']
        lines += ['fi\n']
        s = '\n'.join(lines)

        ##
        ## initiate impute2 command
        ##
        s += '%s \\' %(self.fp_software_IMPUTE2)

        ##
        ## Required arguments
        ## http://mathgen.stats.ox.ac.uk/impute/required_arguments.html
        ##
        ## http://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-g
        s += '-g %s \\' %(fn_genotype,)
        ## http://mathgen.stats.ox.ac.uk/impute/required_arguments.html#-m
        ## http://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-m
        ## Ask Deepti whether a different recombination map file should be used
        s += '-m $map \\'
        ## http://mathgen.stats.ox.ac.uk/impute/basic_options.html#-int
        s += '-int $((($LSB_JOBINDEX-1)*%i+1)) $posmax \\' %(5000000)

        ##
        ## Basic options
        ## http://mathgen.stats.ox.ac.uk/impute/basic_options.html
        ##

        ##
        ## Output file options
        ## http://mathgen.stats.ox.ac.uk/impute/output_file_options.html
        ##
        s += '-o out_IMPUTE2/impute2.gen.$CHROMOSOME.$LSB_JOBINDEX \\'

        ##
        ## Input file options
        ##
        ## http://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-l
        s += '-l %s \\' %(s_legend,)
        ## http://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-h
        s += '-h %s \\' %(s_haps,)

        self.write_shell('shell/IMPUTE2.sh',s,)

        ##
        ## execute shell script
        ##
        ## bsub -w "done(BEAGLE%s[*])" %(chrom)

        return s


    def write_shell(self,fp,lines,):

        if type(lines) != list:
            print type(lines)
            stop

        s = '\n'.join(lines)+'\n\n'
        fd = open(fp,'w')
        fd.write(s)
        fd.close()
        os.system('chmod +x %s' %(fp))

        return


    def BEAGLE_parse_input(self,chrom,):

        l_fps = []
        for fp_in in [self.fp_BEAGLE_phased,self.fp_BEAGLE_markers,]:
            fd = open(fp_in,'r')
            lines = fd.readlines()
            fd.close()
            for line in lines:
                l = line.strip().split()
                if l[0] == chrom:
                    l_fps += [l[1]]
                    break
                elif l[0] == '$CHROM':
                    l_fps += [l[1].replace('$CHROM',chrom)]
                    break
                continue

        fp_phased = l_fps[0]
        fp_markers = l_fps[1]

        return fp_phased, fp_markers


    def BEAGLE_divide(self,chrom,):

        ## To save disk space I should rewrite this function to act directly on ProduceBeagleInput.bgl instead of ProduceBeagleInput.$CHROM.bgl
        ## I can't remember why I split up the bgl file per chromosome in the first place...

        ##
        ## size and edge
        ##
        size = self.i_BEAGLE_size*1000000
        edge = self.i_BEAGLE_edge*1000

        fp_in = 'out_ProduceBeagleInput/ProduceBeagleInput.%s.bgl' %(chrom)
        ## genotype likelihood file
        fp_out_prefix = 'in_BEAGLE/%s' %(chrom)

        fp_phased, fp_markers = self.BEAGLE_parse_input(chrom,)

        ##
        ## parse genotype likelihood file header
        ##
        ## append header
        fd = open('out_ProduceBeagleInput/ProduceBeagleInput.bgl','r')
        header = fd.readline()
        fd.close()

##        ## parse header
##        header_phased = fd_phased.readline()

        ##
        ## open files
        ##
        fd_phased = open(fp_phased,'r')
        fd_markers = open(fp_markers,'r')
        fd_bgl = open(fp_in,'r') ## with header
        fd_bgl.readline() ## skip header

        ##
        ## parse phased header and skip phased header
        ##
        header_phased = fd_phased.readline()

        ##
        ## initiate lines out
        ##
        lines_out_phased1 = [header_phased]
        lines_out_phased2 = [header_phased]
        lines_out_markers1 = []
        lines_out_markers2 = []

        ##
        ## set smallest fragment size
        ##
        min_bps = 400000
##        min_bps = size ## tmp!!!
##        min_bps = 480000 ## tmp!!!

        ##
        ## get position of first SNP from genotype likelihood file
        ##
        lines_out1 = [header]
        lines_out2 = [header]
        for line in fd_bgl:
            lines_out1 += [line]
            l = line.strip().split()
            position = int(l[0].split(':')[1])
            pos_init1 = position-position%size
            pos_term1 = self.calc_pos_term(pos_init1,position,size,min_bps)
            break

        d_index2pos = {}
        index = 1 ## LSF does not allow 0 for LSB_JOBINDEX! "Bad job name. Job not submitted." ## e.g. bsub -J"test[0-2]" echo A
        pos_markers = 0
        bool_append_markphas = False
        bool_append_to_prev = False
        bool_append_to_next = False
        pos_prev = 0
        bool_BOF2 = False
        bool_EOF1 = False
        ## loop over BEAGLE genotype likelihood file lines
        for line in fd_bgl:
            l = line.strip().split()
            ## assume markers to be formatted like CHROM:POS
            position = int(l[0].split(':')[1])

            ##
            ## loop to determine
            ## whether markers and phased should be appended or not
            ##
            while True:
                if pos_markers == position:
                    bool_append_markphas = True
                    break
                elif pos_markers > position:
                    bool_append_markphas = False
                    break
                else:
                    line_markers = fd_markers.readline()
                    line_phased = fd_phased.readline()
                    l_markers = line_markers.split()
                    pos_markers = int(l_markers[1])
                    continue
            if pos_markers == position and bool_append_markphas == False:
                print pos_markers, position, bool_append_markphas
                stop1
            if pos_markers != position and bool_append_markphas == True:
                print pos_markers, position, bool_append_markphas
                stop2

            ##
            ## avoid multiple comparisons of large integers
            ## by means of booleans instead of nesting
            ##
            if bool_BOF2 == False and position > pos_term1-edge:
                bool_BOF2 = True
                pos_init2 = pos_term1
                ## centromere
                if position > pos_term1+edge:
                    bool_EOF1 = True
                    if pos_prev-pos_init1 < min_bps:
                        bool_append_to_prev = True
            elif bool_BOF2 == True and bool_EOF1 == False and position > pos_term1+edge:
                bool_EOF1 = True
            ## centromere
            elif bool_EOF1 == False and position > pos_term1+edge:
                print bool_BOF2 == True
                pos_init2 = position-position%size
                bool_BOF2 = True
                bool_EOF1 = True
                print 'position', position
                print 'pos_prev', pos_prev
                print 'pos_init1', pos_init1
                print 'pos_term1', pos_term1
                print 'pos_init2', pos_init2
                print 'size', size
                print edge
                print 'bool_BOF2', bool_BOF2
                print 'bool_EOF1', bool_EOF1
                print len(lines_out1)
                print len(lines_out2)
                print bool_append_to_previous
                stop
            else:
                if bool_BOF2 == True and bool_EOF1 == False:
                    pass
                elif bool_BOF2 == True or bool_EOF1 == True:
                    print bool_BOF2, bool_EOF1
                    print position
                    print pos_init1
                    print pos_term1
                    stop
##            ## else reset starting position and continue
####            elif pos_prev-pos_init1 < min_bps: ## slow
##            else: ## faster
##                bool_append_to_previous = True
####                pos_init1 = position-position%size
##                ## centromere
##                if position > pos_init1+size:
##                    bool_BOF2 = True
##                    bool_EOF1 = True
##                    pos_init2 = position-position%size
##                    print 'centromere',
##                    pos_init1, pos_prev, pos_init2, position
##                    stop
##                pass
####            if pos_prev-pos_init1 >= min_bps:

            ## first position after centromere
            if [chrom,position] in [
                ['1',144888371 ,],
                ['2',95329645,],
                ['3', 93506547 ,],
                ['4', 52682571 ,],
                ['5', 49441621 ,],
                ['6' ,61880512 ,],
                ['7', 61060256 ,],
                ['8', 46883398 ,],
                ['9', 66334018 ,],
                ['10', 42397373 ,],
                ['11', 54794727 ,],
                ['12', 37857751 ,],
                ['13', 19168767 ,],
                ['14', 19099878 ,],
                ['15', 20004821 ,],
                ['16', 46436844 ,],
                ['17', 25308930 ,],
                ['18', 18520354 ,],
                ['19', 27740659 ,],
                ['20', 29444239 ,],
                ['21', 9412629 ,],
                ['21', 14594887 ,],
                ['22', 16554965 ,],
                ['X', 61731996 ,],
                ['Y', 2657176 ,],
                ['Y', 13196611 ,],
                ['Y', 28505204 ,],
                ]:
                print '@@@@@', chrom,position,bool_BOF2,bool_EOF1
                if bool_BOF2 == False or bool_EOF1 == False:
                    stop
            ## last position before centromere
            if [chrom,position] in [
                ['1',121351046,],
                ['2',92304874,],
                ['3', 90502860,],
                ['4', 49632969,],
                ['5', 46405262,],
                ['6' ,58775200,],
                ['7', 58053551,],
                ['8', 43820454,],
                ['9', 44905287,],
                ['10', 39154363,],
                ['11', 51576318,],
                ['12',  34856284,],
                ['13',  19168767,],
                ['14',  19099878,],
                ['15',  20004821,],
                ['16',  35281785,],
                ['17',  22244499,],
                ['18',  15384732,],
                ['19',  24600804,],
                ['20',  26316584,],
                ['21',  9412629,],
                ['21',  10980490,],
                ['22',  16554965,],
                ['X',  58563466,],
                ['Y',  2657176,],
                ['Y',  10077460,],
                ['Y',  24874071,],
                ]:
                print '@@@@@', chrom,position,bool_BOF2,bool_EOF1
                if bool_BOF2 == True or bool_EOF1 == True:
                    stop

            ## BOF2
            if bool_BOF2 == True:
                lines_out2 += [line]
                if bool_append_markphas == True:
                    lines_out_markers2 += [line_markers]
                    lines_out_phased2 += [line_phased]
            if bool_EOF1 == False:
                if position>pos_term1+edge:
                    print pos_init1, position
                    print pos_term1
                    print bool_BOF2
                    print (bool_BOF2 == True and bool_EOF1 == False and position > pos_term1+edge)
                    stop
                lines_out1 += [line]
                if bool_append_markphas == True:
                    lines_out_markers1 += [line_markers]
                    lines_out_phased1 += [line_phased]
            ## EOF1
            else:
                ## append to previous file if less than 1000 variants
                if bool_append_to_prev == False and len(lines_out1) < 1000:
                    ## 2nd append to prev
                    if pos_init1-size == d_index2pos[index-1][0]:
                        bool_append_to_prev = True
                    ## append to next
                    else:
                        bool_append_to_next = True
                if bool_append_to_prev == True:
                    print 'append prev'
                    mode = 'a'
                    ## remove headers
                    lines_out1 = lines_out1[1:]
                    lines_out_phased1 = lines_out_phased1[1:]
                    index -= 1
                    pos_init1 = d_index2pos[index][0]
                elif bool_append_to_next == True:
                    print 'append next'
##                    pos_term1 += size
##                    pos_prev = position
##                    bool_append_to_next = False
##                    lines_out_phased2 = [header_phased]
##                    lines_out_markers2 = []
##                    ## do not write/append lines yet
##                    ## merge with next fragment instead
##                    continue
                else:
                    mode = 'w'
                print 'pos_curr %9i, pos_init %9i, pos_term %9i, i %3i, n %5i' %(
                    position,pos_init1,pos_term1,index,len(lines_out1))
                if bool_append_to_next == False:
                    ## write/append lines
                    fd_out = open('%s.%i.like' %(fp_out_prefix,index,),mode)
                    fd_out.writelines(lines_out1)
                    fd_out.close()
                    fd = open('%s.%i.markers' %(fp_out_prefix,index,),mode)
                    fd.writelines(lines_out_markers1)
                    fd.close()
                    fd = open('%s.%i.phased' %(fp_out_prefix,index,),mode)
                    fd.writelines(lines_out_phased1)
                    fd.close()
                    ## replace current lines with next lines
                    lines_out1 = lines_out2
                    lines_out_phased1 = lines_out_phased2
                    lines_out_markers1 = lines_out_markers2
                ## reset next lines
                lines_out2 = [header]
                lines_out_phased2 = [header_phased]
                lines_out_markers2 = []
                if bool_append_to_next == False:
                    ## append index and positions to dictionary
                    d_index2pos[index] = [pos_init1,pos_term1,]
                    index += 1
                ## set pos init/term
                if bool_append_to_next == False:
                    pos_init1 = position-position%size
                    pos_term1 = self.calc_pos_term(pos_init1,position,size,min_bps)
                else:
                    pos_term1 += size
                ## reset booleans
                bool_BOF2 = False
                bool_EOF1 = False
                bool_append_to_prev = False
                bool_append_to_next = False

            pos_prev = position

            ## continue loop over genotype likelihoods
            continue

        ## write remaining few lines to output files

        ## less than minimum number of base pairs
        ## or minimum number of variants
        if position-pos_init1 < min_bps or len(lines_out1) < 1000:
            mode = 'a'
            index -= 1
            pos_init1 = d_index2pos[index][0]
            print index, d_index2pos[index]
            print pos_init1,pos_term1
            stoptmp
            ## remove headers
            lines_out1 = lines_out1[1:]
            lines_out_phased1 = lines_out_phased1[1:]
        else:
            mode = 'w'
            if len(lines_out1) < 1000:
                print len(lines_out1)
                stop
        print 'pos_curr %9i, pos_init %9i, pos_term %9i, i %3i, n %5i' %(
            position,pos_init1,pos_term1,index,len(lines_out1))

        ##
        ## write/append lines
        ##
        fd_out = open('%s.%i.like' %(fp_out_prefix,index,),mode)
        fd_out.writelines(lines_out1)
        fd_out.close()

        fd = open('%s.%i.markers' %(fp_out_prefix,index,),mode)
        fd.writelines(lines_out_markers1)
        fd.close()

        fd = open('%s.%i.phased' %(fp_out_prefix,index,),mode)
        fd.writelines(lines_out_phased1)
        fd.close()

        ## append position to dictionary
        d_index2pos[index] = [pos_init1,pos_term1,]

        ##
        ## close all files after looping over last line
        ##
        fd_bgl.close()
        fd_markers.close()
        fd_phased.close()

        return d_index2pos


    def calc_pos_term(self,pos_init1,position,size,min_bps):

        ## acrocentric chromosome telomeres (13,14,15,21,22,Y)
        ## and all centromeres
        if size-min_bps < position-pos_init1:
            if pos_init1 != 0:
                print pos_init1, pos_term1
            pos_term1 = pos_init1+2*size
            if pos_init1 != 0:
                print pos_init1, pos_term1
                print pos_init1, position, pos_term1
                stoptmp
        ## non-acrocentric chromosome
        ## and all cemtromeres
        else:
            pos_term1 = pos_init1+1*size

        return pos_term1


    def BEAGLE(self,l_chroms,):

        ## http://faculty.washington.edu/browning/beagle/beagle_3.3.2_31Oct11.pdf

        ## http://www.broadinstitute.org/gatk/guide/article?id=43
        ## Interface with BEAGLE imputation software - GSA
        ## IMPORTANT: Due to BEAGLE memory restrictions,
        ## it's strongly recommended that BEAGLE be run on a separate chromosome-by-chromosome basis.
        ## In the current use case, BEAGLE uses RAM in a manner approximately proportional to the number of input markers.

        ## dg11: "Then run beagle imputation
        ## (you may have to chunk up for imputation
        ## and use the "known" option with the 1000G reference)"
        ## ask Deepti what the "known" option is...

        ## fgrep Mem stdout/BEAGLE/BEAGLE.22.* | sort -k5n,5 | tail -n1 | awk '{print $5}'
        ## 2885 with lowmem=true
        memMB = 4000

        ##
        ## 1) check input existence
        ##
        fp_in = 'out_ProduceBeagleInput/ProduceBeagleInput.bgl'
        bool_exit = self.check_in('ProduceBeagleInput',[fp_in,],)

        ##
        ## 2) touch
        ##
        bool_return = self.touch('BEAGLE')
        if bool_return == True: return

        ##
        ## write shell script
        ##
        self.BEAGLE_write_shell_script(memMB)

        ##
        ## chunk up for imputation due to memory requirements...
        ##
##        ## split bgl by chromosome
##        ## this should be part of the BEAGLE_divide function
##        ## to avoid looping over the same lines twice (a bit stupid at the moment)
##        cmd = 'cat out_ProduceBeagleInput/ProduceBeagleInput.bgl '
##        cmd += ''' | awk 'BEGIN{FS=":"} NR>1'''
##        cmd += '''{print>"out_ProduceBeagleInput/ProduceBeagleInput."$1".bgl"}' '''
##        print cmd
##        os.system(cmd)
        ## split further into smaller fragments
        d_indexes = {}
        for chrom in l_chroms:
            if chrom in ['1',]: continue ## tmp!!!
##            if chrom != '1': continue ## tmp!!!
##            if chrom != '22': continue ## tmp!!!
            d_chrom_lens = self.parse_chrom_lens()
            d_index2pos = self.BEAGLE_divide(chrom,)
            d_indexes[chrom] = d_index2pos
            continue
        stoptmpend
##            os.remove(
##                'out_ProduceBeagleInput/ProduceBeagleInput.%s.bgl' %(chrom))

        s = ''
        for chrom,d_index2pos in d_indexes.items():
##            for index,pos in d_index2pos.items():
            for index,[pos1,pos2,] in d_index2pos.items():
##                s += '%s %i %i\n' %(chrom,index,pos,)
                s += '%s %i %i %i\n' %(chrom,index,pos1,pos2,)
##        fd = open('BEAGLE_divide_indexes.txt','w')
        fd = open('BEAGLE_divide_indexes.txt','a') ## tmp!!!
        fd.write(s)
        fd.close()

        ##
        ## execute shell script
        ##
        for chrom in l_chroms:
##            if chrom != '1': continue ## tmp!!!
            J = '%s%s[%i-%i]' %('BEAGLE',chrom,1,max(d_indexes[chrom].keys()),)
            std_suffix = '%s/%s.%s.%%I' %('BEAGLE','BEAGLE',chrom)
            cmd = self.bsub_cmd(
                'BEAGLE',J,memMB=memMB,std_suffix=std_suffix,chrom=chrom,)
            print cmd
            os.system(cmd)

        return


    def BEAGLE_write_shell_script(self,memMB,):

        fp_out = 'out_BEAGLE/$CHROMOSOME/$CHROMOSOME.$LSB_JOBINDEX.bgl'

        ## initiate shell script
        lines = ['#!/bin/bash\n']

        ## parse chromosome from command line
        lines += ['CHROMOSOME=$1\n']

        ## init cmd
        lines += ['cmd="']

        ##
        ## initiate BEAGLE
        ##
        ## Could not find the main class: Djava.io.tmpdir=out_BEAGLE.  Program will exit.
##        lines += ['java -Xmx%im java.io.tmpdir=out_BEAGLE -jar %s \\' %(
        lines += ['java -Xmx%im -jar %s \\' %(
            memMB,self.fp_software_beagle)]

        lines += self.body_BEAGLE(fp_out,)

        ##
        ## terminate BEAGLE
        ##
        lines += [';']

        ##
        ## gunzip the output files after runs to completion
        ##
        l = []
        for suffix in ['dose','phased','gprobs',]:
            fp = '%s' %(fp_out)
            fp += '.ProduceBeagleInput.${CHROMOSOME}.bgl.${LSB_JOBINDEX}.%s.gz' %(
                suffix,
                )
            l += ['gunzip %s' %(fp)]
        s = ';'.join(l).replace('$CHROMOSOME','${CHROMOSOME}')
        lines += [s]

        ## term cmd
        lines += self.term_cmd('BEAGLE',[fp_out],)

        ## write shell script
        self.write_shell('shell/BEAGLE.sh',lines,)

        return


    def body_BEAGLE(self,fp_out,):

        fp_phased, fp_markers = self.BEAGLE_parse_input('$CHROMOSOME',)

        fp_like = 'in_BEAGLE/$CHROMOSOME.$LSB_JOBINDEX.like'
        fp_phased = 'in_BEAGLE/$CHROMOSOME.$LSB_JOBINDEX.phased'
        fp_markers = 'in_BEAGLE/$CHROMOSOME.$LSB_JOBINDEX.markers'

        lines = []

##like=<unphased likelihood data file> where <unphased likelihood data file> is
##the name of a genotype likelihoods file for unphased, unrelated data
##(see Section 2.2). You may use multiple like arguments if data from different
##cohorts are in different files.
        lines += [' like=%s \\' %(fp_like)]
####arguments for phasing and imputing data ...
#### Arguments for specifying files
## phased=<phased unrelated file> where <phased unrelated file> is the name of a
## Beagle file containing phased unrelated data (see Section 2.1).
## You may use multiple phased arguments if data from different cohorts are in
## different files.
##        lines += [' phased=in_BEAGLE/ALL.chr$CHROMOSOME.phase1_release_v3.20101123.filt.renamed.bgl \\']
        lines += [' phased=%s \\' %(fp_phased)]
####  unphased=<unphased data file>                     (optional)
####  phased=<phased data file>                         (optional)
####  pairs=<unphased pair data file>                   (optional)
####  trios=<unphased trio data file>                   (optional)
####  like=<unphased likelihood data file>              (optional)
##markers=<markers file> where <markers file> is the name of the markers file containing
## marker identifiers, positions, and alleles described in Section 2.4.
## The markers argument is optional if you specify only one Beagle file,
## and is required if you specify more than one Beagle file.
##        s += ' markers=/lustre/scratch107/projects/uganda/users/tc9/in_BEAGLE/ALL.chr%s.phase1_release_v3.20101123.filt.markers ' %(chrom)
        lines += [' markers=%s \\' %(fp_markers)]
####missing=<missing code> where <missing code> is the character or sequence of characters used to represent a missing allele (e.g. missing=-1 or missing=?).
#### The missing argument is required.
##        s += ' missing=? '

##nsamples=<number of samples> where <number of samples> is positive integer
##giving the number of haplotype pairs to sample for each individual during each
##iteration of the phasing algorithm. The nsamples argument is optional. The
##default value is nsamples=4. If you are phasing an extremely large sample
## (say > 4000 individuals), you may want to use a smaller nsamples parameter
## (e.g. 1 or 2) to reduce computation time.  If you are phasing a small sample
## (say < 200 individuals), you may want to use a larger nsamples parameter
## (say 10 or 20) to increase accuracy.
        lines += [' nsamples=%i \\' %(int(self.i_BEAGLE_nsamples))]

        lines += [' niterations=10 \\'] ## default 10

        lines += [' omitprefix=true \\'] ## default false

        lines += [' verbose=false \\'] ## default false

        lines += [' lowmem=true \\'] ## default false

        ## non-optional output prefix
        lines += [' out=%s \\' %(fp_out)]

        return lines


    def write_nonarray_scripts(self,d_shell,):

        ##
        ## write shell scripts for each non-array step
        ##
##        d_bsub_lines = {}
        for i_step in range(len(self.l_steps)):

            step = self.l_steps[i_step]

            if self.bool_sequential == True:

                ## no additional steps
                if i_step+1 == len(self.l_steps):
                    bsub_line_next = None
                ## this is not the final step; i.e. next step exists
                else:
                    step_next = self.l_steps[i_step+1]
                    if step_next == 'BEAGLE':
                        if step_next in self.d_memory.keys():
                            s_queue = self.d_queues[step_next]
                            memory_MB = self.d_memory[step_next]
                        else:
                            s_queue = 'normal'
                            memory_MB = 4000
                        bsub_line_next = self.generate_bsub_line(
                            step_next,
                            job_array_intervals = 24,
                            job_array_title = '%s' %(step_next[:6],),
                            bool_wait = False,
                            s_queue = s_queue,
                            memory_MB = memory_MB,
                            )
            else:
                bsub_line_next = None

            s = d_shell[step]
            s += '\n\nfi\n\n'
            self.write_shell_script(
                step,
                s,
                line_additional = bsub_line_next,
                )
##            print d_shell[step]
##            stop

        return


    def write_array_scripts_chromosome(self,d_shell,):

        ## useful if per-chromosome basis... not whole-genome and not sub-chromosome fragments...
        if step == 'ApplyRecalibration':
            s = '#!/bin/bash\n'
        else:
            s = '''#!/bin/bash
LSB_JOBINDEX=$LSB_JOBINDEX
if [ ${LSB_JOBINDEX} -eq 23 ]
then
CHROMOSOME="X"
elif [ ${LSB_JOBINDEX} -eq 24 ]
then
CHROMOSOME="Y"
else
CHROMOSOME=$LSB_JOBINDEX
fi
echo $CHROMOSOME
'''

        return


    def write_array_scripts_subchromosome(self,d_array,):
        
        ##
        ## write shell script for UnifiedGenotyper (job array and first step)
        ##
        for step in self.l_steps_array_subchromosome:

            s = '#!/bin/sh\n\n'
            for chrom in ['X','Y',]+range(1,22+1,):
                chrom = str(chromosome)
                s += '%s\n' %(d_array[step][chrom])
                s += 'fi\n'

        return


    def parse_chrom_lens(self):

        d_chrom_lens = {}

        ## 1000G chromosome ranges
        fn = '%s.fai' %(self.fp_FASTA_reference_sequence)
        fd = open(fn)
        lines = fd.readlines()
        fd.close()
        for line in lines:
            l = line.strip().split()
            chrom = l[0]
            chrom_len = int(l[1])
##            if not chromosome in ['X','Y',]:
##                chrom = int(chromosome)
            d_chrom_lens[chrom] = chrom_len
            ## break, when we reach the Y chromosome
            if chrom == 'Y':
                break
        
        return d_chrom_lens


    def mkdirs(self,l_chroms,):

        l_dn = [
            'touch',
            'stdout','stderr',
            'shell',
            'out_UnifiedGenotyper',
            'out_VariantRecalibrator',
            'out_ApplyRecalibration',
            'out_ProduceBeagleInput',
            'in_BEAGLE','in_IMPUTE2',
            'out_BEAGLE','out_IMPUTE2',
            ]

        ## create subdirs
        for dn in l_dn:
            if not os.path.isdir(dn):
                os.mkdir(dn)
            if dn != 'touch' and dn[:3] == 'out':
                if not os.path.isdir(os.path.join('touch',dn)):
                    os.mkdir(os.path.join('touch',dn))
        for std in ['out','err',]:
            for dn in l_dn:
                if dn[:3] != 'out': continue
                if not os.path.isdir('std%s/%s' %(std,dn[4:],)):
                    os.mkdir('std%s/%s' %(std,dn[4:],))

        for chrom in l_chroms:
            if not os.path.isdir('out_BEAGLE/%s' %(chrom)):
                os.mkdir('out_BEAGLE/%s' %(chrom))
            if not os.path.isdir('touch/out_BEAGLE/%s' %(chrom)):
                os.mkdir('touch/out_BEAGLE/%s' %(chrom))

        return


    def check_in(self,analysis_type,l_fp_in,):
        
##        fd = open('%s.touch' %(analysis_type),'r')
##        s = fd.read()
##        fd.close()
##        l_fp_out = s.split('\n')

        l = os.listdir(os.path.join('touch','out_%s' %(analysis_type)))
        l_fp_out = [os.path.join('out_%s' %(analysis_type),fn) for fn in l]

        bool_exit = False
        if len(set(l_fp_in)-set(l_fp_out)) > 0:
            print '%s and possibly other files not generated.' %(
                list(set(l_fp_in)-set(l_fp_out))[0])
            print '%s has not run to completion. Exiting.' %(analysis_type)
            bool_exit = True
            sys.exit()

        return bool_exit


    def ProduceBeagleInput(self,):

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_beagle_ProduceBeagleInput.html

        ## "After variants are called and possibly filtered,
        ## the GATK walker ProduceBeagleInputWalker will take the resulting VCF as input,
        ## and will produce a likelihood file in BEAGLE format."

        T = analysis_type = 'ProduceBeagleInput'
        ##fula_20120704/stdout/ProduceBeagleInput/ProduceBeagleInput.out:    Max Memory :      2154 MB
        ##uganda_20130113/stdout/ProduceBeagleInput/ProduceBeagleInput.out:    Max Memory :      2168 MB
        ##zulu_20121208/stdout/ProduceBeagleInput/ProduceBeagleInput.out:    Max Memory :      2171 MB
        memMB = 4000
        ##fula_20120704/stdout/ProduceBeagleInput/ProduceBeagleInput.out:    CPU time   :   2483.24 sec.
        ##zulu_20121208/stdout/ProduceBeagleInput/ProduceBeagleInput.out:    CPU time   :   4521.89 sec.
        ##uganda_20130113/stdout/ProduceBeagleInput/ProduceBeagleInput.out:    CPU time   :   4948.94 sec.
        queue = 'normal'

        bool_return = self.touch(analysis_type)
        if bool_return == True: return

        fp_out = 'out_ProduceBeagleInput/ProduceBeagleInput.bgl'

        fp_in = 'out_ApplyRecalibration/ApplyRecalibration.recalibrated.filtered.vcf'
        bool_exit = self.check_in('ApplyRecalibration',[fp_in],)

        lines = ['#!/bin/bash']

        lines += self.init_cmd(analysis_type,memMB,)

        ## GATKwalker, required, out
        lines += ['--out %s \\' %(fp_out,)]
        ## GATKwalker, required, in
        lines += ['--variant %s \\' %(fp_in)]

##        l_fp_out = ['out_%s/%s.%i.bgl' %(T,T,chrom) for chrom in l_chrom]
        lines += self.term_cmd(analysis_type,[fp_out],)

        self.write_shell('shell/ProduceBeagleInput.sh',lines,)

        ## execute shell script
        J = 'PBI'
        cmd = self.bsub_cmd(analysis_type,J,memMB=memMB,)
        os.system(cmd)

        return


    def ApplyRecalibration(self,l_chroms,d_chrom_lens,):

        '''
Validated human SNP data suggests that the Ti/TV should be ~2.1 genome-wide and ~2.8 in exons (ref ???)
http://www.broadinstitute.org/gsa/wiki/index.php/QC_Methods#SNP_callset_metrics
http://www.broadinstitute.org/gsa/wiki/index.php/Variant_quality_score_recalibration#Ti.2FTv-free_recalibration
http://www.broadinstitute.org/gsa/wiki/images/b/b2/TiTv_free_VQSR.pdf
"For new approach: just cut at 99% sensitivity"
http://www.broadinstitute.org/gsa/wiki/images/a/ac/Ngs_tutorial_depristo_1210.pdf
http://www.broadinstitute.org/gsa/wiki/images/b/bc/Variant_Recalibrator_Sanger_June_2010.pdf
##
http://www.broadinstitute.org/gsa/wiki/images/e/eb/FP_TITV.jpg
'''

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_ApplyRecalibration.html

        T = analysis_type = 'ApplyRecalibration'
        ##fula_20120704/stdout/ApplyRecalibration/ApplyRecalibration.out:    Max Memory :      2187 MB
        ##zulu_20121208/stdout/ApplyRecalibration/ApplyRecalibration.out:    Max Memory :      2210 MB
        ##uganda_20130113/stdout/ApplyRecalibration/ApplyRecalibration.out:    Max Memory :      2233 MB
        memMB = 4000
        ##fula_20120704/stdout/ApplyRecalibration/ApplyRecalibration.out:    CPU time   :   2316.11 sec.
        ##zulu_20121208/stdout/ApplyRecalibration/ApplyRecalibration.out:    CPU time   :   3448.45 sec.
        ##uganda_20130113/stdout/ApplyRecalibration/ApplyRecalibration.out:    CPU time   :   3575.38 sec.
        queue = 'normal'

        ##
        ## 1) check input existence
        ##
        fp_in_recal = 'out_VariantRecalibrator/VariantRecalibrator.recal'
        fp_in_tranches = 'out_VariantRecalibrator/VariantRecalibrator.tranches'
##        T_prev = self.seqsteps[self.seqsteps.index(analysis_type)-1]
##        l_fp_in = self.d_out[T_prev]
        bool_exit = self.check_in('VariantRecalibrator',[fp_in_recal,fp_in_tranches,],)

        ##
        ## 2) touch
        ##
        bool_return = self.touch(analysis_type)
        if bool_return == True: return
           
        ##
        ##
        ##
        l_vcfs_in = self.get_fps_in(
            l_chroms,d_chrom_lens,self.bps_per_interval,)

        fp_out_prefix = 'out_ApplyRecalibration/ApplyRecalibration.recalibrated.filtered'
        fp_out_vcf = '%s.vcf' %(fp_out_prefix)
        fp_out_idx = '%s.vcf.idx' %(fp_out_prefix)
        if os.path.isfile(fp_out_vcf):
            return

        ##
        ## initialize shell script
        ##
        lines = ['#!/bin/bash']

        cmd = self.determine_TS_level(fp_in_tranches)
        lines += [cmd]


        lines += self.init_cmd(analysis_type,memMB,)

        ## GATKwalker, required, in
        lines += ['--input %s \\' %(vcf) for vcf in l_vcfs_in]
        ## GATKwalker, required, out
        lines += ['--out %s \\' %(fp_out_vcf,)]
        ## GATKwalker, required, in
        lines += ['--recal_file %s \\' %(fp_in_recal)]
        lines += ['--tranches_file %s \\' %(fp_in_tranches)]
        ## GATKwalker, optional, in
        ## ts_filter_level should correspond to targetTruthSensitivity prior to novelTiTv dropping (see tranches plot)
        ## default is 99.00
##        lines += ['--ts_filter_level 99.0 \\']
        lines += ['--ts_filter_level $ts_filter_level \\']

        lines += self.term_cmd(analysis_type,[fp_out_vcf,fp_out_idx,],)

        self.write_shell('shell/ApplyRecalibration.sh',lines,)

        ##
        ## execute shell script
        ##
        J = 'AR'
        cmd = self.bsub_cmd(analysis_type,J,memMB=memMB,)
        os.system(cmd)

        return


    def determine_TS_level(self,fp_tranches,):

        cmd = 'ts_filter_level=$('
        ## loop over lines
        cmd += 'cat %s' %(fp_tranches)
        ## init awk
        cmd += " | awk '"
        ## BEGIN
        cmd += ' BEGIN{FS=","}'
        ## skip header and loop over rows
        cmd += ' NR>3{'
        ## if below novelTiTv ratio
        cmd += ' if($5<2.1) {'
        ## if first line
        cmd += ' if(NR==4) {print $1; exit}'
        ## else if not first line
        cmd += ' else {'
        cmd += ' novelTiTv2=$5; targetTruthSensitivity2=$1;'
        cmd += ' dy=(novelTiTv2-novelTiTv1);'
        cmd += ' dx=(targetTruthSensitivity2-targetTruthSensitivity1);'
        cmd += ' slope=dy/dx;'
        cmd += ' intercept=novelTiTv1-slope*targetTruthSensitivity1;'
        cmd += ' print (2.1-intercept)/slope; exit}'
        ## else if above novelTiTv ratio
        cmd += ' } else {targetTruthSensitivity1=$1;novelTiTv1=$5}'
        ## continue loop over rows
        cmd += ' }'
        ## term awk
        cmd += " '"
        cmd += ')'

        return cmd


    def get_fps_in(self,l_chroms,d_chrom_lens,bps_per_interval,):

        l_vcfs_in = []
        for chrom in l_chroms:
            intervals = int(math.ceil(
                d_chrom_lens[chrom]/bps_per_interval))
            for interval in range(1,intervals+1,):
                fp_in = 'out_%s/%s.%s.%i.vcf' %(
                    'UnifiedGenotyper','UnifiedGenotyper',
                    chrom,interval,
                    )
                l_vcfs_in += [fp_in]

        return l_vcfs_in


    def VariantRecalibrator(self,l_chroms,d_chrom_lens):

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html

        T = analysis_type = 'VariantRecalibrator'
        ##fula_20120704/stdout/VariantRecalibrator/VariantRecalibrator.out:    Max Memory :     12840 MB
        ##zulu_20121208/stdout/VariantRecalibrator/VariantRecalibrator.out:    Max Memory :     13710 MB
        ##uganda_20130113/stdout/VariantRecalibrator/VariantRecalibrator.out:    Max Memory :     13950 MB
        memMB = 16000
        ##fula_20120704/stdout/VariantRecalibrator/VariantRecalibrator.out:    CPU time   :   9635.38 sec.
        ##zulu_20121208/stdout/VariantRecalibrator/VariantRecalibrator.out:    CPU time   :  10015.61 sec.
        ##uganda_20130113/stdout/VariantRecalibrator/VariantRecalibrator.out:    CPU time   :  11657.10 sec.
        queue = 'normal'

        ##
        ## 1) check input existence
        ##
        l_vcfs_in = self.get_fps_in(
            l_chroms,d_chrom_lens,self.bps_per_interval,)
        bool_exit = self.check_in('UnifiedGenotyper',l_vcfs_in,)

        ##
        ## 2) touch
        ##
        bool_return = self.touch(analysis_type)
        if bool_return == True: return

        fp_tranches = 'out_VariantRecalibrator/VariantRecalibrator.tranches'
        fp_recal = 'out_VariantRecalibrator/VariantRecalibrator.recal'

        ##
        ## initiate GATK walker
        ##
        lines = self.init_cmd(analysis_type,memMB,)

        ## GATKwalker, required, in
        lines += ['--input %s \\' %(vcf) for vcf in l_vcfs_in]
        ## GATKwalker, required, out
        lines += ['--recal_file %s \\' %(fp_recal)]
        lines += ['--tranches_file %s \\' %(fp_tranches)]
        ## GATKwalker, required, in (ask Deepti about this...)
        lines += [
            '--use_annotation QD \\',
            '--use_annotation HaplotypeScore \\',
            '--use_annotation MQRankSum \\',
            '--use_annotation ReadPosRankSum \\',
            '--use_annotation MQ \\',
            '--use_annotation FS \\',
            '--use_annotation DP \\',
            ]
        ## GATKwalker, optional, in
        ## Ugandan QCed to be added..!
        fd = open('%s' %(self.fp_resources),'r')
        lines_resources = fd.readlines()
        fd.close()
        lines += ['%s \\' %(line.strip()) for line in lines_resources]

        l_TStranches = [100,]
        l_TStranches += [99+i/10. for i in range(9,0,-1,)]
        l_TStranches += [90+i/2. for i in range(18,-1,-1,)]
        s_TStranches = ''
        for TStranche in l_TStranches:
            s_TStranches += '--TStranche %.1f ' %(TStranche)
        lines += ['%s \\' %(s_TStranches)]

        lines += self.term_cmd(analysis_type,[fp_tranches,fp_recal,],)

        self.write_shell('shell/VariantRecalibrator.sh',lines,)

        J = 'VR'
        cmd = self.bsub_cmd(analysis_type,J,memMB=memMB,queue=queue,)
        os.system(cmd)

        return


    def init_cmd(self,analysis_type,memMB,):

        s = 'cmd="'
        ## run GATK
        s += 'java '
        ## set maximum heap size
        s += '-Xmx%im ' %(memMB)
        s += '-jar %s \\' %(self.fp_GATK)
        lines = [s]

        ## CommandLineGATK, required, in
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--analysis_type
        lines += ['--analysis_type %s \\' %(analysis_type)]
        ## CommandLineGATK, optional, in
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--reference_sequence
        lines += ['--reference_sequence %s \\' %(self.fp_FASTA_reference_sequence)]

        return lines


    def touch(self,analysis_type):

        bool_return = False
        fn_touch = '%s.touch' %(analysis_type)
        if os.path.isfile(fn_touch):
##            if time.time()-os.path.getmtime(fn_touch) < 60*60:
                if self.verbose == True:
                    print 'in progress or completed:', analysis_type
                bool_return = True
        os.system('touch %s' %(fn_touch))

        return bool_return


    def UnifiedGenotyper(self,l_chroms,d_chrom_lens,):

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html

        T = analysis_type = 'UnifiedGenotyper'
        ## fula_20120704/stdout/UnifiedGenotyper.16.9.out:    CPU time   :   9420.62 sec.
        ## zulu_20121208/stdout/UnifiedGenotyper/UnifiedGenotyper.10.12.out:    CPU time   :  12106.00 sec.
        queue = 'normal'
        ## zulu_20121208/stdout/UnifiedGenotyper/UnifiedGenotyper.16.5.out:    Max Memory :      1945 MB
        memMB = 2000

        ##
        ## 1) touch
        ##
        bool_return = self.touch(analysis_type)
        if bool_return == True: return

        fp_out = 'out_%s/%s.$CHROMOSOME.$LSB_JOBINDEX.vcf' %(
            analysis_type,analysis_type,)

        ## initiate shell script
        lines = ['#!/bin/bash']

        ## commands prior to GATK command
        lines += self.init_UnifiedGenotyper(
            l_chroms,d_chrom_lens,)

        ## initiate GATK command
        lines += self.init_cmd('UnifiedGenotyper',memMB,)

        ## append GATK command options
        lines += self.body_UnifiedGenotyper(fp_out,)

        ## terminate shell script
        lines += self.term_cmd(analysis_type,[fp_out],)

        ## write shell script
        self.write_shell('shell/UnifiedGenotyper.sh',lines,)

        ## execute shell script
        for chrom in l_chroms:
            fp_stdout = 'stdout/%s.UnifiedGenotyper.%s.1.out' %(
                self.project,chrom,)
            if os.path.isfile(fp_stdout):
                continue
            fp_out_chrom = fp_out.replace('$LSB_JOBINDEX','1')
            fp_out_chrom = fp_out_chrom.replace('$CHROMOSOME',chrom)
            if os.path.isfile(fp_out_chrom):
                continue
            intervals = int(math.ceil(
                d_chrom_lens[chrom]/self.bps_per_interval))

            J = '%s%s[%i-%i]' %('UG',chrom,1,intervals,)
            std_suffix = '%s/%s.%s.%%I' %(T,T,chrom)
            cmd = self.bsub_cmd(
                analysis_type,J,std_suffix=std_suffix,chrom=chrom,
                memMB=memMB,queue=queue,)
            os.system(cmd)

        return


    def bsub_cmd(
        self,
        analysis_type,
        J,
        queue='normal',memMB=4000,
        std_suffix=None,
        chrom=None,
        ):

        if not std_suffix:
            std_suffix = '%s/%s' %(analysis_type,analysis_type,)

        cmd = 'bsub -J"%s%s" -q %s' %(self.name,J,queue,)
        cmd += ' -G %s' %(self.project)
        cmd += " -M%i000 -R'select[mem>%i] rusage[mem=%i]'" %(
            memMB,memMB,memMB,)
        cmd += ' -o stdout/%s.out' %(std_suffix)
        cmd += ' -e stderr/%s.err' %(std_suffix)
        cmd += ' ./shell/%s.sh' %(analysis_type,)
        if chrom != None:
            cmd += ' %s' %(chrom)

        return cmd


    def body_UnifiedGenotyper(self,fp_out,):

        lines = []

        ##
        ## required
        ##
        ## File to which variants should be written
        lines += [' --out %s \\' %(fp_out)]

        ##
        ## append input files
        ##
        l = os.listdir(self.fp_bams)
        for s in l:
            if s[-4:] != '.bam':
                continue
            ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--input_file
            lines += [' --input_file %s/%s \\' %(
                self.fp_bams,s,)]

        ##
        ## CommandLineGATK, optional
        ##

        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--intervals
        lines += [
            '--intervals $CHROMOSOME:$((($LSB_JOBINDEX-1)*%i+1))-$posmax \\' %(
                int(self.bps_per_interval),)]

        ##
        ## UnifiedGenotyper, optional
        ##
        ## dbSNP file. rsIDs from this file are used to populate the ID column of the output.
        lines += [' --dbsnp %s \\' %(self.fp_vcf_dbsnp)]
        ## Selecting an appropriate quality score threshold
        ## A common question is the confidence score threshold
        ## to use for variant detection. We recommend:
        ## Deep (> 10x coverage per sample) data:
        ## we recommend a minimum confidence score threshold of Q30.
        ## Shallow (< 10x coverage per sample) data:
        ## because variants have by necessity lower quality
        ## with shallower coverage we recommend
        ## a minimum confidence score of Q4 in projects
        ## with 100 samples or fewer
        ## and Q10 otherwise.
        lines += [' -stand_call_conf 4 \\'] ## ask Deepti if OK
        lines += [' -stand_emit_conf 4 \\'] ## ask Deepti if OK
        lines += [' --output_mode EMIT_VARIANTS_ONLY \\'] ## default value EMIT_VARIANTS_ONLY

        ## http://www.broadinstitute.org/gatk/gatkdocs/#VariantAnnotatorannotations
        s_annotation = ''
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_DepthOfCoverage.html
        s_annotation += ' --annotation DepthOfCoverage'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_FisherStrand.html
        s_annotation += ' -A FisherStrand'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_HaplotypeScore.html
        s_annotation += ' -A HaplotypeScore'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_MappingQualityRankSumTest.html
        s_annotation += ' -A MappingQualityRankSumTest'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_QualByDepth.html
        s_annotation += ' -A QualByDepth'
        ## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_RMSMappingQuality.html
        s_annotation += ' -A RMSMappingQuality'
        s_annotation += ' -A ReadPosRankSumTest'
        lines += [s_annotation]

        return lines


    def term_cmd(self,analysis_type,l_fp_out,):

        if type(l_fp_out) != list:
            print l_fp_out
            stop

        ## cont cmd
        lines = [';']
        for fp_out in l_fp_out:
            fp_out = fp_out
            lines += ['echo %s >> %s.touch;' %(fp_out,analysis_type,)]
            lines += ['touch touch/%s;' %(fp_out,)]
        s = '/software/bin/python-2.7.3'
        s += ' %s/GATK_pipeline2.py' %(os.path.dirname(sys.argv[0]))
        for k,v in vars(self.namespace_args).items():
            s += ' --%s %s' %(k,v)
        lines += [s]
        ## term cmd
        lines += ['"']

        lines += ['echo $cmd']
        lines += ['eval $cmd']

        return lines


    def init_UnifiedGenotyper(self,l_chroms,d_chrom_lens,):

        lines = []

        lines += ['\n## parse chromosome from command line']
        lines += ['CHROMOSOME=$1']

        lines += ['\n## define arrays']
        lines += ['CHROMOSOMES=(%s)' %(' '.join(l_chroms))]
        lines += ['LENCHROMOSOMES=(%s)' %(' '.join(
            [str(d_chrom_lens[chrom]) for chrom in l_chroms],
            ),)]

        lines += ['\n## find chromosome index in array of chromosomes']
        lines += ['for ((index=0; index<${#CHROMOSOMES[@]}; index++))\ndo']
        lines += ['if [ "${CHROMOSOMES[$index]}" = "$CHROMOSOME" ]\nthen']
        lines += ['break\nfi\ndone\n']

        lines += ['\n## find length of current chromosome passed via the command line']
        lines += ['LENCHROMOSOME=${LENCHROMOSOMES[$index]}\n']

        ##
        ## do not allow interval to exceed the length of the chromosome
        ## otherwise it will raise an I/O error (v. 1.4-15)
        ##
        lines += ['posmax=$((($LSB_JOBINDEX+0)*%i))' %(int(
            self.bps_per_interval
            ))]
        lines += ['if test $posmax -gt $LENCHROMOSOME']
        lines += ['then posmax=$LENCHROMOSOME']
        lines += ['fi\n']

        return lines


    def BeagleOutputToVCF(self,chrom,):

        ## this function is currently not called

        lines = self.init_cmd(analysis_type,memMB,)

        ## GATKwalker, required, in
        lines += [
            '--beaglePhased:BEAGLE out_BEAGLE/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.phased \\' %(
                chromosome,chromosome,
                ),
            ]
        lines += [
            '--beagleProbs:BEAGLE out_BEAGLE/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.gprobs \\' %(
                chromosome,chromosome,
                ),
            ]
        lines += [
            '--beagleR2:BEAGLE out_BEAGLE/BeagleOutput.%s.bgl.ProduceBeagleInput.%s.bgl.r2 \\' %(
                chromosome,chromosome,
                ),
            ]
        ## GATKwalker, required, out
        fp_out = 'out_GATK/%s.%s.vcf' %(analysis_type,chromosome,)
        lines += ['--out %s \\' %(fp_out,)]
        ## GATKwalker, required, in
        lines += [
            '--variant out_GATK/ApplyRecalibration.recalibrated.filtered.%s.vcf \\' %(
                chromosome
                )
            ]

        s = '\n'.join(lines)+'\n\n'

        return s


    def parse_arguments(self,):

        ## http://docs.python.org/2/library/argparse.html#module-argparse
        ## New in version 2.7
        ## optparse deprecated since version 2.7
        parser = argparse.ArgumentParser()

        ##
        ## add arguments
        ##

        parser.add_argument(
            '--bam','--bams','--bamdir','--fp_bams',
            dest='fp_bams',
            help='Path to directory containing improved BAMs',
            metavar='FILE',default=None,
            required = True,
            )

        parser.add_argument(
            '--GATK','--fp_GATK',
            dest='fp_GATK',
            help='File path to GATK (e.g. /software/varinf/releases/GATK/GenomeAnalysisTK-1.4-15-gcd43f01/GenomeAnalysisTK.jar)',
            metavar='FILE',default=None,
            required = True,
            )

        ## wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.3/b37/human_g1k_v37.fasta.gz
        parser.add_argument(
            '--FASTA','--reference','--reference-sequence','--reference_sequence','--fp_FASTA_reference_sequence',
            dest='fp_FASTA_reference_sequence',
            help='File path to reference sequence in FASTA format (e.g. /lustre/scratch111/resources/vrpipe/ref/Homo_sapiens/1000Genomes/human_g1k_v37.fasta)',
            metavar='FILE',default=None,
            required = True,
            )

        parser.add_argument(
            '--project','-P','-G',
            dest='project',
            help='Project',
            metavar='STRING',default=None,
            required = True,
            )

        parser.add_argument(
            '--arguments','--args','--options','--opts','--fp_arguments',
            dest='fp_arguments',
            metavar='FILE',default=None,
            required = False,
            )

        parser.add_argument(
            '--name','--dataset',
            dest='name',
            help='Name of dataset',
            metavar='STRING',default=None,
            required = False,
            )

        ##
        ## VariantRecalibrator resources
        ##

        parser.add_argument(
            '--resources','--VariantRecalibrator','--fp_resources',
            dest='fp_resources',
            help='File path to a file with -resource lines to append to GATK VariantRecalibrator',
            metavar='FILE',default=None,
            required = False,
            )

        ## wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.3/b37/hapmap_3.3.b37.vcf.gz
##        parser.add_argument(
##            '--hapmap',
##            dest='fp_resource_hapmap',
##            help='File path to hapmap vcf to be used by VariantCalibrator (e.g. /lustre/scratch107/projects/uganda/users/tc9/in_GATK/hapmap_3.3.b37.sites.vcf)',
##            metavar='FILE',default=None,
##            required = False,
##            )
##

        ## wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.3/b37/1000G_omni2.5.b37.vcf.gz
##        parser.add_argument(
##            '--omni',
##            dest='fp_resource_omni',
##            help='File path to omni vcf to be used by VariantCalibrator (e.g. /lustre/scratch107/projects/uganda/users/tc9/in_GATK/1000G_omni2.5.b37.sites.vcf)',
##            metavar='FILE',default=None,
##            required = False,
##            )

        ## dbSNP file. rsIDs from this file are used to populate the ID column of the output.
        ## wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.3/b37/dbsnp_137.b37.vcf.gz
        s_help = 'File path to dbsnp vcf to be used by VariantCalibrator'
        s_help += ' (e.g. /lustre/scratch107/projects/uganda/users/tc9/in_GATK/dbsnp_135.b37.vcf)'
##        s_help += '\nYou can get the vcf file with this command:'
##        s_help += '\ncurl -u gsapubftp-anonymous: ftp.broadinstitute.org/bundle/1.5/b37/dbsnp_135.b37.vcf.gz -o dbsnp_135.b37.vcf.gz; gunzip dbsnp_135.b37.vcf.gz'
        parser.add_argument(
            '--dbsnp','--fp_vcf_dbsnp',
            dest='fp_vcf_dbsnp',
            help=s_help,
            metavar='FILE',default=None,
            required = True,
            )

        ##
        ## BEAGLE
        ##
        parser.add_argument(
            '--beagle','--BEAGLE','--BEAGLEjar','--fp_software_beagle',
            dest='fp_software_beagle',
            help='File path to BEAGLE .jar file (e.g. /nfs/team149/Software/usr/share/beagle_3.3.2.jar)',
            metavar='FILE',default=None,
            required = True,
            )

        parser.add_argument(
            '--i_BEAGLE_size',
            dest='i_BEAGLE_size',
            help='Size (Mbp) of divided parts.',
            metavar='FILE',default=2,
            required = False,
            )

        parser.add_argument(
            '--i_BEAGLE_edge',
            dest='i_BEAGLE_edge',
            help='Window size (kbp) at either side of the divided part to avoid edge effects.',
            metavar='FILE',default=150,
            required = False,
            )

        s_help = 'Path to a file with 2 columns\n'
        s_help += 'Column 1: Chromosome (e.g. 22 or $CHROM)\n'
        s_help += 'Column 2: Phased file to be divided (e.g. ../in_BEAGLE/'
        s_help += 'ALL.chr22.phase1_release_v3.20101123.filt.renamed.bgl or ALL.$CHROM...'
        parser.add_argument(
            '--fp_BEAGLE_phased',
            dest='fp_BEAGLE_phased',
            help=s_help,
            metavar='FILE',default=None,
            required = True,
            )

        s_help = 'Path to a file with 2 columns\n'
        s_help += 'Column 1: Chromosome (e.g. 22 or $CHROM)\n'
        s_help += 'Column 2: Markers file to be divided (e.g. ../in_BEAGLE/'
        s_help += 'ALL.chr22.phase1_release_v3.20101123.filt.renamed.markers or ALL.$CHROM...'
        parser.add_argument(
            '--fp_BEAGLE_markers',
            dest='fp_BEAGLE_markers',
            help=s_help,
            metavar='FILE',default=None,
            required = True,
            )

        parser.add_argument(
            '--i_BEAGLE_nsamples',
            dest='i_BEAGLE_nsamples',
            help='',
            metavar='FILE',default=20,
            required = False,
            )

        ##
        ## IMPUTE2
        ##
        parser.add_argument(
            '--impute2','--IMPUTE2','--IMPUTE2jar','--fp_software_impute2',
            dest='fp_software_impute2',
            help='File path to BEAGLE .jar file (e.g. /software/hgi/impute2/v2.2.2/bin/impute2)',
            metavar='FILE',default=None,
            required = True,
            )

        parser.add_argument(
            '--hap','--impute2-hap','--fp_impute2_hap',
            dest='fp_impute2_hap',
            help='File path to directory containing hap files used by IMPUTE2 (e.g. /nfs/t149_1kg/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr$CHROMOSOME_impute.hap.gz)',
            ## Ask Deepti where this is downloaded from
            ## and whether it is always split by chromosome
            metavar='FILE',default=None,
            required = True,
            )

        parser.add_argument(
            '--legend','--impute2-legend','--fp_impute2_legend',
            dest='fp_impute2_legend',
            help='File path to directory containing legend files used by IMPUTE2 (e.g. /nfs/t149_1kg/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr$CHROMOSOME_impute.legend.gz)',
            ## Ask Deepti where this is downloaded from
            ## and whether it is always split by chromosome
            metavar='FILE',default=None,
            required = True,
            )

        parser.add_argument(
            '--map','--impute2-map','--fp_impute2_map',
            dest='fp_impute2_map',
            help='File path to directory containing map files used by IMPUTE2 (e.g. /nfs/t149_1kg/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr$CHROMOSOME_combined_b37.txt)',
            ## Ask Deepti where this is downloaded from
            ## and whether it is always split by chromosome
            metavar='FILE',default=None,
            required = True,
            )

##        ## IMPUTE2 input files
##        self.fp_impute2_map = '/nfs/t149_1kg/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr$CHROMOSOME_combined_b37.txt'

        ##
        ##
        ##
        
        ## http://docs.python.org/2/library/argparse.html#argparse.ArgumentParser.parse_args
        ## parse arguments to argparse NameSpace
        self.namespace_args = namespace_args = parser.parse_args()        

        ## http://docs.python.org/2/library/functions.html#vars
        for k,v in vars(namespace_args).items():
            setattr(self,k,v)
            continue


        if self.fp_GATK is None and self.fp_options is None:
            parser.error('--GATK or --arguments')

        bool_not_found = False
        s_arguments = ''
        for k,v in vars(namespace_args).items():
            s_arguments += '%s %s\n' %(k,v)
            if k[:3] != 'fp_':
                continue
            fp = v
            ## argument not specified
            if fp == None or fp == 'None': continue
            if fp == self.fp_bams:
                f = os.path.isdir
            else:
                f = os.path.isfile
            if not f(fp):
                print 'file path does not exist:', fp
                bool_not_found = True
        if bool_not_found == True:
            sys.exit(0)

        if self.fp_arguments == None or self.fp_arguments == 'None':
            self.fp_arguments = '%s.arguments' %(self.project)
            fd = open(self.fp_arguments,'w')
            fd.write(s_arguments)
            fd.close()
        else:
            fd = open(self.fp_arguments,'r')
            lines = fd.readlines()
            fd.close()
            for line in lines:
                l = line.strip().split()
                k = l[0]
                v = l[1]
                setattr(self,k,v)

        if self.name == None:
            self.name = os.path.basename(self.fp_bams)

        print 'dirname', os.path.dirname(self.fp_bams)
        print 'basename', os.path.basename(self.fp_bams)

        return


    def __init__(self,):

        ##
        ## parse command line arguments
        ##
        self.parse_arguments()

        ## A node on the cluster does approximately 4Mbp in 8 hours;
        ## i.e. UnifiedGenotyper, ugandan dataset, chromosome1=249Mbp=3.1weeks=504hours,
        ## so do intervals of 4Mbp to avoid cluster time out after 12 hours
        ## this should be set as a flag with optparse as it differs from dataset to dataset...
        ## this one varies from chromosome to chromosome...
        ## this variable should be dependent on the number of samples (and the coverage / number of variants... ???)
        ## todo: collect stats...
        ##
        ## the same variable is also used for IMPUTE2, so I decided to lower it from 10Mbp to 5Mbp
        self.bps_per_interval = 10.*10**6
        self.bps_per_interval_IMPUTE2 = 5.*10**6

        ## sequential order in which to run GATK walkers and other steps
        self.l_steps = [
            'UnifiedGenotyper',
            'CombineVariants',
            'VariantRecalibrator',
            'ApplyRecalibration',
            ## Imputation, BEAGLE
            'ProduceBeagleInput',
            'BEAGLE', ## not a GATK step
            ## Imputation, IMPUTE2
            'IMPUTE2', ## not a GATK step
            ]

        ## which jobs need to be run as arrays because of special memory/processor requirements?
        self.l_steps_array_subchromosome = [
            'UnifiedGenotyper',
            'IMPUTE2',
            ]
        self.l_steps_array_chromosome = [
            'ApplyRecalibration',
            'BEAGLE',
            ]

        self.d_file_outputs = {
            'ApplyRecalibration':[
                'ApplyRecalibration.recalibrated.filtered.$CHROMOSOME.vcf',
                'ApplyRecalibration.recalibrated.filtered.$CHROMOSOME.vcf.idx',
                ],
            'ProduceBeagleInput':['ProduceBeagleInput.$CHROMOSOME.bgl',],
            'BEAGLE':[
                'BeagleOutput.$CHROMOSOME.bgl.ProduceBeagleInput.$CHROMOSOME.bgl.phased',
                'BeagleOutput.$CHROMOSOME.bgl.ProduceBeagleInput.$CHROMOSOME.bgl.r2',
                'BeagleOutput.$CHROMOSOME.bgl.ProduceBeagleInput.$CHROMOSOME.bgl.gprobs',
                'BeagleOutput.$CHROMOSOME.bgl.ProduceBeagleInput.$CHROMOSOME.bgl.dose',
                'BeagleOutput.$CHROMOSOME.bgl.ProduceBeagleInput.$CHROMOSOME.bgl.r2.idx',
                'BeagleOutput.$CHROMOSOME.bgl.ProduceBeagleInput.$CHROMOSOME.bgl.gprobs.idx',
                'BeagleOutput.$CHROMOSOME.bgl.ProduceBeagleInput.$CHROMOSOME.bgl.phased.idx',
                ],
            'IMPUTE2':[
                'chromosome$CHROMOSOME.impute2.gen.$LSB_JOBINDEX_summary',
                'chromosome$CHROMOSOME.impute2.gen.$LSB_JOBINDEX',
                ],
            }

        ##
        ## binary options
        ##

        ## run everything in sequential;
        ## i.e. wait for previous steps to finish
        ## and then auto submit next step
        self.bool_sequential = True

        ## create empty shell script if file already exists
        ## should really skip the step entirely to avoid submitting 0 second job to cluster and get unpopular and acquire poor cluster stats on behalf of the ms23 group...
        self.bool_skip_if_output_exists = True

        ## not currently used, but will be implemented...
        self.bool_send_mail_upon_job_completion = True

        self.verbose = True

        return


if __name__ == '__main__':
    self = main()
    self.main()
