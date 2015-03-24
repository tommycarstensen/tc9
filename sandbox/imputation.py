#!/bin/env python3

## Tommy Carstensen
## Wellcome Trust Sanger Institute
## May-September 2013, March 2015

import os
import gzip
import subprocess
import sys
import argparse
import heapq
import re


def main():

    args = argparser()

    nsamples = int(len(gzip.open(args.hap[-1]).readline().split())/2)
    IMPUTE(args, nsamples)

    return


def IMPUTE(args, n_samples):

    queue = 'normal'

    d_index2variants, pos_max_legend = count_variants(args)

    if args.blcr:
        write_brestart()

    if args.merge_ref_panels_output_ref:
        pos_max = pos_max_legend
        d_index2variants_g = d_index2variants
    else:
        if args.gen:
            path_in = args.gen
        else:
            path_in = args.known_haps_g
        print('finding largest position in', path_in, 'in order to determine interval range')
        with gzip.open(path_in) as f:
            d_index2variants_g = {}
            for line in f:
                pos = int(line.split()[2])
                index = 1+(pos-1)//args.int_size
                try:
                    d_index2variants_g[index] += 1
                except KeyError:
                    d_index2variants_g[index] = 1
                continue
            pos_max_gen = int(line.split()[2])
##        path_strand = '{}.strand' %(args.out)
##        pos_max_haps = write_strand(path_strand)
        pos_max = max(pos_max_legend, pos_max_gen)

    for index, variants in d_index2variants.items():

        if variants == 0:
            continue
        ## Avoid this error:
        ## There are no type 2 SNPs after applying the command-line settings
        ## for this run, which makes it impossible to perform imputation.
        if not index in d_index2variants_g.keys():
            continue
        ## Skip if less than 1 type 2/3 SNP per 50kbp
        if d_index2variants_g[index] < args.int_size/50000:
            continue
        ## Skip if only 1 type 2/3 SNP
        if d_index2variants_g[index] == 1:
            continue

        path_out = os.path.join(args.out, str(index))
        os.makedirs(os.path.dirname(path_out), exist_ok=True)
        if args.merge_ref_panels_output_ref:
            ext = 'hap'
        else:
            ext = 'gen'
        ## output exists?
        if os.path.isfile('{}.{}.gz'.format(path_out, ext)):
            if os.path.getsize('{}.{}.gz'.format(path_out, ext)):
                bool_continue = True
                ## check that associated summary and warnings file was created
                for suffix in ['summary', 'warnings']:
                    path = '{}.{}'.format(path_out, suffix)
                    if not os.path.isfile(path):
                        print(path)
                        bool_continue = False
                        stop1
                        break
                if bool_continue:
                    continue

        ## memory requirements
        if not args.merge_ref_panels_output_ref:
            if 'omni' in path_in: ## tmp!!! make dependent on total count of haplotypes in panel 0 and panel 2
                intersect = 1750
                if 'p3' in ''.join(args.hap) and 'ug2g' in ''.join(args.hap):
                    slope = 17500
                else:
                    slope = 15500
                if not args.blcr:
                    queue = 'basement'
                else:
                    queue = 'long'
            else:
                intersect = 1750
                ## slope is dependent on the number of panel 2 samples
                slope = 4248 # 2184/1092 panel 0 haplotypes/samples and ~100 panel 2 samples
                slope = 11000 # 6,780 panel 0 haplotypes and ~100 panel 2 samples
        else:
            slope = 9250 # 2184/1092 and 4596/2298 panel 0 and 1 haplotypes/samples and no panel 2 samples
            slope = 15000 # 2504/5008 (1000G phase 3) and 2298/4596 (ug2g+agv) samples/haplotypes
            intersect = 2250
##        if n_samples > 1092:
##            slope *= ((n_samples/1092)**2)/8
##            intersect += 300
        memMB = intersect+int(slope*int(variants)/100000.)

        nsamples1 = int(len(gzip.open(args.hap[0]).readline().split())/2)
        nsamples2 = int(len(gzip.open(args.hap[-1]).readline().split())/2)
        print(
            args.chrom, index, 'nsamples', nsamples1, nsamples2, 'variants', variants,
            'memMB', memMB, 'intersect', intersect, 'slope', slope)

        pathLSF = 'LSF/{}'.format(args.out)
        os.makedirs(pathLSF, exist_ok=True)

        ## job running?
        if os.path.isfile('{}/{}.out'.format(pathLSF, index)):
            print('exists', '{}/{}.out'.format(pathLSF, index))
            continue

        cmd = 'bsub'
        cmd += ' -G {}'.format(args.project)
        if args.blcr:
            if queue == 'normal':
                blcrkill = 600
            elif queue == 'long':
                blcrkill = 2760
            cmd += ' -k "{} method=blcrkill {}"'.format(
                os.path.join(os.getcwd(), 'checkpoint'), blcrkill)
        cmd += " -R 'select[mem>%i] rusage[mem=%i]' -M%i" %(memMB, memMB, memMB)
        cmd += ' -q %s' %(queue)
        cmd += ' -e {}/{}.err'.format(pathLSF, index)
        cmd += ' -o {}/{}.out'.format(pathLSF, index)
##        cmd += ' -J "IMPUTE.{}.{}"'.format(args.chrom, index)
        if args.merge_ref_panels_output_ref:
            cmd += ' -J "IMPUTE.merge.{}.{}.{}.{}"'.format(intersect, slope, args.chrom, index)  ## tmp!!!
        else:
            cmd += ' -J "IMPUTE.{}.{}.{}.{}"'.format(intersect, slope, args.chrom, index)  ## tmp!!!
        if args.blcr:
            cmd += ' cr_run'

        cmd += ' {}'.format(args.impute2)
        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#required_args
        cmd += ' -m %s' %(args.map)
        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#input_options
        cmd += ' -h {}'.format(' '.join(args.hap))
        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#input_options
        cmd += ' -l {}'.format(' '.join(args.legend))

        if not args.merge_ref_panels_output_ref:
            if args.known_haps_g:
                ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-known_haps_g
                cmd += ' -known_haps_g {}'.format(args.known_haps_g)
                ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-use_prephased_g
                cmd += ' -use_prephased_g'
            else:
                cmd += ' -g {}'.format(args.gen)
#                    cmd += ' -prob_g -pgs_prob' ## only relevant if prob!=0,1
##                ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-strand_g
##                cmd += ' -strand_g %s' %(path_strand)
        else:
            ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#ex7
            ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#merging_panels
            cmd += ' -merge_ref_panels_output_ref {}'.format(path_out)

        ## BASIC OPTIONS
        ## http://mathgen.stats.ox.ac.uk/impute/basic_options.html#-int
        cmd += ' -int {:d} {:d}'.format(
            args.int_size*(index-1)+1,
            min(pos_max, args.int_size*index))
        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-buffer
        cmd += ' -buffer {}'.format(args.buffer)
##            ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-allow_large_regions
##            cmd += ' -allow_large_regions'
##            ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-include_buffer_in_output
##            cmd += ' -include_buffer_in_output'
        ## http://mathgen.stats.ox.ac.uk/impute/basic_options.html#-Ne
        cmd += ' -Ne {}'.format(args.Ne)
##            ## http://mathgen.stats.ox.ac.uk/impute/basic_options.html#-call_thresh
##            cmd += ' -call_thresh 0.9'
##            ## http://mathgen.stats.ox.ac.uk/impute/basic_options.html#-nind
##            cmd += ' -nind %i'
        if args.verbose:
            ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-verbose
            cmd += ' -verbose'

        ## MCMC OPTIONS
        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-iter
        cmd += ' -iter {}'.format(args.iter)
        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-burnin
        cmd += ' -burnin {}'.format(args.burnin)
        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-k
        cmd += ' -k {}'.format(args.k)
        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-k_hap
        cmd += ' -k_hap {}'.format(args.k_hap)

        ## FILTERING OPTIONS:
        if args.sample_g:
            cmd += ' -sample_g {}'.format(args.sample_g)
        if args.sample_g_ref:
            cmd += ' -sample_g_ref {}'.format(args.sample_g_ref)

        ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-o
        cmd += ' -o %s' %(path_out)
        if False:
            if not bool_chrX and not bool_prephased:
                ## https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-pgs
                cmd += ' -pgs' ## r2_type2 will not be calculated without it
        ## https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-o_gz
        cmd += ' -o_gz'
#            ## http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-os
#            cmd += ' -os 2'

##            cmd += ' -pgs_prob' ## not relevant when discrete values after phasing with SHAPEIT ... does SHAPEIT make use of probabilities?

        if args.chrX:
            cmd += ' -chrX'
        if args.Xpar:
            cmd += ' -Xpar'

        if not args.merge_ref_panels_output_ref:
            cmd += ' -i {}.info'.format(path_out)
        cmd += ' -r {}.summary'.format(path_out)
        cmd += ' -w {}.warnings'.format(path_out)
        if args.no_sample_qc_info:
            cmd += ' -no_sample_qc_info'

        if not args.blcr:
            subprocess.call(cmd, shell=True)
        else:
            s = subprocess.check_output(cmd, shell=True).decode()
            print(s)
            jobID = int(re.match('.*?<(.*?)>', s).group(1))
            print(jobID)
            cmd_brestart = 'bsub -G {}'.format(args.project)
            cmd_brestart += ' -o brestart.out -e brestart.err'
            cmd_brestart += ' -q small -w "ended({})"'.format(jobID)
            cmd_brestart += ' bash brestart.sh {:d} {} {:d}'.format(
                jobID, args.project, memMB)
            print(cmd_brestart)
            subprocess.call(cmd_brestart, shell=True)

    return


def write_brestart():

    ## clean up this ugly function!!!

    with open('brestart.sh', 'w') as f:
        f.write('sleep 30\n')
        ## internal field separator
        f.write("IFS=$'\\n'\n")
        ## parse input variables
        f.write('jobID=$1\n')
        f.write('project=$2\n')
        f.write('memMB=$3\necho memMB $memMB\n')
        ## set pwd
        f.write('pwd=$(pwd)\n')
        ## parse bhist
        f.write('bhist=$(bhist -l $jobID)\n')
        ## Checkpoint succeeded
        f.write('''cpsucc=$(echo $bhist | sed 's/ *//g' | grep Checkpointsucceeded | wc -l)\n''')
        ##  exit code 13, TERM_CHKPNT, Job killed after checkpointing
        f.write('''exit13=$(echo $bhist | sed 's/ *//g' | grep "Exitedwithexitcode13" | wc -l)\n''')  # could also be 13x
        ##  exit code 130, TERM_MEMLIMIT
        f.write('''exit130=$(echo $bhist | sed 's/ *//g' | grep "Exitedwithexitcode130" | wc -l)\n''')
        ##  exit code 143, SIGTERM
        f.write('''exit143=$(echo $bhist | sed 's/ *//g' | grep "Exitedwithexitcode143" | wc -l)\n''')
        ## exit code 140, run limit
        f.write(
            '''exit140=$(echo $bhist | grep TERM_RUNLIMIT | wc -l)\n''')
        ## exit code 16, pid taken
        f.write('''exit16=$(echo $bhist | sed 's/ *//g'| grep Exitedwithexitcode16 | wc -l)\n''')
        ## Checkpoint failed
        f.write('''cpfail=$(echo $bhist | sed 's/ *//g'|''')
        f.write(''' grep "Checkpointfailed" | wc -l)\n''')
        ## Done successfully
        f.write('''donesuc=$(echo $bhist | sed 's/ *//g'| grep "Donesuccessfully" | wc -l)\n''')
        ## exit if done succesfully
        f.write('if [ $donesuc -eq 1 ]; then echo $bhist >> bhist_success.txt; exit; fi\n')
        ## exit if not checkpoint succeeded and not PID taken
        f.write('if [ $exit143 -eq 0 -a $cpsucc -eq 0 -a $exit13 -eq 0')
        f.write(' -a $exit16 -eq 0 ]; then echo $bhist')
        f.write(' >> bhist_unexpectederror.txt; exit; fi\n')
        ## exit if ran out of memory
        f.write('if [ $exit130 -eq 1 ]; then echo $bhist')
        f.write(' >> bhist_TERM_MEMLIMIT.txt; exit; fi\n')
        ## restart job and capture jobID
        f.write(
            's=$(brestart -G $project -M$memMB $pwd/checkpoint/$jobID)\n')
        f.write('''jobID=$(echo $s | awk -F "[<>]" '{print $2}')\n''')
        ## report if checkpoint failed
        f.write('if [ $cpfail -ne 0 ]; then echo $s')
        f.write(' >> checkpointfailed_brestartout.txt; fi\n')
        ## be verbose
        f.write('echo s $s\n')
        f.write('echo jobID $jobID\n')
        f.write('echo memMB $memMB\n')
        ## bsub this chaperone restart script again
        f.write('bsub')
        f.write(" -R 'select[mem>'$memMB'] rusage[mem='$memMB']'")
        f.write(" -M$memMB \\\n")
        f.write(' -o brestart.out -e brestart.err \\\n')
        f.write(' -G $project -q small -w "ended($jobID)" \\\n')
        f.write(' bash brestart.sh $jobID $project $memMB\n')

    return


def argparser():

    parser = argparse.ArgumentParser()
    group_in = parser.add_mutually_exclusive_group(required=True)
    group_out = parser.add_mutually_exclusive_group(required=True)

    ## Other arguments.
    parser.add_argument('--project', required=True)
    parser.add_argument('--chrom', required=True, type=str)
    parser.add_argument('--blcr', action='store_true', default=False)
    parser.add_argument('--impute2', required=True)
    parser.add_argument('--int_size', default=2000000, type=int)

##    ## REQUIRED ARGUMENTS

    ## INPUT FILE OPTIONS
    group_in.add_argument('--gen')
    parser.add_argument('--map', required=True)
    parser.add_argument('--hap', nargs='+')
    parser.add_argument('--legend', nargs='+')

    ## OUTPUT FILE OPTIONS
    group_out.add_argument('--out')
    parser.add_argument(
        '--no_sample_qc_info', action='store_true', default=False)

    ## BASIC OPTIONS
    parser.add_argument('--buffer', default=250, type=int)
    parser.add_argument('--Ne', default=17469, type=int)
    parser.add_argument('--verbose', action='store_true', default=False)

##    ## STRAND ALIGNMENT OPTIONS

    ## FILTERING OPTIONS
    parser.add_argument('--sample_g')
    parser.add_argument('--sample_g_ref')
##    parser.add_argument('--exclude_samples_g_ref')

    ## MCMC OPTIONS
    parser.add_argument('--iter', default=30, type=int)
    parser.add_argument('--burnin', default=10, type=int)
    parser.add_argument('--k', default=80, type=int)
    parser.add_argument('--k_hap', default=500, type=int)

##    ## PRE-PHASING OPTIONS
##    parser.add_argument('--use_prephased_g')

    ## PANEL MERGING OPTIONS
    group_out.add_argument('--merge_ref_panels_output_ref')

    ## CHROMOSOME X OPTIONS
    parser.add_argument('--chrX', action='store_true', default=False)
    parser.add_argument('--Xpar', action='store_true', default=False)

    ## EXPERT OPTIONS

    args = parser.parse_args()
    if args.merge_ref_panels_output_ref and not args.out:
        args.out = args.merge_ref_panels_output_ref

    if args.merge_ref_panels_output_ref:
        assert len(args.hap) == 2 and len(args.legend) == 2
    else:
        assert len(args.hap) == 1 and len(args.legend) == 1
        if args.gen:
            assert os.path.exists(args.gen)
        else:
            assert os.path.exists(args.known_haps_g)

    assert os.path.exists(args.impute2)
    assert os.path.exists(args.map)
    for file in args.hap+args.legend:
        assert os.path.exists(file)

    if args.chrX:
        if not args.merge_ref_panels_output_ref:
            assert os.path.exists(args.sample_g)

    return args


def open_file(fn):

    if fn[-3:] == '.gz':
        print(fn)
        print(os.path.realpath(fn))
        f = gzip.open(os.path.realpath(fn),'rt')
    else:
        f = open(fn)

    return f


def write_strand(path_haps,path_strand):

    with open(path_strand, 'w') as file_strand:
        
        with open_file(path_haps) as file_haps:
            for line in file_haps:
                pos = int(line.split()[2])
                file_strand.write('%i +\n' %(pos))

    return pos


def key_generator(fd):
    for line in fd:
##        yield int(line.split()[1]), line
        l = line.split()
        yield (int(l[1]),l[2],l[3]), line


def count_variants_python35(args):

    def keyfunc(s):
        return int(s.split()[1])

    with open('file1') as fd1, open('file2') as fd2:
        fd1.readline()
        fd2.readline()
        for line in heapq.merge(fd1, fd2, key=keyfunc):
            print(line)

    return


def count_variants(args):

    print('counting variants in legend file(s)', ' '.join(args.legend), file=sys.stderr)

    ## file1 and file2 can be identical if only one file is provided
    file1 = args.legend[0]
    file2 = args.legend[-1]

    with gzip.open(file1, 'rt') as fd1, gzip.open(file2, 'rt') as fd2:
        fd1.readline()
        fd2.readline()
        it1 = key_generator(fd1)
        it2 = key_generator(fd2)
        prev_index = None
        d_index2variants = {}
        for key, line in heapq.merge(it1, it2):
            pos = key[0]
            index = 1+(pos-1)//args.int_size
            if index != prev_index:
                set_variants = set([key])
                prev_index = index
                print(index, file=sys.stderr)
            else:
                set_variants.add(key)
            d_index2variants[index] = len(set_variants)

    pos_max = pos

    return d_index2variants, pos_max


if __name__ == '__main__':
    main()
