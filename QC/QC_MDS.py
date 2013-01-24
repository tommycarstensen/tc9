#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import math, os, sys, inspect
sys.path.append('/nfs/users/nfs_t/tc9/github/sandbox')
import gnuplot

def main():

    d_options = parse_options()

    run_PLINK(d_options)

    plot_MDS(d_options)

    return


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


def plot_MDS(d_options):

    l_cmds = []

    bfile_in = d_options['bfile']
    bfile_out = os.path.split(bfile_in)[1]

    if not os.path.isfile('%s.mds' %(bfile_out)):
        return

    ## count number of samples
    n_samples = int(os.popen('cat %s.fam | wc -l' %(bfile_in)).read())
    if d_options['remove'] != None:
        execmd('cat %s | sort -k1,1 > remove.sorted' %(d_options['remove']))
        execmd('cat %s.fam | sort -k1,1 > fam.sorted' %(d_options['bfile']))
        execmd('join fam.sorted remove.sorted > %s.fam.joined' %(bfile_out))
        n_samples -= int(os.popen('cat %s.fam.joined | wc -l' %(bfile_out)).read())
        os.remove('remove.sorted')
        os.remove('fam.sorted')

    ## sort
    execmd('cat %s.mds | awk \'NR>1\' | sort -k1,1 > %s.mds.sorted' %(
        bfile_out,bfile_out,))
    execmd('sort -k1,1 samples2pops.dic > samples2pops.dic.sorted')
    ## join samples
    cmd = 'join -a1 -e "Unknown" -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2'
    cmd += ' %s.mds.sorted samples2pops.dic.sorted' %(bfile_out)
    cmd += ' > %s.mds.joined' %(bfile_out)
    execmd(cmd)
    if (
        int(os.popen('cat %s.mds.joined | wc -l' %(bfile_out)).read())+1
        !=
        int(os.popen('cat %s.mds | wc -l' %(bfile_out)).read())
        ):
        print int(os.popen('cat %s.mds.joined | wc -l' %(bfile_out)).read())+1
        print int(os.popen('cat %s.mds | wc -l' %(bfile_out)).read())
        stop
    ##
    cmd = 'cat %s.mds.joined' %(bfile_out)
    cmd += " | awk '{"
    cmd += 'print $1,$2,$3,$4,$5,$6,$7 > $8".mds"'
    cmd += "}'"
    execmd(cmd)
    l_pops = os.popen(
        "cat %s.mds.joined | awk '{print $8}' | sort -u" %(bfile_out)
        ).read().strip().split('\n')
    
    ## define colors for set of populations
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
        [0,0,0,],
        [85,85,85,],
        [170,170,170,],
        ]

    l_pts = [5,7,9,11]

    line_plot = 'set key out\n'
    line_plot += 'plot '
    for i in xrange(len(l_pops)):
        pop = l_pops[i]
        color = "".join(map(chr, l_colors[i%len(l_colors)])).encode('hex')
        pt = l_pts[i%len(l_pts)]
        line_plot += '"%s.mds" u 4:5 pt %i ps 2 lc rgb "#%s" t "%s", ' %(
            pop,pt,color,pop)
    line_plot = line_plot[:-2]

    c1 = 1
    c2 = 2
    gnuplot.scatter_plot_2d(
        '%s.mds' %(bfile_out),
        line_plot = line_plot,
        xlabel = 'C%i' %(c1),
        ylabel = 'C%i' %(c2),
        title='%s' %(bfile_out),
        prefix_out='mds.2D.%s.%i.%i' %(bfile_out,c1,c2),
        bool_execute = False,
        bool_remove = False,
        )

    return


def parse_options():

    if not '--bfile' in sys.argv:
        print '--bfile option missing'
        sys.exit(0)

    d_options = {
        'bfile':None,
        'remove':None,'keep':None,
        'extract':None,'exclude':None,
        }
    for option in d_options.keys():
        if not '--%s' %(option) in sys.argv:
            continue
        d_options[option] = sys.argv[sys.argv.index('--%s' %(option))+1]

    return d_options


def run_PLINK(d_options):

    bfile_in = d_options['bfile']
    bfile_out = os.path.split(bfile_in)[1]
    l_cmds = []

    ##
    ## 1) http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#freq
    ##
    cmd = 'if [ ! -s %s.frq ]; then\n\n' %(bfile_out)
    cmd += 'plink \\\n'
    cmd += '--freq \\\n'
    cmd += '--out %s \\\n' %(bfile_out,)
    for option_key,option_value in d_options.items():
        if option_value == None:
            continue
        cmd += '--%s %s \\\n' %(option_key,option_value,)
    cmd += '\nfi'
    l_cmds += [cmd]

    ##
    ## 2) http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune
    ##
    cmd = 'if [ ! -s %s.prune.in ]; then\n\n' %(bfile_out)
    cmd += 'plink \\\n'
    cmd += '--indep-pairwise 50 5 0.2 --maf 0.05 \\\n'
##    s += '--read-freq %s \\\n' %(bfile,)
    cmd += '--out %s \\\n' %(bfile_out,)
    for option_key,option_value in d_options.items():
        if option_value == None:
            continue
        cmd += '--%s %s \\\n' %(option_key,option_value,)
    cmd += '\nfi'
    l_cmds += [cmd]

    ##
    ## 3) http://pngu.mgh.harvard.edu/~purcell/plink/ibdibs.shtml#genome
    ##
    if not os.path.isfile('%s.genome' %(bfile_out)):
        l_cmds += genome(d_options)

    ##
    ## 4) http://pngu.mgh.harvard.edu/~purcell/plink/strat.shtml#cluster
    ##
    cmd = 'if [ ! -s %s.mds ]; then\n\n' %(bfile_out)
    cmd += 'plink \\\n'
    cmd += '--cluster \\\n'
    cmd += '--mds-plot 10 \\\n'
##    cmd += '--exclude ldregions.SNPs \\\n' %(bfile,)
    cmd += '--read-genome %s.genome \\\n' %(bfile_out,)
    for option_key,option_value in d_options.items():
        if option_key == 'extract':
            continue
        if option_value == None:
            continue
        cmd += '--%s %s \\\n' %(option_key,option_value,)
    cmd += '--extract %s.prune.in \\\n' %(bfile_out,)
    cmd += '--out %s \\\n' %(bfile_out,)
    cmd += '\nfi'
    l_cmds += [cmd]

    ##
    ## plot
    ##
    cmd = 'gnuplot mds.2D.%s.1.2.plt' %(bfile_out)
    l_cmds += [cmd]

    ##
    fd = open('%s.sh' %(bfile_out),'w')
    fd.write('\n\n'.join(l_cmds)+'\n\n')
    fd.close()
    execmd('chmod +x %s.sh' %(bfile_out))
    cmd = "bsub -R 'select[mem>4000] rusage[mem=4000]' -M4000000"
    cmd += ' -o %s.out -e %s.err ./%s.sh' %(bfile_out,bfile_out,bfile_out,)
    execmd(cmd)

    return


def genome(d_options):

    l_cmds = []

    bfile_in = d_options['bfile']
    bfile_out = os.path.split(bfile_in)[1]

    n_samples = int(os.popen('cat %s.fam | wc -l' %(bfile_in)).read())

    if d_options['remove'] != None:
        execmd('cat %s | sort -k1,1 > remove.sorted' %(d_options['remove']))
        execmd('cat %s.fam | sort -k1,1 > fam.sorted' %(d_options['bfile']))
        execmd('join fam.sorted remove.sorted > %s.fam.joined' %(bfile_out))
        n_samples -= int(os.popen('cat %s.fam.joined | wc -l' %(bfile_out)).read())
        execmd('join -v1 fam.sorted remove.sorted > %s.fam.joined' %(bfile_out))
        os.remove('remove.sorted')
        os.remove('fam.sorted')

    multiple = max(200,min(400,int(math.ceil(n_samples/5.))))

    if not os.path.isfile('%s.genome' %(bfile_out)):
        if d_options['remove'] != None:
            cmd = '\ncat %s.fam.joined' %(bfile_out)
        else:
            cmd = '\ncat %s.fam' %(bfile_out)
        cmd += " | awk '{print $1,$2}' "
        cmd += ' | split -d -a 2 -l %i - %s.fam.' %(
            multiple,bfile_out,
            )
        execmd(cmd)

    cmd = 'n=%i' %(n_samples)
    l_cmds += [cmd]
    ## ceiling
    cmd = 'cnt=$(((${n}+%i-1)/%i))' %(multiple, multiple,)
    l_cmds += [cmd]

    ##
    ## loop 1 - execute (and rerun)
    ##
    cmd = '\n##\n## loop 1\n##\n'

    ## begin if
    cmd += 'if [ ! -f %s.00.00.log ]; then\n'

    ## begin loop
    cmd += 'for i in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
    cmd += 'for j in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
    cmd += 'if [ $j -gt $i ]\nthen continue\nfi\n\n'
    cmd += 'if [ -s %s.$i.$j.log ]\nthen continue\nfi\n\n' %(bfile_out)
    cmd += 'if [ -s %s.genome ]\nthen continue\nfi\n\n' %(bfile_out)
    l_cmds += [cmd]

    ## init cmd
    cmd = "cmd='\n"
    ## http://pngu.mgh.harvard.edu/~purcell/plink/ibdibs.shtml#genome
    cmd += 'plink \\\n'
    cmd += '--genome \\\n'
    cmd += '--genome-lists %s.fam.$i %s.fam.$j \\\n' %(bfile_out,bfile_out,)
    cmd += '--out %s.$i.$j \\\n' %(bfile_out)
    for option_key,option_value in d_options.items():
        if option_key == 'extract':
            continue
        if option_value == None:
            continue
        cmd += '--%s %s \\\n' %(option_key,option_value,)
    cmd += '--extract %s.prune.in' %(bfile_out)

    cmd += ';'
    cmd += './%s.sh' %(bfile_out)
    ## term cmd
    cmd += "'\n"
    cmd = cmd.replace('\\\n','')
    cmd = cmd.replace('$i','$0')
    cmd = cmd.replace('$j','$1')

    cmd_LSF = '''bsub -R 'select[mem>4000] rusage[mem=4000]' -M4000000'''
    cmd += '''cmd="%s bash -c '$cmd' $i $j"\n''' %(cmd_LSF)
    cmd += 'echo $cmd\n'
    cmd += 'eval $cmd\n'

    ## end loop
    cmd += '\ndone\ndone\n'

    ## end if
    cmd += 'fi\n'
    
    l_cmds += [cmd]

    ##
    ## loop 2 - count
    ##
    cmd = '\n##\n## loop 2\n##\n'
    cmd += 'nlines=0\n'
    cmd += 'break=0\n\n'
    cmd += 'for i in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
    cmd += 'if [ $break -eq 1 ]\nthen\nbreak\nfi\n'
    cmd += 'for j in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
    cmd += 'if [ $j -gt $i ]\nthen continue\nfi\n\n'
    cmd += 'if [ ! -s %s.$i.$j.genome ]\nthen\n' %(bfile_out)
    cmd += 'break=1\nbreak\nfi\n\n'
    cmd += 'let nlines=$nlines+$('
    cmd += 'cat %s.$i.$j.genome | wc -l' %(bfile_out)
    cmd += ')-1\n'
    cmd += '\ndone\ndone\n'
    l_cmds += [cmd]

    cmd += 'echo actual $nlines expected $(($n*($n-1)/2))\n'
    cmd += 'if [ $nlines -ne $(($n*($n-1)/2)) ];then\nexit\nfi\n'
    cmd += 'if [ -f %s.genome ];then\nexit\nfi\n' %(bfile_out)
    l_cmds += [cmd]

    ##
    ## loop 3 - concatenate
    ##
    s = '                   FID1                   IID1                   FID2                   IID2 RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO'
    cmd = '\n##\n## loop 3 - concatenate genome files\n##\n'
    cmd += 'echo "%s" > %s.genome\n' %(s,bfile_out,)
    ## loop .genome files
    cmd += 'for i in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
    cmd += 'for j in $(seq -f "%02g" 0 $(($cnt-1)))\ndo\n'
    cmd += 'if [ $j -gt $i ]\nthen continue\nfi\n\n'
    ## append to .genome file
    cmd += 'cat %s.$i.$j.genome' %(bfile_out)
    cmd += " | awk 'NR>1' >> %s.genome\n" %(bfile_out)
    cmd += '\ndone\ndone\n\n'
    l_cmds += [cmd]

    return l_cmds


if __name__ == '__main__':
    main()
