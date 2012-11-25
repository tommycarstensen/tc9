#!/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012


## built-ins
import os, re
## MSandhu
import sys
sys.path.append('/nfs/users/nfs_t/tc9/github/ms23/GATK_pipeline')
import GATK_pipeline
sys.path.append('/nfs/users/nfs_t/tc9/github/tc9/misc')
import gnuplot

class main:


    def main(self,):

        instance_GATK = GATK_pipeline.main()
        d_chromosome_lengths = instance_GATK.parse_chromosome_ranges()

        l_fn = os.listdir('stdout')

##        for fn in l_fn:
##            print fn
##            if (
##                fn[:len('UnifiedGenotyper')] == 'UnifiedGenotyper'
##                and
##                fn[len('UnifiedGenotyper')] != '.'
##                ):
##                old = os.path.join('stdout',fn)
##                new = os.path.join('stdout',fn.replace('UnifiedGenotyper','UnifiedGenotyper.'))
##                os.rename(old,new)
##            if '99' in fn:
##                fn_new = '.'.join([fn.split('.')[0],fn.split('.')[1],fn.split('.')[3],fn.split('.')[2],])
##                old = os.path.join('stdout',fn)
##                new = os.path.join('stdout',fn.replace(fn_new,''))
##                print old
##                print new
##                stop
##                os.rename(old,new)
##                continue
##            if (
##                fn[:len('UnifiedGenotyper')] == 'UnifiedGenotyper'
##                and
##                fn[-1] != 't'
##                ):
##                old = os.path.join('stdout',fn)
##                fn_new = '.'.join([fn.split('.')[0],fn.split('.')[1],fn.split('.')[3],fn.split('.')[2],])
##                new = os.path.join('stdout',fn_new)
##                os.rename(old,new)
##        stop

        d_resources = {'CPU':{},'Memory':{},}
        l_fn.sort()
        for fn in l_fn:
            if os.path.isdir(os.path.join('stdout',fn)):
                continue
##            print fn
            index1 = fn.index('.')
            step = fn[:index1]
            if '.' in fn[index1+1:]:
                index2 = index1+fn[index1+1:].index('.')+1
                chromosome = fn[index1+1:index2]
            else:
                chromosome = ''

            if chromosome == '23': chromosome = 'X'
            if chromosome == '24': chromosome = 'Y'

            fd = open('stdout/%s' %(fn),'r')
            lines = fd.readlines()
            fd.close()

            ## it would be faster to do rindex instead of regex,
            ## but in a few cases rubbish was appended to the farm log files
            keyword1 = re.compile(r'    Max Memory :')
            keyword2 = re.compile(r'    CPU time   :')
            l_mem = []
            l_cpu = []
            for line in lines:
                result1 = keyword1.search(line)
                result2 = keyword2.search(line)
                if result1 or result2:
                    v = float(line.split(':')[1].replace('sec.','').replace('MB',''))
                    if result1: l_mem += [v]
                    else: l_cpu += [v]

            cpu = max(l_cpu)
            mem = l_mem[l_cpu.index(cpu)]
##            if 'ApplyRecalibration' in fn and cpu > 1:
            if 'VariantRecalibrator' in fn and cpu > 1:
                print '%4i %4i %s' %(int(mem), int(cpu), chromosome), fn

            ## ignore if took less than a minute
            if cpu < 60:
                if os.path.getsize(os.path.join('stdout',fn)) < 2200:
                    print os.path.getsize(os.path.join('stdout',fn)), fn
                    stop
                    os.remove(os.path.join('stdout',fn))
                    continue
                continue

            for k_resource, v_resource in [
                ['CPU',cpu,],
                ['Memory',mem,],
                ]:
                if not step in d_resources[k_resource].keys():
                    d_resources[k_resource][step] = {}
                if not chromosome in d_resources[k_resource][step].keys():
                    d_resources[k_resource][step][chromosome] = []
                elif step not in ['UnifiedGenotyper','IMPUTE2',]:
                    print step, chromosome, k_resource, v_resource
                d_resources[k_resource][step][chromosome] += [v_resource]

        for k_resource in d_resources.keys():
            for step in d_resources[k_resource].keys():
                if 'Downsample' in step: continue
                if 'samtools' in step: continue
                l_y = []
                for chromosome in d_resources[k_resource][step].keys():
                    y = usage = d_resources[k_resource][step][chromosome]
                    if k_resource == 'CPU':
                        y = sum(y)/3600.
                    elif k_resource == 'Memory':
                        y = (sum(y)/len(y))
                    else:
                        print k_resource
                        stop
##                    lines += ['%s %s\n' %(x,y,)]
                    l_y += [y]
                if k_resource == 'CPU':
                    print k_resource, step, sum(l_y), len(l_y)
                else:
                    print k_resource, step, sum(l_y)/len(l_y), len(l_y)

        d_labels = {'Memory':'Mb','CPU':'hours'}
        for k_resource in d_resources.keys():
            for step in d_resources[k_resource].keys():
                l_x = []
                l_y = []
                if len(d_resources[k_resource][step].keys()) <= 3:
                    continue
                for chromosome in d_resources[k_resource][step].keys():
                    if chromosome == '':
                        continue
                    if chromosome[0] == '_':
                        continue
                    x = chromosome_length = d_chromosome_lengths[chromosome]/(10**6)
                    y = usage = d_resources[k_resource][step][chromosome]
                    if k_resource == 'CPU':
                        y = sum(y)/3600.
                    elif k_resource == 'Memory':
                        y = (sum(y)/len(y))
                    else:
                        print k_resource
                        stop
##                    lines += ['%s %s\n' %(x,y,)]
                    l_x += [x]
                    l_y += [y]
                if len(l_x) == 0:
                    print step
                    continue
##                fd = open('gnuplot_%s_%s.data' %(k_resource,step,),'w')
##                fd.writelines(lines)
##                fd.close()
                prefix = '%s_%s' %(k_resource,step,)
                print 'plotting', k_resource, step
                if k_resource == 'CPU':
                    print k_resource, step, sum(l_y), len(l_y)
                else:
                    print k_resource, step, sum(l_y)/len(l_y), len(l_y)
                gnuplot.scatter_plot_2d(
                    prefix,l1=l_x,l2=l_y,
                    ylabel='%s (%s)' %(k_resource,d_labels[k_resource]),
                    xlabel='chromosome length (Mbp)',
                    title=prefix.replace('_',' '),
                    )

        return


    def __init__(self,):

        self.fp_FASTA_reference_sequence = '/lustre/scratch111/resources/vrpipe/ref/Homo_sapiens/1000Genomes/human_g1k_v37.fasta'

        return


if __name__ == '__main__':
    instance_main = main()
    instance_main.main()
