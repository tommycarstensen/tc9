#!/software/bin/python

## this python script should be made redundant by adding a few lines to the shell scripts...

import time, math, os
import GATK_pipeline

def main():

    instance_main = GATK_pipeline.main()
    d_chromosome_lengths = instance_main.parse_chromosome_ranges()

    hours = 0

    while hours < 72:
        bool_break = True
        for chromosome in ['X','Y',]+range(1,22+1,):
            chromosome = str(chromosome)
            intervals = int(
                math.ceil(
                    d_chromosome_lengths[chromosome]/instance_main.bps_per_interval
                    )
                )
            for array in range(1,intervals+1,):
                for suffix in ['','.idx',]:
                    fn = 'out_GATK/UnifiedGenotyper.%s.vcf.%s.idx' %(
                        chromosome,array,
                        )
                    ## file missing
                    if not os.path.isfile(fn):
                        print fn, 'missing'
                        bool_break = False
                        break
                    ## file empty
                    elif os.path.getsize(fn) == 0:
                        print fn, 'empty'
                        bool_break = False
                        break
                    else:
                        continue
                if bool_break == False:
                    break
            if bool_break == False:
                break

        if bool_break == False:
            print 'I have slept for %.2f hours, and now I will sleep some more.' %(hours)
            time.sleep(15*60)
            hours += .25
        else:
            break

    print '--------------------PROCEED TO NEXT STEP--------------------'

    return

if __name__ == '__main__':
    main()
