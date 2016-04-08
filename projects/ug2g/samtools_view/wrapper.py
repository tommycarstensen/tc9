#!/software/bin/python

import os
import time

project='gdap-wgs'

d_dirs = {
'egypt':'/lustre/scratch114/projects/gdap-wgs/vrpipe-egypt/dataroot',
'fula':'/lustre/scratch114/projects/gdap-wgs/vrpipe-fula/dataroot',
'ethiopia':'/lustre/scratch114/projects/gdap-wgs/vrpipe-ethiopia/dataroot',
}

d_coverage_in = {
'uganda':4.28523, ## http://hgi-www.internal.sanger.ac.uk/qc_grind/samples_view.pl?db=hgip_vrtrack_mercury_prod&proj_id=10
'zulu':4.33162, ## http://hgi-www.internal.sanger.ac.uk/qc_grind/samples_view.pl?db=hgip_vrtrack_mercury_prod&proj_id=34
'fula':9.36184, ## http://hgi-www.internal.sanger.ac.uk/qc_grind/samples_view.pl?db=hgip_vrtrack_mercury_prod&proj_id=12
#'ethiopia':7.03864, ## http://hgi-www.internal.sanger.ac.uk/qc_grind/samples_view.pl?db=hgip_vrtrack_mercury_prod&proj_id=56
'ethiopia':7.20675,
'egypt':8.29459,
}

d_coverage_out = {
'fula':[4],
'ethiopia':[4],
'egypt':[4],
}

d_tsvs={
'ethiopia':'SEQCAP_WGS_Ethiopia_Genome_Project_low_coverage.tsv',
'fula':'SEQCAP_WGS_Low_coverage_sequencing_of_the_Fula_from_Gambia.tsv',
'egypt':'SEQCAP_WGS_EGYPT_LOW_COVERAGE.tsv',
}

memMB = 4900
memMB = 99
#path_picard = '/software/varinf/releases/PicardTools/picard-tools-1.77'
#path_GATK = '/nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar' ## 2014aug28
path_samtools = '/usr/bin/samtools'
path_samtools = '/lustre/scratch113/teams/sandhu/software/samtools/samtools'
path_pwd = os.getcwd()
path_out = path_pwd
path_out = 'out'
path_ref='/lustre/scratch111/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa'
path_java='/software/team149/opt/jre1.7.0_45/bin/java'

## loop over pops
for pop in list(d_coverage_out.keys()):
    if not os.path.isdir('%s/%s' %(path_out,pop)):
        os.mkdir('%s/%s' %(path_out,pop))

    ## loop over desired coverages
    for coverage_out in d_coverage_out[pop]:
        print(pop, coverage_out)
        os.system('touch %s.%s.in' %(coverage_out,pop,))
        os.remove('%s.%s.in' %(coverage_out,pop,))
        if not os.path.isdir('%s/%s/%s' %(path_out,pop,coverage_out)):
            os.mkdir('%s/%s/%s' %(path_out,pop,coverage_out))
        #if len(os.listdir('%s/%s' %(pop,coverage_out))) > 0:
        #    continue
        #l_fn = os.listdir(d_dirs[pop])

        cmd = '''awk 'BEGIN{FS="\\t"} FNR>1{if($8!="failed") print $4}' '''
        cmd += ' ../samtools_flagstat/%s' %(d_tsvs[pop])
        l_fn = os.popen(cmd).read().strip().split('\n')

        ## loop over bams
        for fn in l_fn:

            prefix_out = '%s/%s/%s/%s' %(path_out,pop,coverage_out,fn)

            ## input does not exist
            if not os.path.isfile(os.path.join(d_dirs[pop], fn)+'.bam'):
                continue

            ## input not indexed
            if not os.path.isfile(os.path.join(d_dirs[pop], fn)+'.bam.bai'):
                print('indexing', fn)
                os.system('bsub -G gdap-wgs -o o -e e samtools index %s.bam' %(os.path.join(d_dirs[pop],fn)))
                continue

            cmd = 'bsub -G gdap-wgs'
            cmd += ' -o o'
            cmd += ' -e e'
            cmd += ' -M%i' %(memMB)
            cmd += " -R'select[mem>%i] rusage[mem=%i]'" %(memMB,memMB)
#            cmd += ' -n2 -R"span[hosts=1]"'

            ## Picard DownsampleSam
##            cmd += ' %s -Xmx%im -Djava.io.tmpdir=tmp' %(path_java, memMB-100)
            #cmd += ' -jar %s/DownsampleSam.jar' %(path_picard)
            #cmd += ' INPUT=%s/%s' %(d_dirs[pop],fn)
            #cmd += ' OUTPUT=%s/%s' %(pop,fn)
            #cmd += ' PROBABILITY=%f' %(float(coverage_out)/d_coverage_in[pop])

            ## GATK PrintReads
##            cmd += ' %s -Xmx%im -Djava.io.tmpdir=tmp' %(path_java, memMB-100)
#            cmd += ' -jar %s' %(path_GATK)
#            cmd += ' -T PrintReads'
#            cmd += ' -o %s.bam' %(prefix_out)
#            cmd += ' -I %s/%s.bam' %(d_dirs[pop],fn)
#            cmd += ' -dfrac %f' %(float(coverage_out)/d_coverage_in[pop])
#            cmd += ' -R %s' %(path_ref)

            cmd += ' {}'.format(path_samtools)
            cmd += ' view -s {}'.format(float(coverage_out)/d_coverage_in[pop])
            cmd += ' -o %s.bam -b -h' %(prefix_out,)
            cmd += ' {}/{}.bam'.format(d_dirs[pop], fn)

            if not os.path.isfile('%s.bam' %(prefix_out)):
                print('downsampling', prefix_out)
                os.system(cmd)
                continue

            ## output not indexed
            if not os.path.isfile('%s.bam.bai' %(prefix_out)):
                print('downsampling and indexing in progress', prefix_out)
                if time.time()-os.path.getmtime('%s.bam' %(prefix_out)) > 300:
                    os.system('bsub -G {} -o index.out -e index.err {} index {}.bam'.format(project, path_samtools, prefix_out))

            ## continue loop over bams
            continue

        ## continue loop over desired coverages
        continue

    ## continue loop over pops
    continue
