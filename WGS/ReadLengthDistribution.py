import os

##for chromosome in range(1,22+1,)+['X','Y',]:
##    fp = 'out_GATK/ReadLengthDistribution%s.tbl' %(chromosome)
##    if os.path.isfile(fp):
##        continue
##    s = 'bsub \
##    -o tmp.out -e tmp.err \
##    -M1000000 -R\'select[mem>1000] rusage[mem=1000]\' \
##    java -Xmx1g \
##    -jar /software/varinf/releases/GATK/GenomeAnalysisTK-1.4-15-gcd43f01/GenomeAnalysisTK.jar \
##    -T ReadLengthDistribution \
##    -I /lustre/scratch111/projects/uganda/release/20120610/chromosome_bams/chrom%s.bam \
##    -R /lustre/scratch111/resources/vrpipe/ref/Homo_sapiens/1000Genomes/human_g1k_v37.fasta \
##    -o %s \
##    ' %(chromosome, fp,)
##    os.system(s)

l_samples = []
d_read_lengths_normalized = {}
for chromosome in range(1,22+1,)+['X','Y',]:
    fp = 'out_GATK/ReadLengthDistribution%s.tbl' %(chromosome)
    if not os.path.isfile(fp):
        continue
    if chromosome in ['X','Y',]:
        continue
    fd = open(fp,'r')
    lines = fd.readlines()
    fd.close()
    l_samples = lines[1].split()[1:]
    l_read_lengths = [int(read_length) for read_length in lines[-2].split()[1:]]
    read_length_average = float(sum(l_read_lengths))/len(l_read_lengths)
    for i in range(100):
        read_length = l_read_lengths[i]
        sample = l_samples[i]
        ## normalize read length across samples
        read_length_normalized = read_length/read_length_average
        if read_length_normalized < 0.55 or read_length_normalized > 1.55:
            print chromosome, sample, read_length_normalized, read_length
        if not sample in d_read_lengths_normalized.keys():
            d_read_lengths_normalized[sample] = {}
        d_read_lengths_normalized[sample][chromosome] = [read_length_normalized,read_length,read_length_average,]

for sample in d_read_lengths_normalized.keys():
    l_read_lengths_normalized = [l[0] for l in d_read_lengths_normalized[sample].values()]
    ## normalize read length across chromosomes
    read_length_normalized_average = sum(l_read_lengths_normalized)/len(l_read_lengths_normalized)
    for chromosome,l in d_read_lengths_normalized[sample].items():
        read_length_normalized = l[0]
        read_length = l[1]
        ## average read length for a chromosome across samples
        read_length_average = l[2]
        read_length_normalized_normalized = read_length_normalized/read_length_normalized_average
        if read_length_normalized_normalized < 0.9 or read_length_normalized_normalized > 1.1:
            print sample, chromosome,
##            print round(read_length_normalized_average,2),
##            print round(read_length_normalized,2), round(read_length_normalized_normalized,2),
            print read_length, '%7i' %(int(read_length_average))
##            print '%5.1f%%' %(round(100*read_length_normalized_normalized,0))
