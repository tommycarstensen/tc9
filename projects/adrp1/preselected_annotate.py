import pyensembl
import pysam
from natsort import natsorted
import os

ensembl = pyensembl.EnsemblRelease(release=75)


def main():

    d_variant2module = parse_modules('UKBB_modules.txt')

    d_intersect = intersect_affy6_omni25()

    path_out = 'preselected_annotated.txt'
    i = 0
    while os.path.isfile(path_out):
        i += 1
        path_out = f'preselected_annotated_{i}.txt'

    path_ref = '/lustre/scratch115/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa'

    with open('preselected.txt') as fr, open(path_out, 'w') as fw:
        ## Write single line header to output file.
        print(
            '#affyID', 'chrom', 'pos', 'ref', 'alt', '71mer',
            'rep_cnt', 'feature_cnt', 'rsID', 'gene_name(s)',
            'AF_supAFR', 'CSQ', 'module(s)', 'intersectionAffy6Omni25',
            sep='\t', file=fw,
            )
        ## Skip single line header in input file.
        for line in fr:
            break
        ## Loop over data lines of input file.
        for i, line in enumerate(fr):
            print('\n{}'.format(i), line.rstrip())
            l = line.rstrip().split('\t')
            (
                modules, affyID, chrom, pos, ref, alt, kmer,
                rep_cnt, feature_cnt,) = l[:9]
            ## Convert position from string to integer.
            pos = int(pos)

            ## Do chromosome Y from scratch with Yali.
            if chrom in ('Y', 'MT'):
                continue

            ## Parse ref and alt if missing.
            if ref == '' and alt == '':
                ref, alt = get_ref_alt_from_dbSNP(
                    chrom, pos,
                    '/lustre/scratch115/teams/sandhu/resources/ftp.ncbi.nih.gov/snp/organisms/human_9606_b149_GRCh37p13/VCF/All_20161121.norm.vcf.gz',
                    )
            assert alt[0] != '"'

            assert ref[0] in ('-', read_reference(chrom, pos, path_ref)), (ref, affyID, chrom, pos, ref, alt, read_reference(chrom, pos, path_ref))

            if rep_cnt == '':
                rep_cnt = 'NA'
            if feature_cnt == '':
                feature_cnt = 'NA'

            if not affyID.isdigit():
                assert affyID in ('NA', '',), affyID
                affyID = 'NA'
            if len(kmer) < 75:
                assert kmer in ('NA', '', 'N/A', '.', '71mer',), kmer
                kmer1, kmer2 = get_kmer(chrom, pos, ref, alt)
                kmer = kmer1
                print('kmer', kmer)
            assert kmer in get_kmer(chrom, pos, ref, alt), get_kmer(chrom, pos, ref, alt)
            ## Replace separator / with word "and".
            modules = modules.replace('Jade/Rasmus', 'Jade and Rasmus')
            modules = set(modules.split('/'))
            modules -= set(['ADME', 'Pharmacogenetics'])
            modules -= set(['Autoimmune', 'Inflammatory'])
            try:
                modules |= d_variant2module[':'.join([
                    chrom, str(pos), ref, alt])]
            except KeyError:
                pass
            modules = '; '.join(sorted(modules))
            ## https://github.com/hammerlab/pyensembl
            gene_names = parse_gene_names(chrom, pos)

            if len(l) > 9 and l[9] != 'zzzyyyxxx':
                rsID = l[9]
                print(l[9:])
                print(modules)
            else:
                rsID = None

            Continue = parse_VR(chrom, pos, ref, alt, modules)
            if Continue:
                continue

            ## Parse the rsID from dbSNP if not on the line.
            if not rsID:
                rsID = parse_dbSNP(
                    chrom, pos, ref, alt, line, modules)
                ## Skip if not in dbSNP.
                if rsID == None:
                    continue

            ## MVNcall output with triallelics not quite ready yet...
            AF_supAFR, CSQ = parse_annotations(chrom, pos)

            if pos in d_intersect[chrom]:
                intersect = True
            else:
                intersect = False

            print(
                affyID, chrom, pos, ref, alt, kmer, rep_cnt, feature_cnt,
                rsID, gene_names, AF_supAFR, CSQ, modules, intersect,
                file=fw, sep='\t',
                )

            print(i, line)

            continue

        pass

    return


def intersect_affy6_omni25():

    d_intersect = {}
    with open('intersection_Affy6_Omni2.5.txt') as f:
        for line in f:
            chrom, pos = line.rstrip().split()
            try:
                d_intersect[chrom].add(int(pos))
            except KeyError:
                d_intersect[chrom] = set([int(pos)])

    return d_intersect


def get_ref_alt_from_dbSNP(chrom, pos, path_vcf):

    tbx = pysam.TabixFile(path_vcf)
    for row in tbx.fetch(chrom, pos - 1, pos, parser=pysam.asTuple()):
        if len(row[3]) == 1 and 1 in map(len, row[4].split(',')):
            break
    else:
        stop_not_found_in_dbSNP

    assert ',' not in row[4], row

    return row[3], row[4]


def get_kmer(
    CHROM, POS, REF, ALT,
    size=71,
    path='/lustre/scratch115/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa',
    ):

    size_original = size

    ## Parse from beginning of 71mer, so subtract half the size.
    POS -= (size // 2)

    ## Insertion.
    if REF == '-':
        POS += 1
        size -= 1

    ## Deletion.
    elif ALT == '-' and len(REF) > 1:
        size += len(REF) - 1

    ## MNP
    elif len(REF) == len(ALT) and len(REF) > 1:
        size += len(REF) - 1

    ## Deletion or MNP.
    elif len(REF) > 1:
        size += len(REF) - 1

    read = read_reference(CHROM, POS, path, size=size)

    half1 = read[:size_original//2]
    half2 = read[-(size_original//2):]
#    ## Fudge solution...
#    if len(REF) > 1 and ALT == '-':
#        half1 = half1[:35]
#        half2 = half2[-35:]

#    ## Fudge assertion...
#    assert len(half1) == 35 and len(half2) == 35, (len(half1), len(half2))

    kmer1 = half1+'['+REF+'/'+ALT.replace(',','/')+']'+half2
    kmer2 = half1+'['+ALT.replace(',','/')+'/'+REF+']'+half2

    return kmer1, kmer2


def read_reference(CHROM, POS, path, size=1):

    ## This should be read from the fai file instead, but I'm lazy...
    cnt = cnt_chars_per_line_excluding_newline = 60

    ## It's a bit stupid opening the .fai each time, but I'm lazy...
    d_fai = {}
    with open(path+'.fai') as file_fai:
        for line_fai in file_fai:
            l = line_fai.rstrip().split()
            chrom = l[0]
            byte_length = int(l[1])
            byte_start = int(l[2])
            bytes_per_line_excl_line_break = int(l[3])
            bytes_per_line_incl_line_break = int(l[4])
            d_fai[chrom] = {'length': byte_length, 'start': byte_start}

    with open(path) as fd_ref:
        row1 = (POS - 1) // cnt
        row2 = (POS - 1 + size) // cnt
        bytesize = size + row2 - row1
        col = (POS - 1) % cnt
        byte_init = d_fai[CHROM]['start']
        offset = byte_init + (cnt + 1) * row1 + col
        fd_ref.seek(offset)
        read = fd_ref.read(bytesize).replace('\n', '')

    assert len(read) == size

    return read


def parse_annotations(chrom, pos):

    AF_supAFR = CSQ = 'NA'

    if chrom == 'X':
        if pos <= 2699520:
            replace = 'X_PAR1'
        elif pos >= 154931044:
            replace = 'X_PAR2'
        else:
            replace = 'X_nonPAR'
    else:
        replace = chrom
    path_vcf = '../../../SHAPEIT/out_annotate_2016Dec28/{}.minAC1.no_mask.without_related.vcf.gz'.format(replace)
    tbx = pysam.TabixFile(path_vcf)
    for row in tbx.fetch(chrom, pos - 1, pos, parser=pysam.asTuple()):
        for _ in row[7].split(';'):
            if _ == 'DB':
                continue
            k, v = _.split('=')
            if k == 'AF_supAFR':
                AF_supAFR = v
            elif k == 'CSQ':
                CSQ = v

    return AF_supAFR, CSQ


def parse_modules(path):

    d_variant2module = {}
    with open(path) as f:
        for line in f:
            module, chrom, pos, ref, alt = line.rstrip().split('\t')
            variant = ':'.join([chrom, pos, ref, alt])
            try:
                d_variant2module[variant].add(module)
            except KeyError:
                d_variant2module[variant] = set([module])

    return d_variant2module


def parse_gene_names(chrom, pos):

    gene_names = ensembl.gene_names_at_locus(contig=chrom, position=pos)
    ## https://en.wikipedia.org/wiki/Overlapping_gene
#    assert len(gene_names) <= 2, gene_names
    gene_names = '; '.join(natsorted(gene_names))

    if gene_names == '':
        gene_names = 'NA'

    return gene_names


def parse_VR(chrom, pos, ref, alt, modules):

    if chrom == 'X':
        path = '../../../pipeline_UG3.4/out_VariantRecalibrator_chromX'
    else:
        path = '../../../pipeline_UG3.4/out_VariantRecalibrator_chroms1-22'
    ## SNP or MNP
    if ref != '-' and alt != '-' and len(ref) in map(len, alt.split(',')):
        path_vcf = path + '/SNP.recal.gz'
    elif alt == '-' or len(alt) > 1 or ref == '-' or len(ref) > 1:
    #else:
        path_vcf = path + '/INDEL.recal.gz'
    else:
        print(ref, alt)
        print(modules)
        stop
    tbx = pysam.TabixFile(path_vcf)
    ## Break if found.
    for row in tbx.fetch(chrom, pos - 1, pos, parser=pysam.asTuple()):
        break
    ## Else not found.
    else:
        ## Accept that eQTLs are not present in the AGR? Discuss with Deepti.
        ## Deepti OKed this decision on 27Feb2017.
        print(
            chrom, pos, ref, alt,
            f'bcftools view -H -r {chrom}:{pos} {path_vcf}')

        ## Cystic fibrosis will for obvious reasons not be present,
        ## but don't filter it out.
        if any([
            'Cystic fibrosis' in modules,
            'Hepatitis' in modules,
            'Tuberculosis' in modules,
            'Thiopurine Methyltransferase' in modules,
            ]):
            return False
        ## SNP
        elif all([
            len(ref) == 1,
            ref != '-',
            alt != '-',
            set(map(len, alt.split(','))) == set([1]),
            ]):
            assert any([
                'Other rare coding variants' in modules,
                'Rare, possibly disease causing, mutations' in modules,
                'eQTL' in modules,
                'KIR' in modules,
                'Cardiometabolic' in modules,
                'HLA' in modules,
                ]),  modules
            ## Check if co-located indel, in which case keep the SNP.
            path_vcf = path_vcf.replace('SNP', 'INDEL')
            tbx = pysam.TabixFile(path_vcf)
            for row in tbx.fetch(chrom, pos - 1, pos, parser=pysam.asTuple()):
                with open('colocated_indel.txt', 'a') as f:
                    f.write('\t'.join(row)+'\n')
                return False
        ## Indel
        else:
            assert any([
                'Other rare coding variants' in modules,
                'Rare, possibly disease causing, mutations' in modules,
                'Protein truncating variants' in modules,
                'eQTL' in modules,
                'Neurological disorders' in modules,
                'ClinVar' in modules,
                ]), modules
            ## Check if co-located SNP, in which case keep the indel.
            path_vcf = path_vcf.replace('INDEL', 'SNP')
            tbx = pysam.TabixFile(path_vcf)
            for row in tbx.fetch(chrom, pos - 1, pos, parser=pysam.asTuple()):
                with open('colocated_SNP.txt', 'a') as f:
                    f.write('\t'.join(row)+'\n')
                return False
        ## Continue loop over preselected variants,
        ## if not found in the AGR prior to filtering.
        return True

    return False


def parse_dbSNP(chrom, pos, ref, alt, line, modules):

#    path_vcf = '/lustre/scratch114/teams/sandhu/resources/ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/All_20160601.vcf.gz'
    ## HLA and KIR were not on the alternate contigs in dbSNP142,
    ## which makes life easier.
    if 'HLA' in modules or 'KIR' in modules:
        path_vcf = '/lustre/scratch115/resources/variation/Homo_sapiens/grch37/dbsnp_142.vcf.gz'
    else:
        path_vcf = '/lustre/scratch115/teams/sandhu/resources/ftp.ncbi.nih.gov/snp/organisms/human_9606_b149_GRCh37p13/VCF/All_20161121.norm.vcf.gz'

    tbx = pysam.TabixFile(path_vcf)
    row = None
    for row in tbx.fetch(chrom, pos - 1 - 1, pos, parser=pysam.asTuple()):
        print('dbSNP', row)
        if any([
            ## SNP.
            row[3] == ref and alt in row[4].split(','),
            ## SNP in dbSNP and MNP in preselected.txt
            all([
                row[3] == ref[0],
                len(set(alt[0].split(',')) & set(row[4].split(','))) > 0,
                len(ref) in map(len, alt.split(',')),
                len(row[3]) == 1,
                ]),
            ## Insertion.
            all([
                int(row[1]) == pos,
                len(row[3]) == 1,
                ref == '-',
                ## One or more ALTs overlap (e.g. rs3835252).
                len(set(x[1:] for x in row[4].split(',')) & set(alt.split(','))) >= 1,
                ]),
            ## Deletion.
            all([
                int(row[1]) == pos,
                len(row[4]) == 1,
                alt == '-',
                row[3][:1] == ref,
                len(row[3]) > 1,
                len(row[4]) == 1,
                ]),
            ## Deletion.
            all([
                int(row[1]) + 1 == pos,
                len(row[3]) == len(ref) + 1,
                set(map(len, row[4].split(','))) == set([1]),
                alt == '-',
                row[3][1:] == ref,
                ]),
            ]):
            rsID = row[2]
            break
    ## Not found in dbSNP.
    else:
        print('\n\n')
        print(line)
        print(chrom, pos)
        with open('not_in_dbSNP.txt', 'a') as fa:
            fa.write(line)
        print('row', row[3], row[4])
        print('ref, alt', ref, alt)
        ## A lot of the neurological disorder markers seem to be erroneous.
        if modules == 'Neurological disorders':
            rsID = None
        elif all([
            'Other rare coding variants' in modules,
            ref == '-' or alt == '-' or set(map(len, alt.split(','))) != set([1]),
            ]):
            rsID = None
        elif'Rare, possibly disease causing, mutations' in modules:
            rsID = None
        else:
            print(modules)
            print('stop_not_in_dbSNP')
            exit()

    return rsID

if __name__ == '__main__':
    main()
