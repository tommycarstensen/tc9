import pyensembl
import pysam
from natsort import natsorted

ensembl = pyensembl.EnsemblRelease(release=75)


def main():

    d_variant2module = parse_modules('UKBB_modules.txt')

    with open('preselected.txt') as fr, open('preselected_annotated.txt', 'w') as fw:
        ## Skip single line header.
        for line in fr:
            break
        for i, line in enumerate(fr):
            l = line.rstrip().split('\t')
            (
                modules, affyID, chrom, pos, ref, alt, kmer,
                rep_cnt, feature_cnt,) = l[:9]
            if affyID in ('NA', '',):
                print(affyID)
                stop1
            if kmer in ('NA', '', 'N/A', '.') or len(kmer) < 71:
                print(kmer)
                stop2
            modules = set(modules.split('/'))
            try:
                modules |= d_variant2module[':'.join([chrom, pos, ref, alt])]
            except KeyError:
                pass
            modules = '; '.join(sorted(modules))
            ## https://github.com/hammerlab/pyensembl
            pos = int(pos)
            print(i, line)
            gene_names = parse_gene_names(chrom, pos)

            if len(l) > 9:
                print(l[9:])
                print(modules)
                stop

            ## Do chromosome Y from scratch with Yali.
            if chrom in ('Y', 'MT'):
                continue

            Continue = parse_VR(chrom, pos, ref, alt, modules)
            if Continue:
                continue

            rsID, triallelic_or_indel = parse_dbSNP(
                chrom, pos, ref, alt, line, modules)
            if rsID == None:
                continue

            triallelic_or_indel = parse_AR(
                chrom, pos, triallelic_or_indel=triallelic_or_indel)

            ## MVNcall output with triallelics not quite ready yet...
            if not triallelic_or_indel:
                AF_supAFR, CSQ = parse_annotations(chrom, pos)
            else:
                AF_supAFR = CSQ = 'NA'

            print(
                affyID, chrom, pos, ref, alt, kmer, rep_cnt, feature_cnt,
                rsID, gene_names, AF_supAFR, CSQ, modules, file=fw, sep='\t',
                )

            continue

        pass

    return


def parse_AR(chrom, pos, triallelic_or_indel=False):

    path_vcf = '../../../pipeline_UG3.4_recall_union/out_ApplyRecalibration/{}.vcf.gz'.format(chrom)
    tbx = pysam.TabixFile(path_vcf)
    ## Break if found.
    for row in tbx.fetch(chrom, pos - 1, pos, parser=pysam.asTuple()):
        if len(row[4].split(',')) > 1:
            triallelic = triallelic_or_indel = True
        elif len(row[3]) > 1 or len(row[4]) > 1:
            indel = triallelic_or_indel = True
        print('**************', row[4], row[3])

    return triallelic_or_indel


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

        ## SNP
        if all([
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
        else:
            assert any([
                'Other rare coding variants' in modules,
                'Rare, possibly disease causing, mutations' in modules,
                'Protein truncating variants' in modules,
                'eQTL' in modules,
                'Neurological disorders' in modules,
                ]),  modules
            ## Check if co-located SNP, in which case keep the indel.
            path_vcf = path_vcf.replace('INDEL', 'SNP')
            tbx = pysam.TabixFile(path_vcf)
            for row in tbx.fetch(chrom, pos - 1, pos, parser=pysam.asTuple()):
                return False
        ## Continue loop over preselected variants,
        ## if not found in the AGR prior to filtering.
        return True

    return False


def parse_dbSNP(chrom, pos, ref, alt, line, modules):

    triallelic_or_indel = False
    rsID = None

    path_vcf = '/lustre/scratch115/teams/sandhu/resources/ftp.ncbi.nih.gov/snp/organisms/human_9606_b149_GRCh37p13/VCF/All_20161121.norm.vcf.gz'
    path_vcf = '/lustre/scratch114/teams/sandhu/resources/ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/All_20160601.vcf.gz'
#    path_vcf = '/lustre/scratch115/resources/variation/Homo_sapiens/grch37/dbsnp_142.vcf.gz'
    tbx = pysam.TabixFile(path_vcf)
    row = None
    for row in tbx.fetch(chrom, pos - 1 - 1, pos, parser=pysam.asTuple()):
        print(row)
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
                row[4][1:] in alt.split(','),
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
                len(row[4]) == 1,
                alt == '-',
                row[3][1:] == ref,
                ]),
            ]):
            rsID = row[2]
            if len(row[4].split(',')) > 1:
                triallelic = triallelic_or_indel = True
            elif len(row[3]) > 1 or len(row[4]) > 1:
                indel = triallelic_or_indel = True
            break
    else:
        print('\n\n')
        print(line)
        print(chrom, pos)
        with open('not_in_dbSNP.txt', 'a') as fa:
            fa.write(line)
        print('row', row[3], row[4])
        print('ref, alt', ref, alt)
        if all([
            modules == 'Neurological disorders',
            len(ref) == 1,
            len(alt) == 1,
            len(row[3]) == 1,
            len(row[4]) == 1,
            ]):
            return rsID, triallelic_or_indel
        print(modules)
#        if line.split('\t')[1] in ('89016437'):
#            return rsID, triallelic_or_indel
        stop_not_in_dbSNP

    return rsID, triallelic_or_indel

if __name__ == '__main__':
    main()
