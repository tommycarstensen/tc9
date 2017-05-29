#!/bin/env python3

import requests
import pandas as pd
import pysam
import numpy as np
import collections
import random
import re
import os
import time


def main():

    ## Define the nonPAR range.
    nonPAR = (2649521, 59034049)

    ## Create a set of positions that were called in the African samples.
    n_samples, df_AGR = read_vcf()

    ## Parse various tables with Pandas.
    df_isogg, df_affy = parse_tables()

    ## Be verbose and print some set sizes.
    print('ISOGG', len(set(df_isogg['pos'])))
    print('Affy', len(set(df_affy['pos'])))
    print('AGR', len(set(df_AGR['pos'])))

    ## Yali
    #1. Get the major branch SNPs from these two files
    #https://raw.githubusercontent.com/23andMe/yhaplo/master/input/representative.SNPs.isogg.2015tree.txt
    #https://raw.githubusercontent.com/23andMe/yhaplo/master/input/representative.SNPs.additional.txt
    #
    #Chromosome coordinates for these sites could be get from this file
    #https://raw.githubusercontent.com/23andMe/yhaplo/master/input/isogg.2016.01.04.txt
    df_step1 = df_major = step1(df_isogg)

    ## Yali
    #2. Include any SNPs (if allowed) for haplogroup A, B, E (Africa Y chromosome haplogroups) from file isogg.2016.01.04.txt <https://github.com/23andMe/yhaplo/blob/master/input/isogg.2016.01.04.txt> (excluding indels)
    #
    #2) Select at least 3 diagnostic SNPs for all the African Y haplogroups (Haplogroup A, B, E and African R1b) to the finest resolution we have for the moment according to ISOGG.
    df_step2 = step2(df_isogg, df_affy, nonPAR, df_major, df_AGR)

    ## Yali
    #3) We should also include up to 3 diagnostic SNPs for all of other non-African Y haplogroup to the resolution of the three letter level, to identify male admixture from non-African populations.
    df_step3 = step3(df_isogg, df_affy, nonPAR, df_major, df_AGR)

    ## Yali
    #4. Add all of the African R1b haplogroup SNPs identified from 1000G project below
    df_step4 = step4()

    ## Concatenate the data frames and write it all to a file.
    concatenate_and_print(df_step1, df_step2, df_step3, df_step4)

    return


def read_vcf(ACminor_min=2):

    ## Create a pysam VCF object from the chromY VCF file.
    vcf = pysam.VariantFile(
        '../../../../pipeline_UG3.4_recall_union/out_ApplyRecalibration/Y.vcf.gz')

    ## Subset the samples to African samples only.
    with open('../../../../../../metadata/samples.AFR_exclACBASW.tsv') as f:
        samples_AFR = set(vcf.header.samples) & set(f.read().splitlines())
    vcf.subset_samples(samples_AFR)

    n_samples = len(vcf.header.samples)

    ## Don't waste time reading the entire VCF each time.
    ## We just want the positions.
    if not os.path.isfile('df_AGR.csv'):
        set_AGR = set()
        print('Looping over VCF records.')
        for rec in vcf.fetch('Y'):
            ## Count genotype occurences in record.
            x = collections.Counter(
                (s['GT'] for s in rec.samples.values()))
            ## Check that minor allele count is greater than minimum.
            if n_samples - x[(None,)] - x[(0,)] < ACminor_min:
                continue
            if x[(0,)] < ACminor_min:
                continue
            set_AGR.add(rec.pos)
        df_AGR = pd.DataFrame(list(set_AGR), columns=('pos',))
        df_AGR.to_csv('df_AGR.csv', index=False)
        print('Fetched {} records.'.format(len(set_AGR)))
    else:
        df_AGR = pd.read_csv('df_AGR.csv')

    return n_samples, df_AGR


def concatenate_and_print(df_step1, df_step2, df_step3, df_step4):

    print(len(df_step1))
    print(len(df_step2))
    print(len(df_step3))
    print(len(df_step4))

    ## Concatenate data frames from steps 1 through 3.
    df_concat = pd.concat((df_step1, df_step2, df_step3)).drop_duplicates(
        subset=('pos', 'Mutation'), keep='first')
##    df_concat['ref'] = df_concat.Mutation.str[:1]
##    df_concat['alt'] = df_concat.Mutation.str[-1:]
    df_concat.loc[df_concat['RefSNP ID'] == 'rs2032623', 'Mutation'] = '-->T'
    df_concat.loc[df_concat['RefSNP ID'] == 'rs2032642', 'Mutation'] = 'T->-'
    df_concat.loc[df_concat['RefSNP ID'] == 'rs2032651', 'Mutation'] = '-->T'
    df_concat.loc[df_concat['RefSNP ID'] == 'rs2032678', 'Mutation'] = 'TTCTC->-'
    ## Write ref and alt columns based on mutation column.
    df_concat['ref'] = df_concat.Mutation.apply(lambda x: x.split('->')[0])
    df_concat['alt'] = df_concat.Mutation.apply(lambda x: x.split('->')[1])

    df_concat = df_concat.append(df_step4)

    ## Assert that there are no duplicate positions such as M35.1 and M35.2.
    assert len(df_concat) == len(df_concat['pos'].unique())

    ## Check that ref and alt are not swapped.
    path = '/lustre/scratch115/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa'
    byte_offset = get_byte_offset(path)
    for i, row in df_concat.iterrows():
        ref = read_ref(path, row.pos, byte_offset)
        if row.ref == ref or row.ref == '-' or row.alt == '-':
            continue
        assert row.alt == ref, row
        ## Switch ref and alt.
        df_concat.set_value(i, 'ref', row.alt)
        df_concat.set_value(i, 'alt', row.ref)

    df_concat.to_csv(
        'chromY_{}.tsv'.format(time.strftime("%Y%m%d-%H%M%S")),
        index=False, sep='\t')

    print('Total SNPs selected:', len(df_concat))

    return


def get_byte_offset(path):

    d_fai = {}
    with open(path + '.fai') as file_fai:
        for line_fai in file_fai:
            l = line_fai.rstrip().split()
            chrom = l[0]
            byte_length = int(l[1])
            byte_start = int(l[2])
            bytes_per_line_excl_line_break = int(l[3])
            bytes_per_line_incl_line_break = int(l[4])
            ## Map chromosome to byte offset.
            d_fai[chrom] = {'length': byte_length, 'start': byte_start}

    return d_fai['Y']['start']


def read_ref(path, pos, byte_init, chrom='Y', size=1, cnt=60):

    with open(path) as fd_ref:
        row1 = (pos - 1) // cnt
        row2 = (pos - 1 + size) // cnt
        bytesize = size + row2 - row1
        col = (pos - 1) % cnt
        byte_offset = byte_init + (cnt + 1) * row1 + col
        fd_ref.seek(byte_offset)
        read = fd_ref.read(bytesize).replace('\n', '')

    return read


def step4():

    print()

    #4a. Add all of the African R1b haplogroup SNPs identified from 1000G project below
    df_African_R1b = pd.read_csv(
        'African_R1b.tsv', sep='\t',
        names=['pos', 'ref', 'alt'],
        )
    df_African_R1b['haplogroup'] = 'R1b'

    ##4b. R1b haplogroup from Chad
    #echo "R1b1c V88 4862861 C T 1" | tr " " "\t"
    df_African_R1b = df_African_R1b.append(pd.DataFrame(
        [[4862861, 'C', 'T', 'R1b1c']],
        columns=['pos', 'ref', 'alt', 'haplogroup']))

    print('step4 SNPs', len(df_African_R1b))

    return df_African_R1b


def step2(df_isogg, df_affy, nonPAR, df_major, df_AGR):

    print()

    ##df_affy.pos = df_affy.pos.astype(str)
    df_isogg_ABE = df_isogg[
        (df_isogg.haplogroup.str.match('^[A,B,E]') == True) &
        (df_isogg['Mutation'].str.contains('indel') == False) &
        (df_isogg['Mutation'].str.len() == 4) &
        (df_isogg.haplogroup.str.contains('Investigation') == False)].drop_duplicates(
            subset=('pos', 'Mutation', 'haplogroup'), keep='first')

    set_affy = set(df_affy['pos'])
    set_AGR = set(df_AGR['pos'])

    ## Get an idea of set sizes.
    set_ABE = set(df_isogg_ABE['pos'])
    print('step2 haplogroups', len(df_isogg_ABE['haplogroup'].unique()))
    print('AGR', len(set_AGR))
    print('ABE', len(set_ABE))
    print('Affy', len(set_affy))
    print('Affy U ABE', len(set(set_affy & set_ABE) - set_AGR))
    print('Affy U AGR', len(set(set_affy & set_AGR) - set_ABE))
    print('AGR U ABE', len(set(set_AGR & set_ABE) - set_affy))
    print('ABE U Affy U AGR', len(set_AGR & set_ABE & set_affy))

    df_selected = pick_from_haplogroups(
        df_isogg_ABE, set_affy, set_AGR, df_major, nonPAR, n=5)

    return df_selected


def pick_from_haplogroups(
    df_select_from, set_affy, set_AGR, df_selected, nonPAR,
    n=3,  # Number of SNPs needed from each haplogroup.
    ):

    ## Keep count of how many SNPs from each SNP priority set are added.
    cnt00 = cnt0 = cnt1 = cnt2 = cnt3 = cnt4 = 0
    for haplogroup in df_select_from['haplogroup'].unique():
##        cnt00 = cnt0 = cnt1 = cnt2 = cnt3 = cnt4 = 0
        df_haplogroup = df_select_from.query(f'haplogroup == "{haplogroup}"')
##        print(
##            haplogroup,
##            df_haplogroup.shape[0],
##            len(df_haplogroup['pos'].unique()))
##        assert df_haplogroup.shape[0] == len(df_haplogroup['pos'].unique())
        ## Candidates are those belonging to the haplogroup
        ## and not already selected as part of the representative SNPs.
        set_haplogroup = set(df_haplogroup['pos'])
        set_selected = set(
            df_selected.query(f'haplogroup == "{haplogroup}"')['pos'])
        set_haplogroup -= set_selected
        cnt00 += len(set_selected)
        while True:
            ## Don't add anything, if the haplogroup is already represented by a sufficient number of SNPs.
            if len(set_selected) >= n:
                print(haplogroup, set_selected, set_haplogroup)
                stop
                break
##            ## Trivial case if haplogroup does not contain more SNPs (n) than needed.
##            if len(set_haplogroup | set_selected) <= n:
##                set_selected.update(set_haplogroup)
##                cnt0 += len(set_haplogroup)
##                break
            ## Priority 1: ISOGG and Affymetrix and AGR intersection.
            t = tuple(set_haplogroup & set_affy & set_AGR)
            set_ = random.sample(t, k=min(len(t), n - len(set_selected)))
            set_selected.update(set_)
            cnt1 += len(set_)
            if len(set_selected) == n:
                break
            ## Priority 2: ISOGG and Affymetrix intersection.
            t = tuple((set_haplogroup & set_affy) - set_AGR)
            set_ = random.sample(t, k=min(len(t), n - len(set_selected)))
            set_selected.update(set_)
            cnt2 += len(set_)
            if len(set_selected) == n:
                break
            ## Priority 3: ISOGG and AGR intersection.
            t = tuple((set_haplogroup & set_AGR) - set_affy)
            set_ = random.sample(t, k=min(len(t), n - len(set_selected)))
            set_selected.update(set_)
            cnt3 += len(set_)
            if len(set_selected) == n:
                break
            ## Priority 4: ISOGG complement between PAR regions.
            t = tuple(pos for pos in set(set_haplogroup - (set_affy | set_AGR)) if (pos > nonPAR[0] and pos < nonPAR[1]))
            set_ = random.sample(t, k=min(len(t), n - len(set_selected)))
            set_selected.update(set_)
            cnt4 += len(set_)
            if len(set_selected) == n:
                break
            break
##        print(
##            haplogroup, len(set_haplogroup),
##            cnt00, cnt0, cnt1, cnt2, cnt3, cnt4, sep='\t')
        assert len(set_selected) <= n
        assert len(set_selected) > 0
##        ## Do not append the markers that were already selected in step 1.
##        set_selected -= set(
##            df_selected.query(f'haplogroup == "{haplogroup}"')['pos'])
        ## Append SNPs selected from haplogroup to data frame of selected SNPs.
        df_selected = df_selected.append(
            df_haplogroup.query('pos in @set_selected'))
        ## Explicitly continue loop over haplogroups.
        continue

    ## Be a bit verbose and print some statistics.
    print(
        cnt00, cnt0, cnt1, cnt2, cnt3, cnt4,
        sum((cnt00, cnt0, cnt1, cnt2, cnt3, cnt4)),
        sum((cnt00, cnt0, cnt1, cnt2, cnt3, cnt4)) - cnt00,
        )
    print('All selected haplogroups', len(df_selected['haplogroup'].unique()))
    print(
        'Selected haplogroups this step',
        len(df_select_from['haplogroup'].unique()))
    print(
        'Haplogroups that were not selected from in this step',
        set(df_selected['haplogroup'].unique()) - set(df_select_from['haplogroup'].unique()))
##    print(set(df_selected['haplogroup'].unique()) & set(df_select_from['haplogroup'].unique()))
    print()

    return df_selected


def step1(df_isogg):

    print()

    ## Parse representative markers (SNPs and indels) of all haplogroups.
    df_repr = parse_repr()

    ## Merge all of ISOGG with the set of representative markers.
    df_major = pd.merge(
        df_isogg, df_repr, on='SNP', how='inner', suffixes=('', '_'))
    ## Select all branches from A, B, E, but not C-D, F-Z branches.
    df_major = df_major[df_major.haplogroup.str.match('^[C-D,F-Z][1-3]') == False]

    print('step1 SNPs', len(df_major), len(df_major['pos'].unique()))

    return df_major


def step3(df_isogg, df_affy, nonPAR, df_major, df_AGR):

    print()

    ## Select "to the resolution of the three letter level"
    depth_max = 3

    ##df_affy.pos = df_affy.pos.astype(str)
    df_isogg_notABE_lvl3 = df_isogg[
        (df_isogg['haplogroup'].str.match('^[A,B,E]') == False) &
##        (df_isogg.haplogroup.str.match('^R1b') == False) &
        (df_isogg['Mutation'].str.contains('indel') == False) &
        (df_isogg['Mutation'].str.len() == 4) &
        (df_isogg['haplogroup'].str.len() <= depth_max) &
        (df_isogg['haplogroup'].str.contains('Investigation') == False)].drop_duplicates(
            subset=('pos', 'Mutation'), keep='first')
    print(
        'Haplogroups being sampled from',
        list(sorted(df_isogg_notABE_lvl3['haplogroup'].unique())),
        sep='\n',
        )

    set_affy = set(df_affy['pos'])
    set_AGR = set(df_AGR['pos'])

    ## Get an idea of set sizes.
    set_notABE = set(df_isogg_notABE_lvl3['pos'])
    print('step3 haplogroups', len(df_isogg_notABE_lvl3['haplogroup'].unique()))
    print('AGR', len(set_AGR))
    print('notABE', len(set_notABE))
    print('Affy', len(set_affy))
    print('Affy U notABE', len(set(set_affy & set_notABE) - set_AGR))
    print('Affy U AGR', len(set(set_affy & set_AGR) - set_notABE))
    print('AGR U notABE', len(set(set_AGR & set_notABE) - set_affy))
    print('notABE U Affy U AGR', len(set_AGR & set_notABE & set_affy))

    df_selected = pick_from_haplogroups(
        df_isogg_notABE_lvl3, set_affy, set_AGR, df_major, nonPAR, n=3)

    return df_selected


def parse_tables():

    df_isogg = parse_isogg()

    df_affy = parse_affy()

    return df_isogg, df_affy


def parse_affy():

##    df_affy = pd.read_csv(
##        '../mtDNA_and_chromY.txt', sep='\t', header=None,
##        names=('module', 'chrom', 'pos', 'ref', 'alt')).query('chrom == "Y"')
    df_affy = pd.read_csv(
        'Candidate_Y_Mito.txt', sep='\t',
##        dtype={'position': int, 'pos': int},
        converters={'position': int},
##        names=('module', 'chrom', 'pos', 'ref', 'alt'),
##        engine='python',
        ).query('chr == "Y"').rename(columns={'position': 'pos'})

    return df_affy


def parse_repr():

    repo = 'https://raw.githubusercontent.com/23andMe/yhaplo/master/'

    url1 = repo + 'input/representative.SNPs.isogg.2015tree.txt'
    r = requests.get(url1)
    table = []
    for line in r.text.splitlines():
        l = line.strip().split()
        ## Position is wrong in ISOGG as pointed out by Anu on April 24th.
        if l[1] == 'Z14889':
            continue
        assert len(l) == 2
        for SNP in re.split('[,/]', l[1]):
            table.append([l[0], SNP])
    df1 = pd.DataFrame(table, columns=('haplogroup', 'SNP'))

    url2 = repo + 'input/representative.SNPs.additional.txt'
    df2 = pd.read_csv(
        url2, delim_whitespace=True, header=None, skipinitialspace=True,
        names=('haplogroup', 'SNP'),
        )

    df_repr = pd.concat([df1, df2], ignore_index=True)

    return df_repr


def parse_isogg():

    ## https://docs.google.com/spreadsheets/d/1jE0w48zwP3H2XV-2FBL3UBic84xQz-qBmdAl2-wJXOY/edit?ts=574acf68#gid=1934392066
    ## Get the lastest SNP content of ISOGG without full haplogroup names,
    ## because SNPs have been removed for not meeting quality guidelines.
    with open('isogg.current.txt') as f:
        set_isogg = f.read().replace('^^', '').splitlines()

    repo = 'https://raw.githubusercontent.com/23andMe/yhaplo/master/'

    url = repo + 'input/isogg.2016.01.04.txt'
    r = requests.get(url)
    table = []
    for i, line in enumerate(r.text.splitlines()):
        l = line.strip().split('\t')

        ## The first line is the header.
        if i == 0:
            header = [s.strip() for s in l]
            continue
        ## Rows with 6 colums.
        elif len(l) == 6:
            l = [s.strip() for s in l]
        ## Rows with 7 columns.
        elif len(l) == 7:
            l = [l[i].strip() for i in (0, 2, 3, 4, 5, 6)]
        ## Skip rows with neither 6 nor 7 columns.
        else:
            continue

        ## Skip if "Removed", "Withdrawn" or "Investigation",
        ## but allow "Private" and "Notes"
        if 'Removed' in l[1] or 'Investigation' in l[1] or 'Withdrawn' in l[1]:
            continue
        if l[1].endswith(' (Private)'):
            l[1] = l[1].replace(' (Private)', '')
        if l[1].endswith(' (Notes)'):
            l[1] = l[1].replace(' (Notes)', '')
        if l[1] == 'Freq. Mut. SNP in R':
            pass
        ## Skip M35.2 (C->T) and retain M35.1 (G->C),
        ## which is the more frequently occuring biallelic SNP.
        ## I could also have dropped duplicates based on position.
        if l[0] == 'M35.2':
            continue
        ## 92R7_1 has been removed from ISOGG
        ## for "not meeting quality guidelines on 22 January 2016".
        if l[0] == '92R7_1':
            continue
        if l[0] not in set_isogg:
            continue
        ## Skip if position is missing.
        if l[4].strip() == '':
            continue
        ## Convert indel ranges to single positions.
        if '..' in l[4]:
            l[4] = l[4].split('..')[0]
        if '-' in l[4]:
            l[4] = l[4].split('-')[0]
        ## Skip if multiple positions.
        if ';' in l[4]:
            continue
        l[4] = int(l[4])
        table.append(l)

    df_isogg = pd.DataFrame(table, columns=header).rename(
            columns={'Y-position (GRCh37)': 'pos', 'Haplogroup': 'haplogroup'})
    ## Change Z13537 from haplogroup H3 at position 14934195 from T->C to G->A.
    df_isogg.loc[df_isogg.SNP == 'Z13537', 'Mutation'] = 'G->A'
##    df_isogg.loc[df_isogg.SNP == 'Z16827', 'Mutation'] = 'G->A'
    ## Incorrect position discovered by Anu on April 24th.
    df_isogg.loc[df_isogg['SNP'] == 'V83', 'pos'] = 6870606

    return df_isogg


if __name__ == '__main__':
    main()
