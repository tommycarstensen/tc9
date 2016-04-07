#!/bin/env python3

import itertools
import os
import re
import urllib.request

# https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#chrX
# To phase X chromosomes, you need to specify the sex of each sample:
# A male is coded as "1"
# A female is coded as "2".
# In the SAMPLE file format, an additional column is needed to specify sex.
# It needs to be called sex or gender in the header as shown below:

# https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#family
# You can also use the GEN/SAMPLE format to specify parent-child
# as long as you add two additional columns in the SAMPLE file
# that specify the relationships between individuals.
# These two additional columns need to be called ID_father and ID_mother
# in the header line of the SAMPLE file
# in order for SHAPEIT to recognize them.


def main():

    d_ped = {}

    with open('../analysis/ped/pedigree.2016-03-15.txt') as f:
        for line in f:
            l = line.rstrip().split()
            sampleID = l[0]
            sampleIDfather = l[1]
            sampleIDmother = l[2]
            d_ped[sampleID] = (sampleIDfather, sampleIDmother)

    with open('../analysis/ped/nama_trio.2016-03-15.txt') as f:
        for line in f:
            l = line.rstrip().split()
            sampleID = l[0]
            sampleIDfather = l[2].replace('NA','0')
            sampleIDmother = l[4].replace('NA','0')
            d_ped[sampleID] = (sampleIDfather, sampleIDmother)

    d_APP2EGAN = {}
    d_EGAN2APP = {}
    d_sex = {}
    with open('../../../metadata/samples.tsv') as f:
        for line in f:
            ID_EGAN, ID_APP, sex, _ = line.rstrip().split('\t', 3)
            if ID_EGAN != 'NA':
                d_APP2EGAN[ID_APP] = ID_EGAN
                d_EGAN2APP[ID_EGAN] = ID_APP
            d_sex[ID_APP] = sex.replace('M','1').replace('F','2')
            d_sex[ID_EGAN] = d_sex[ID_APP]
    with open('../../../metadata/samples.1000g.tsv') as f:
        for line in f:
            ID, pop, superpop, sex = line.rstrip().split('\t')
            d_sex[ID] = sex.replace('male','1').replace('female','2')

    url2013 = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped'
    url2011 = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20111108_samples_pedigree/G1K_samples_20111108.ped'
    for url in (url2011, url2013,):
        with urllib.request.urlopen(url) as response:
            html = response.read().decode('utf8').split('\n')
        for line in html:
            break
        for line in html:
            if line == '':
                continue
            l = line.rstrip().split()
            IID = l[1]
            PID = l[2]
            MID = l[3]
            if PID == '0' and MID == '0':
                continue
            sex = l[4]
            d_sex[IID] = sex
            d_ped[IID] = (PID, MID,)

    with open('NHGRI_1000_Genomes_sample_list_by_population_022813.tsv') as f:
        f.readline()
        for line in f:
            if line.rstrip() == '':
                continue
            try:
                (ID, _, _, _, sex) = line.rstrip().split('\t')
            except ValueError:
                print(line)
                continue
            d_sex[ID] = sex.replace('M', '1').replace('F', '2')

    for ID_APP, (father, mother) in tuple(d_ped.items()):
        try:
            ID_EGAN = d_APP2EGAN[ID_APP]
        except KeyError:
            # not sequenced
            continue
        try:
            father = d_APP2EGAN[father]
        except KeyError:
            ## father not sequenced
            pass
        try:
            mother = d_APP2EGAN[mother]
        except KeyError:
            ## mother not sequenced
            pass
        d_ped[ID_EGAN] = (father, mother)

    d_IID2FID = cluster_samples_by_familyID(d_ped)

    for chrom in itertools.chain.from_iterable((
        range(1, 23), ('X_nonPAR',),
        )):
        print(chrom)

        for ext in ('gen', 'hap'):
            rename_add_ped(
                'out_prepareGenFromBeagle4/input.shapeit.{}.{}.sample'.format(
                    chrom, ext),
                d_ped, d_APP2EGAN, d_sex, d_EGAN2APP, d_IID2FID,
                rename=False)

        for folder in os.listdir('out_SHAPEIT_array'):
#        for folder in (
#            'Baganda_octo_quad_nodups_excrelout_sd3pc1exc_subsamp.postQC',
#            'Banyarwanda_octo_quad_excrelout_sd3pc1exc_subsamp.postQC',
#            'Barundi_subsamp.postQC',
#            'omni2.5-8_20120809_gwa_uganda_gtu_flipped.postQC',
#            'omni2.5-8_gdafpr_20141016.illuminus.nama_flipped_PARfix.postQC',
#            ):
            path = 'out_SHAPEIT_array/{}/chrom{}.sample'.format(folder, chrom)

            rename_add_ped(
                path, d_ped, d_APP2EGAN, d_sex, d_EGAN2APP, d_IID2FID)

    return


def cluster_samples_by_familyID(d_ped):

    i_FID = 1
    d_IID2FID = {}
    d_FID2IID = {}
    for ID_child, (ID_father, ID_mother) in d_ped.items():
        set_FIDs = set()
        set_IIDs = set((ID_child, ID_father, ID_mother))
        set_IIDs.discard('0')
        ## Get all FIDs this unrelated/duo/trio belongs to.
        for IID in set_IIDs:
            try:
                set_FIDs.add(d_IID2FID[IID])
            except KeyError:
                continue
        ## Create a new FID if neither child, father, mother belongs to FID.
        if len(set_FIDs) == 0:
            FID = 'FID{}'.format(i_FID)
#            FID = i_FID
            i_FID += 1
        ## Keep existing FID. Merge unrelated/duo/trio with existing FID.
        elif len(set_FIDs) == 1:
            FID = set_FIDs.pop()
            set_IIDs |= d_FID2IID[FID]
        ## Merge FIDs of previously separate unrelateds/duos/trios.
        else:
            assert ID_father != '0' or ID_mother != '0'
#            print(set_FIDs)
#            print(set_IIDs)
            ## Get all IIDs from all FIDs and assign new FID randomly.
            for FID in set_FIDs:
                set_IIDs |= d_FID2IID[FID]
                del d_FID2IID[FID]
        d_FID2IID[FID] = set_IIDs
        for IID in set_IIDs:
#            if IID == '0':
#                continue
            d_IID2FID[IID] = FID
#            try:
#                d_FID2IID[FID].add(IID)
#            except KeyError:
#                d_FID2IID[FID] set([IID])

#    print(len(set(d_IID2FID.values())))
#    print(len(d_IID2FID.keys()))
#    print(len(d_ped.keys()))

    return d_IID2FID


def rename_add_ped(
    path, d_ped, d_APP2EGAN, d_sex, d_EGAN2APP, d_IID2FID, rename=True):

    if not os.path.isfile(path):
        print('not found', path)
        return

    d_short2long = {}
    with open(path) as f:
        ## Skip two line header.
        for i in range(2):
            f.readline()
        for line in f:
            IDlong = line.rstrip().split()[0]
            IDshort = long2short(IDlong)
            d_short2long[IDshort] = IDlong

    o = ''
    with open(path) as f:
        ## write header
        for i in range(2):
            line = f.readline()
            ## 1000G OMNI chromX
            if line == 'ID_1 ID_2 missing\n':
                line = 'ID_1 ID_2 missing ID_father ID_mother sex\n'
            elif line == '0 0 0\n':
                line = '0 0 0 D D D\n'
            ## PrepareGenFromBeagle4 output.
            elif line == 'ID_1 ID_2 missing father mother\n':
                line = 'ID_1 ID_2 missing ID_father ID_mother sex\n'
            elif line == '0 0 0 D D\n':
                line = '0 0 0 D D D\n'
            ## 1000G
            elif line == 'ID_1\tID_2\tmissing\tID_father\tID_mother\n':
                line = 'ID_1 ID_2 missing ID_father ID_mother sex\n'
            elif line == '0\t0\t0\t0\t0\n':
                line = '0 0 0 D D D\n'
            ## PLINK .bed file derived data; e.g. Igbo
            elif line == 'ID_1 ID_2 missing father mother sex plink_pheno\n':
                line = 'ID_1 ID_2 missing ID_father ID_mother sex\n'
            elif line == '0 0 0 D D D B\n':
                line = '0 0 0 D D D\n'
            assert line.rstrip().split() in (
                ['ID_1', 'ID_2', 'missing', 'ID_father', 'ID_mother', 'sex'],
                ['0', '0', '0', 'D', 'D', 'D'],
                )
            o += line

        ##
        ## write body
        ##
#        d_IID2FID = {}
#        i_FID = 1
#        i_FID = max(d_IID2FID.values())+1
        for line in f:
            t = line.rstrip().split()
            if len(t) == 3:
                (ID_1, ID_2, missing) = t
                sex = ID_father = ID_mother = None
            elif len(t) == 5:
                (ID_1, ID_2, missing, ID_father, ID_mother) = t
                sex = None
            else:
                (ID_1, ID_2, missing, ID_father, ID_mother, sex, plink_pheno,) = t
            ## Convert long ID to short ID.
            IDshort = long2short(ID_2)

            try:
                ID_APP = d_EGAN2APP[IDshort]
            except KeyError:
                ID_APP = IDshort
            assert not ID_APP.startswith('EGAN')

            ID_out = None
            try:
                ID_out = d_APP2EGAN[IDshort]
            except KeyError:
                ID_out = IDshort

            ## Get sex from ID.
            if not sex:
                try:
                    sex = d_sex[IDshort]
                except KeyError:
                    try:
                        sex = d_sex[ID_2]
                    except KeyError:
                        sex = '0'
#                    if (ID_1 == 'NA' and sex == None) or ID_2 not in d_sex.keys():
#                        sex = '0'
#                    else:
#                        sex = d_sex[ID_2]

            ## Get ID of parents.
            if ID_APP in d_ped.keys():

                father = d_ped[ID_APP][0]

                mother = d_ped[ID_APP][1]

                try:
                    father = d_APP2EGAN[father]
                except KeyError:
                    father = '0'

                try:
                    mother = d_APP2EGAN[mother]
                except KeyError:
                    mother = '0'

            else:
                father = '0'
                mother = '0'

            ## Replace parental IDs if not already in the .sample file.
            if ID_father == '0' or not ID_father:
                ID_father = father
            if ID_mother == '0' or not ID_mother:
                ID_mother = mother

#            if IDshort != ID_APP:
#                print(ID_1, ID_2, IDshort, ID_APP, father, mother)

#            FID = None
#            for IID in (ID_out, ID_father, ID_mother):
#                if IID == '0':
#                    continue
#                try:
#                    FID = d_IID2FID[IID]
#                    break
#                except KeyError:
#                    continue
#            if not FID:
#                FID = 'FID{}'.format(i_FID)
#                i_FID += 1
#            for IID in (ID_out, ID_father, ID_mother):
#                d_IID2FID[IID] = FID

            ## Sample is part of pedigree.
            try:
                FID = d_IID2FID[ID_APP]
            ## Sample is not part of pedigree.
            except KeyError:
                if ID_out in d_IID2FID.keys():
                    print(ID_out)
                    stop
                FID = ID_out

            o += ' '.join((FID, ID_out, missing, ID_father, ID_mother, sex))+'\n'

#    print(o)
    prefix = path[:-7]
    with open('{}.motfatID_modified.sample'.format(prefix), 'w') as f:
        f.write(o)

    ## Try to create symlinks, so haps and sample file have same prefix,
    ## which is required by SHAPEIT2.
    try:
        os.symlink(
            '{}.haps.gz', '{}.motfatID_modified.haps.gz'.format(prefix, prefix))
    except FileExistsError:
        pass

    return


def long2short(IDlong):

    ## Check if ID has well number in it and in that case remove it.
    match = re.match(r'.*?\d{6}_[A-Z]\d\d_(.+)', IDlong)
    if match:
        IDshort = '_'.join(IDlong.split('_')[2:])
        IDshort = match.group(1)
    else:
        IDshort = IDlong

    return IDshort


if __name__ == '__main__':
    main()
