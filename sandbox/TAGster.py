#!python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, April-May 2014

import argparse
import os
import sys
import collections
import gzip


def main():

    d_args = argparser()

    set_preselected, set_ignore = parse_input_sets(d_args)

    d_ID2tell, d_ID2cnt = count_and_write_matrix(
        d_args, set_ignore)
    del set_ignore
    len_setS = len(d_ID2tell.keys())  # used for assert in function summarize
    assert len_setS > 0

    d_cnt2IDs = invert_dictionary(d_ID2cnt)

    d_setQ, setT = prepare_sets(
        d_args, d_cnt2IDs, d_ID2cnt, d_ID2tell, set_preselected)
    del set_preselected

    setT = select_tag_SNPs(
        d_args, d_cnt2IDs, setT, d_setQ, d_ID2cnt, d_ID2tell)

    write_output(d_args, setT, d_setQ)

    summarize(setT, d_setQ, d_ID2tell, d_ID2cnt, len_setS)

    clean_up(d_args)

    return


def summarize(setT, d_setQ, d_ID2tell, d_ID2cnt, len_setS):

    print('len(setT)', len(setT))
    for i_pop, setQ in d_setQ.items():
        print('pop', i_pop)
        print('len(setQ)', len(setQ))
        print('len(setQ - setT)', len(setQ - setT))
        print('len(setQ + setT)', len(setQ.union(setT)))
    assert len(d_ID2tell.keys()) == len(d_ID2cnt.keys())
    assert len(setT) == len_setS-len(d_ID2cnt.keys())

    return


def prepare_sets(d_args, d_cnt2IDs, d_ID2cnt, d_ID2tell, set_preselected):

    d_setQ = {
        i: set() for i, _ in enumerate(d_args['in'])}  # set of tagged SNPs
    setT = set()  # set of tagging/tag SNPs
    setT = update_counts_for_preselected_tag_SNPs(
        set_preselected, d_setQ, setT, d_args,
        d_cnt2IDs, d_ID2cnt, d_ID2tell)

    return d_setQ, setT


def invert_dictionary(d_ID2cnt):

    print('invert dictionary')
    sys.stdout.flush()

    d_cnt2IDs = {}  # key = count of SNPs in LD, value = {pop:set of SNPs}
    for ID, count in d_ID2cnt.items():
        try:
            d_cnt2IDs[count].add(ID)
        except KeyError:
            d_cnt2IDs[count] = set([ID])

    return d_cnt2IDs


def write_output(d_args, setT, d_setQ):

    with open('{}.tagSNPs'.format(d_args['out']), 'w') as f:
        for ID in setT:
            f.write('{}\n'.format(ID))
    with open('{}.tagged'.format(d_args['out']), 'w') as f:
        for i_pop, setQ in d_setQ.items():
            for ID in setQ:
                f.write('{}:{}\n'.format(i_pop,ID))

    return


def parse_input_sets(d_args):

    if d_args['pretagged']:
        with open(d_args['pretagged']) as f:
            set_ignore = set(line.rstrip() for line in f)
    else:
        set_ignore = set()

    if d_args['preselected']:
        with open(d_args['preselected']) as f_preselected:
            set_preselected = set(line.rstrip() for line in f_preselected)
    else:
        set_preselected = set()

    ## tag SNPs should not be ignored.
    set_ignore -= set_preselected

    return set_preselected, set_ignore


def update_counts_for_preselected_tag_SNPs(
        set_preselected, d_setQ, setT, d_args,
        d_cnt2IDs, d_ID2cnt, d_ID2tell):

    assert len(d_ID2tell.keys()) == len(d_ID2cnt.keys())
    print(
        'SNPs in union set across populations',
        len(d_ID2tell.keys()), len(d_ID2cnt.keys()))
    with open('{}.LDmatrix'.format(d_args['out'])) as f:
        for ID in set_preselected:
            print(ID,len(d_ID2tell.keys()))
            ## ID was already tagged
#            if not ID in d_ID2cnt.keys():  # added 2014jun15
            if d_ID2cnt[ID] == 0:  # added 2014jun15
                setT.add(ID)
                del d_ID2cnt[ID]
                del d_ID2tell[ID]
                continue
            try:
                d_cnt2IDs[d_ID2cnt[ID]].remove(ID)
            except:
                print(d_ID2cnt[ID])
                stop
            setT = update_counts_and_set_of_tagged_SNPs(
                ID, d_ID2tell, f, setT, d_ID2cnt, d_cnt2IDs, d_setQ, d_args)

    assert len(set_preselected) == len(setT)

    return setT


def clean_up(d_args):

    ## Clean up large temporary files.
    os.remove('{}.LDmatrix'.format(d_args['out']))

    return


def open_file(file):

    if os.path.splitext(file)[1] == '.gz':
        f = gzip.open(file, 'rt')
    else:
        f = open(file, 'r')

    return f


def filter_lines(f, d_args, set_ignore, set_candidate):

    for line in f:

        l = line.rstrip().split()
        ## Parse SNP IDs.
        ID1 = l[0]
        ID2 = l[1]
        ## Parse MAF.
        MAF1 = float(l[5])
        MAF2 = float(l[6])
        ## Determine if SNPs are to be skipped;
        ## if below MAF threshold or set to be ignored.
        bool_ID1_skip = ID1 in set_ignore or MAF1 < d_args['min_MAF']
        bool_ID2_skip = ID2 in set_ignore or MAF2 < d_args['min_MAF']
        if bool_ID1_skip and bool_ID2_skip:
            continue
        elif bool_ID1_skip:
            set_candidate.add(ID2)
            continue
        elif bool_ID2_skip:
            set_candidate.add(ID1)
            continue
        else:
            set_candidate.add(ID1)
            set_candidate.add(ID2)
            pass

        ## Parse positions.
        pos1 = int(l[2])
        pos2 = int(l[3])
        ## Skip if positions are distant from each other.
        if abs(pos2 - pos1) > d_args['max_window']:
            continue
        ## Parse LD.
        LD = float(l[4])
        ## Skip if LD below threshold.
        if LD < d_args['min_LD']:
            continue

        yield ID1, ID2


def append_count(ID, d_ID2cnt, addend=1):

    try:
        d_ID2cnt[ID] += addend
    except KeyError:
        d_ID2cnt[ID] = addend

    return


def append_IDs(ID1, ID2, d_ID2cnt, d_LD):

    ## Append count across populations.
    for ID in set([ID1, ID2]):
        append_count(ID, d_ID2cnt)
    ## Append SNPs in LD to temporary LD dictionary.
    for IDa, IDb in ((ID1, ID2), (ID2, ID1)):
        try:
            d_LD[IDa].add(IDb)
        except KeyError:
            d_LD[IDa] = set([IDb])

    return


def count_and_write_matrix(d_args, set_ignore):

    d_ID2tell = {}  # key = SNP ID, value = pop: byte position
    d_ID2cnt = {}  # key = SNP ID, value = count of SNPs in LD
    x=0
    ## Open pseudo matrix file.
    with open('{}.LDmatrix'.format(d_args['out']), 'w') as f_matrix:
        for i_pop, file_pop in enumerate(d_args['in']):
            ## Read list of files to loop over into memory.
            with open(file_pop) as f_pop:
                files_LD = f_pop.read().rstrip().split('\n')
            ## Loop over input files.
            for file_LD in files_LD:
                print(file_LD)
                sys.stdout.flush()
                ## Declare temporary LD dictionary.
                d_LD = {}
                set_candidate = set()
                with open_file(file_LD) as f_LD:
                    for ID1, ID2 in filter_lines(
                        f_LD, d_args, set_ignore, set_candidate):
                        append_IDs(ID1, ID2, d_ID2cnt, d_LD)
                ## Add SNPs only in LD with 1) SNPs below MAF threshold or
                ## 2) ignored SNPs.
                for ID in set_candidate-set(d_ID2cnt.keys()):
                    append_count(ID, d_ID2cnt, addend=0)
                    d_LD[ID] = set()
                ## Append IDs in temporary LD dictionary to pseudo matrix file.
                for ID, set_IDs in d_LD.items():
                    ## Get pseudo matrix file byte position.
                    tell = f_matrix.tell()
                    ## Append IDs to pseudo matrix file.
                    f_matrix.write('{}\n'.format(' '.join(set_IDs)))
                    ## Append byte position to dictionary.
                    ## Does [ID][i_pop] use less memory than [i_pop][ID] ???
                    try:
                        d_ID2tell[ID][i_pop].append(tell)
                    except KeyError:  # i_pop KeyError
                        try:
                            d_ID2tell[ID][i_pop] = [tell]
                        except KeyError:  # ID KeyError
                            d_ID2tell[ID] = {i_pop: [tell]}
                ## Continue loop over input files.
                continue
            ## Continue loop over populations.
            continue
        ## Close pseudo matrix file.
        pass

    ## Delete temporary LD dictionary from memory.
    del d_LD

    return d_ID2tell, d_ID2cnt


def select_tag_SNPs(d_args, d_cnt2IDs, setT, d_setQ, d_ID2cnt, d_ID2tell):

    print('selecting tag SNPs')
    sys.stdout.flush()

    ## The number of tag SNPs cannnot exceed the number of input SNPs.
    d_args['max_tagSNP'] = min(
        d_args['max_tagSNP'],
        len(d_ID2tell.keys())+len(setT))

    d_len_setS = {i_pop: 0 for i_pop, _ in enumerate(d_args['in'])}
    for ID in d_ID2tell.keys():
        for i_pop in d_ID2tell[ID].keys():
            d_len_setS[i_pop] += 1

    max_coverage = d_args['max_coverage']

    setS0 = set([x for x in d_ID2cnt.keys()])

    count_tagging = most_common_key = max(d_cnt2IDs.keys())
    with open('{}.LDmatrix'.format(d_args['out'])) as f_matrix:
        while len(setT) < d_args['max_tagSNP']:
            ## Find the most common element.
            ## collections.Counter.most_common(1) is too slow.
            while True:
                try:
                    ID = d_cnt2IDs[most_common_key].pop()
                    break
                ## Dictionary d_cnt2IDs exhausted.
                except KeyError:
                    del d_cnt2IDs[most_common_key]
                    try:
                        most_common_key = max(d_cnt2IDs.keys())
                        print('max tag', most_common_key, flush=True)
                    except ValueError:
                        ## Dictionary d_cnt2IDs exhausted.
                        if not d_cnt2IDs:
                            return setT
                        stopshouldnothappen1
                        break
            setT = update_counts_and_set_of_tagged_SNPs(
                ID, d_ID2tell, f_matrix, setT, d_ID2cnt, d_cnt2IDs,
                d_setQ, d_args)
            ## Continue while loop.
            continue

    return setT


def update_counts_and_set_of_tagged_SNPs(
        ID_tagging, d_ID2tell, f, setT, d_ID2cnt, d_cnt2IDs, d_setQ, d_args):

    ## File seek is the second slowest step in this loop.

    for i_pop in d_ID2tell[ID_tagging].keys():
        ## Loop over file position(s) of the tagging SNP.
        for seek_tagging in d_ID2tell[ID_tagging][i_pop]:
            f.seek(seek_tagging)
            ## Loop over tagged SNPs.
            l_IDs_tagged = f.readline().rstrip().split()
            for ID_tagged in l_IDs_tagged:
                ## already a tagged SNP
                if ID_tagged in d_setQ[i_pop]:
                    continue
                ## already a tagging SNP
                if ID_tagged in setT:
                    print(ID_tagged, ID_tagging,l_IDs_tagged)
                    stop_not_expected_should_already_be_in_setQ
                    continue
                d_setQ[i_pop].add(ID_tagged)
                subtract_main(
                    ID_tagged, i_pop, d_ID2tell, f, ID_tagging, d_ID2cnt,
                    d_cnt2IDs)
                ## Continue loop over tagged SNPs.
                continue
            ## Continue loop over file position(s) of the tagging SNPs.
            continue
        ## Add tagging SNP to set of tagged SNPs.
        if not ID_tagging in d_setQ[i_pop]:
            d_setQ[i_pop].add(ID_tagging)
            subtract_main(
                ID_tagging, i_pop, d_ID2tell, f, ID_tagging, d_ID2cnt,
                d_cnt2IDs)
        ## Continue loop over populations.
        continue

    del d_ID2tell[ID_tagging]
    setT.add(ID_tagging)
    del d_ID2cnt[ID_tagging]

    return setT


def subtract_main(
        ID_tagged, i_pop, d_ID2tell, f, ID_tagging, d_ID2cnt, d_cnt2IDs):

    ## CHECK AFTER DEADLINE WHETHER FASTER TO CREATE LOCAL
    ## DICTIONARY WITH COUNTS AND THEN SUBTRACT COUNTS FROM MAIN
    ## DIC AFTER LOOP
    ## Loop over file position(s) of the tagged SNP.
    for seek_tagged in d_ID2tell[ID_tagged][i_pop]:
        f.seek(seek_tagged)
        ## Loop over the SNPs in LD with the tagged SNP.
        l_IDs_LD = f.readline().rstrip().split()
        for ID_LD in l_IDs_LD:
            ## Skip as ID_tagging is deleted from d_ID2cnt anyway.
            if ID_LD == ID_tagging:
                continue
##            if ID_LD == ID_tagged:
##                print(ID_LD, ID_tagged, i_pop)
##                stop_not_expected
            ## "if s_im not in Q_i"
            subtract(ID_LD, d_ID2cnt, d_cnt2IDs, 1, d_ID2tell)
            ## Continue loop over the SNPs
            ## in LD with the tagged SNP.
            continue
        ## Continue loop over file position(s) of the tagged SNP.
        continue

    return


def subtract(ID, d_ID2cnt, d_cnt2IDs, subtrahend, d_ID2tell):

    count = d_ID2cnt[ID]

    if count-subtrahend < 0:  # tmp!!!
        print(ID, d_ID2cnt[ID])  # tmp!!!
        stopshouldnothappen  # tmp!!!

    d_cnt2IDs[count].remove(ID)
    if count-subtrahend > 0:
        try:
            d_cnt2IDs[count-subtrahend].add(ID)
        except KeyError:
            d_cnt2IDs[count-subtrahend] = set([ID])
        d_ID2cnt[ID] -= subtrahend
    else:
        d_ID2cnt[ID] = 0  # added 2014jun15

    return


def argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument('--in', required=True, nargs='+')
    parser.add_argument('--out', required=True)
    parser.add_argument('--min_LD', type=float, default=.8)
    parser.add_argument('--min_MAF', type=float, default=.01)
    parser.add_argument('--max_window', type=int, default=250000)
    parser.add_argument('--preselected')
    parser.add_argument('--pretagged', '--ignore')
    parser.add_argument('--max_tagSNP', type=int, default=1000000)
    parser.add_argument('--max_coverage', type=float, default=1.)

    d_args = vars(parser.parse_args())

    if d_args['preselected'] == 'None':
        d_args['preselected'] = None
    if d_args['pretagged'] == 'None':
        d_args['pretagged'] = None

    return d_args


if __name__ == '__main__':
    main()
