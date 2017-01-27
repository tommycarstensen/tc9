import gzip
import argparse
import fileinput

#http://online.liebertpub.com/doi/pdf/10.1089/cmb.2014.0157

#We further distinguish between three classes of errors:
#1. Flip errors, which are errors in the predicted haplotypes that can be corrected by flipping an isolated 0 to a 1 or vice versa. As an example, consider the correct haplotype pair to be 000j111, while the prediction is 010j101: the second position is then affected by a flip error.
#2. Switch errors are two consecutive SNP positions whose phases have been mistakenly predicted, and which cannot be interpreted as flip errors. For example, if the correct haplotype is 000111j111000 and the predicted haplotype is 000000j111111, then we count one switch error between positions 3 and 4.
#3. Ambiguity errors are SNP positions where flipping 0 and 1 does not lead to decreasing the (w)MEC score. One can expect an error for half of these positions, because they are, effectively, phased randomly.

#http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1180495/pdf/AJHGv73p1162.pdf
# The switch error, which measures the proportion of heterozygote positions whose phase is wrongly inferred relative to the previous heterozygote position.

def main():

    args = argparser()

    ## switch error overall
    ## switch error per sample
    ## flip errors
    ## average distance between switch errors
    ## count homs/hets
    ## count homs/hets vs missing

    (
        l_sampleIds1, l_indexes1, func1,
        l_sampleIds2, l_indexes2, func2,
        ) = parse_sampleIds(args)

    generator1 = func1(args.paths1)
    generator2 = func2(args.paths2)

    assert len(l_indexes1) > 0, 'No samples in common'
    assert len(l_indexes2) > 0, 'No samples in common'
    assert len(l_indexes1) == 1, 'I need to check that the script still works for multiple samples'

    with open(args.out, 'w') as f_out:
        (
            cnt_het, cnt_noswitch, cnt_switch, cnt_sites_intersection,
            cnt_samples_gt_diff, cnt_flip_errors, d_switch,
            ) = loop(
                args, generator1, generator2, l_indexes1, l_indexes2, f_out)

    print()
    print('cnt_het', cnt_het)
    print('cnt_noswitch', cnt_noswitch)
    print('cnt_switch', cnt_switch)
    print('cnt_sites_intersection', cnt_sites_intersection)
    print('cnt_samples_gt_diff', cnt_samples_gt_diff)
    print('cnt_flip_errors', cnt_flip_errors)

    for k in d_switch.keys():
        print(k, d_switch[k], l_sampleIds1[k])

    return


def loop(args, generator1, generator2, l_indexes1, l_indexes2, f_out):

    ## Parse chromosome, positions, ref/alt allele and genotypes.
    chrom1, pos1, a11, a12, gts1 = next(generator1)
    chrom2, pos2, a21, a22, gts2 = next(generator2)
    ## Previous GTs.
    gts_prev1 = [None]*len(gts1)
    gts_prev2 = [None]*len(gts2)
    ## Were the previous GTs flipped. Used to determine whether switch or not.
    flip_prev = [None]*len(l_indexes1)
    ## Position of previous heterozygous GT.
    pos_prev = [None]*len(l_indexes1)
    ## Was the previous position a switch error? Use to determine whether "flip error" of single.
    flip_error_prev = [None]*len(l_indexes1)
    ## Keep count of heterozygous GTs, etc.
    cnt_het = 0
    cnt_noswitch = 0
    cnt_switch = 0
    cnt_sites_intersection = 0
    cnt_samples_gt_diff = 0
    cnt_flip_errors = 0
    d_switch = {i: 0 for i in l_indexes1}

    while True:
        if chrom1 == chrom2:
            if pos1 == pos2:
                pass
            elif pos1 > pos2:
                try:
                    chrom2, pos2, a21, a22, gts2 = next(generator2)
                except StopIteration:
                    break
                continue
            else:
                try:
                    chrom1, pos1, a11, a12, gts1 = next(generator1)
                except StopIteration:
                    break
                continue
        elif args.chroms.index(chrom1) > args.chroms.index(chrom2):
            try:
                chrom2, pos2, a21, a22, gts2 = next(generator2)
            except StopIteration:
                break
            continue
        else:
            try:
                chrom1, pos1, a11, a12, gts1 = next(generator1)
            except StopIteration:
                break
            continue
        if args.verbose:
            print(chrom1, chrom2, pos1, pos2, a11, a12, a21, a22)
        cnt_sites_intersection += 1
        if all([a11 == a21, a12 == a22]):
            flip = False
        elif all([a11 == a22, a12 == a21]):
            flip = True
        else:
            if args.verbose:
                print('diff', chrom1, chrom2, pos1, pos2, a11, a12, a21, a22)
            try:
                chrom1, pos1, a11, a12, gts1 = next(generator1)
            except StopIteration:
                break
            try:
                chrom2, pos2, a21, a22, gts2 = next(generator2)
            except StopIteration:
                break
            continue
        for i, (i1, i2) in enumerate(zip(l_indexes1, l_indexes2)):
            gt1 = gts1[i1]
            gt2 = gts2[i2]
            assert gt1 in (
                ('0', '0'), ('0', '1'), ('1', '0'), ('1', '1'),
                ('./.'), ('.', '1'), ('1', '.')), 'Unexpected genotype'
            assert gt2 in (
                ('0', '0'), ('0', '1'), ('1', '0'), ('1', '1'),
                ('./.'), ('.', '1'), ('1', '.')), 'Unexpected genotype'
            ## Genotype concordance, but not heterozygous.
            if flip == False and gt1 == ('1', '1',) and gt2 == ('1', '1'):
                continue
            elif flip == True and gt1 == ('0', '0',) and gt2 == ('1', '1',):
                continue
            ## Genotype discordance.
            elif gt1 == ('0', '0',) or gt2 == ('0', '0',):
                cnt_samples_gt_diff += 1
                continue
            elif gt1 == ('1', '1',) or gt2 == ('1', '1',):
                cnt_samples_gt_diff += 1
                continue
            ## Genotype missing.
            elif gt1 == ('./.',):
                continue
            elif gt2 == ('./.',):
                continue
            ## Not phased. Separator is / instead of |.
            elif len(gt1) == 1:
                continue
            elif len(gt2) == 1:
                continue
            ## WTF???
            elif gt1 == ('1', '.',):
                continue
            elif gt2 == ('1', '.',):
                continue
            elif gt1 == ('.', '1',):
                continue
            elif gt2 == ('.', '1',):
                continue
            else:
                cnt_het += 1
                try:
                    assert any([gt1 == ('0', '1'), gt1 == ('1', '0')])
                    assert any([gt2 == ('0', '1'), gt2 == ('1', '0')])
                except:
                    print(chrom1, pos1, gt1, gt2)
                    stop
                ## First heterozygous position. Switch error not applicable.
                if flip_prev[i] == None:
                    pass
                elif flip:
                    if args.verbose:
                        print('curr', pos1, gt1, gt2, flip)
                    if any([
                        all([
                            flip_prev[i] == True,
                            gt1 == gts_prev1[i1],
                            gt2 == gts_prev2[i2],
                            ]),
                        all([
                            flip_prev[i] == True,
                            gt1[0] == gts_prev1[i1][1],
                            gt1[1] == gts_prev1[i1][0],
                            gt2[0] == gts_prev2[i2][1],
                            gt2[1] == gts_prev2[i2][0],
                            ]),
                        all([
                            flip_prev[i] == False,
                            gt1 == gts_prev1[i1],
                            gt2[0] == gts_prev2[i2][1],
                            gt2[1] == gts_prev2[i2][0],
                            ]),
                        all([
                            flip_prev[i] == False,
                            gt2 == gts_prev2[i2],
                            gt1[0] == gts_prev1[i1][1],
                            gt1[1] == gts_prev1[i1][0],
                            ]),
                        ]):
                        flip_error_prev[i] = False
                        pass
                    else:
#                        print('curr', pos1, gt1, gt2, flip)
#                        print('prev', pos_prev[i], gts_prev1[i1], gts_prev2[i2], flip_prev[i])
                        cnt_switch += 1
                        f_out.write('{}\t{}\n'.format(pos_prev[i], pos1))
                        if flip_error_prev[i]:
                            cnt_flip_errors += 2
                            print('flip', pos1, pos2, gt1, gt2)
                        flip_error_prev[i] = True
                        d_switch[i1] += 1
                        if flip_prev == False:
                            stop11b
                ## flip == False
                else:
                    if any([
                        ## Neither switching.
                        all([
                            flip_prev[i] == False,
                            gt1 == gts_prev1[i1],
                            gt2 == gts_prev2[i2],
                            ]),
                        ## Switch.
                        all([
                            flip_prev[i] == True,
                            gt1[0] == gts_prev1[i1][1],
                            gt1[1] == gts_prev1[i1][0],
                            gt2 == gts_prev2[i2],
                            ]),
                        ## Switch.
                        all([
                            flip_prev[i] == True,
                            gt1 == gts_prev1[i1],
                            gt2[0] == gts_prev2[i2][1],
                            gt2[1] == gts_prev2[i2][0],
                            ]),
                        ## Both switching.
                        all([
                            flip_prev[i] == False,
                            gt1[0] == gts_prev1[i1][1],
                            gt1[1] == gts_prev1[i1][0],
                            gt2[0] == gts_prev2[i2][1],
                            gt2[1] == gts_prev2[i2][0],
                            ]),
                        ]):
                        flip_error_prev[i] = False
                        pass
                    else:
#                        print('curr', pos1, gt1, gt2, flip)
#                        print('prev', pos_prev[i], gts_prev1[i1], gts_prev2[i2], flip_prev[i])
                        if flip_prev == True:
                            stop22
                        cnt_switch += 1
                        f_out.write('{}\t{}\n'.format(pos_prev[i], pos1))
                        d_switch[i1] += 1
                        if flip_error_prev[i]:
                            cnt_flip_errors += 2
                        print('flip', pos1, pos2, gt1, gt2)
                        flip_error_prev[i] = True
                ## Make current heterozygous genotype the previous one.
                gts_prev1[i1] = gt1
                gts_prev2[i2] = gt2
                flip_prev[i] = flip
                pos_prev[i] = pos1
#                break  # tmp!!!
        
        try:
            chrom1, pos1, a11, a12, gts1 = next(generator1)
        except StopIteration:
            break
        try:
            chrom2, pos2, a21, a22, gts2 = next(generator2)
        except StopIteration:
            break

    return (
        cnt_het, cnt_noswitch, cnt_switch, cnt_sites_intersection,
        cnt_samples_gt_diff, cnt_flip_errors, d_switch,
        )


def parse_sampleIds(args):

    d = {}
    for i, paths in enumerate((args.paths1, args.paths2)):
        if paths[0].endswith('vcf.gz'):
            func = parse_vcf
            with gzip.open(paths[0], 'rt') as f:
                for line in f:
                    if not line.startswith('##'):
                        break
            l_sampleIds = line.rstrip().split()[9:]
        else:
            func = parse_hap
            with open(paths[0][:-len('.haps.gz')]+'.sample') as f:
                for _ in range(2):
                    f.readline()
                l_sampleIds = []
                for line in f:
                    l_sampleIds.append(line.rstrip().split()[1])
        d[i] = {'sampleIds': l_sampleIds, 'func': func}

    l_indexes1 = []
    l_indexes2 = []
    for index1, sampleID1 in enumerate(d[0]['sampleIds']):
        try:
            index2 = d[1]['sampleIds'].index(sampleID1)
        except ValueError:
            continue
        l_indexes1.append(index1)
        l_indexes2.append(index2)

#    l_masks1 = [1 if i in l_indexes1 else 0 for i in range(len(d[0]['sampleIds']))]
#    l_masks2 = [1 if i in l_indexes2 else 0 for i in range(len(d[1]['sampleIds']))]

    return (
        d[0]['sampleIds'], l_indexes1, d[0]['func'],
        d[1]['sampleIds'], l_indexes2, d[1]['func'],
        )


def parse_vcf(paths):

#    with fileinput.input(files=paths, openhook=fileinput.hook_compressed) as f:
    for path in paths:
        with gzip.open(path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                l = line.rstrip().split(sep='\t')
                chrom = l[0].replace('chr','')
                pos = int(l[1])
                a1 = l[3]
                a2 = l[4]
                ## Do tuple to have it identical to the generator parse_hap.
                gts = [tuple(_.split(':')[0].split('|')) for _ in l[9:]]
                yield chrom, pos, a1, a2, gts


def parse_hap(paths):

    for path in paths:
        with gzip.open(path, 'rt') as f:
            for line in f:
                l = line.rstrip().split()
                chrom = l[0]
                pos = int(l[2])
                a1 = l[3]
                a2 = l[4]
                gts = l[5:]
                ## Do list, because zip object is not subscriptable.
                gts = list(zip(l[5::2], l[6::2]))
#                gts = (gt for (i, gt) in enumerate(zip(l[5::2], l[6::2])))
                yield chrom, pos, a1, a2, gts


def argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument('--paths1', required=True, nargs='+')
    parser.add_argument('--paths2', required=True, nargs='+')
    parser.add_argument('--out', required=True)
    parser.add_argument(
        '--chroms', nargs='+',
        default=list(map(str, range(1,23)))+['X'])
    parser.add_argument('--verbose', action='store_true')

    return parser.parse_args()


if __name__ == '__main__':
    main()
