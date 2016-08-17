#!/bin/python3

## Tommy Carstensen
## Wellcome Trust Sanger Institute
## 2016Aug16

import argparse
import contextlib
import gzip
import os


def main():

    args = argparser()

    d_fai = read_fai(args.ref + '.fai')

    l_paths = list(sorted(args.i))

    for path in l_paths:
        if not os.path.isfile('.'.join(path.split('.')[:-2])+'.sample'):
            print('.'.join(path.split('.')[:-2])+'.sample', 'does not exist')
            exit()

    with contextlib.ExitStack() as stack:
        d_files = {}
        for path in l_paths:
            d_files[path] = stack.enter_context(gzip.open(path, 'rt'))
        fd_out = stack.enter_context(gzip.open(args.out, 'wt'))
        fd_ref = stack.enter_context(open(args.ref))
        d_haps2pos = {path: -1 for path in l_paths}
        pos_max = 0
        d_lines = {}
        while True:
            for path in l_paths:
                while d_haps2pos[path] < pos_max:
                    line = d_files[path].readline()
                    if line == '':
                        break
                    l = line.rstrip().split(' ', maxsplit=5)
                    pos = int(line.split()[2])
                    d_haps2pos[path] = pos
                    d_lines[path] = l
                if line == '':
                    break
                pos_max = max(pos_max, pos)
            if line == '':
                break
            if min(d_haps2pos.values()) < pos_max:
                continue
            ## Make sure 2nd while loop is not infinite.
            pos_max += 1

            for path in l_paths:
                A1 = d_lines[path][3]
                A2 = d_lines[path][4]
                ID = d_lines[path][1]
                chrom = d_lines[path][0]
                if A1 != '0' and A2 != '0':
                    break
            assert A1 != '0' and A2 != '0'

            ## Parse reference allele from reference sequence.
            REF = parse_ref(d_fai, fd_ref, chrom, pos)

            if REF == A1:
                line_out = ' '.join((chrom, ID, str(pos), A1, A2))
            else:
                assert A2 == REF
                line_out = ' '.join((chrom, ID, str(pos), A2, A1))
            for path in l_paths:
                if any([
                    d_lines[path][3] == REF and (d_lines[path][4] == A2 or d_lines[path][4] == '0'),
                    d_lines[path][3] == REF and d_lines[path][4] == REF,
                    ]):
                    if d_lines[path][3] == d_lines[path][4]:
                        assert len(d_lines[path][-1])+1 == 2*d_lines[path][-1].count('0')
                    line_out += ' '+d_lines[path][-1]
                    pass
                elif any([
                    d_lines[path][3] == A2 and (d_lines[path][4] == A1 or d_lines[path][4] == '0' or d_lines[path][4] == d_lines[path][3]),
                    ]):
                    assert d_lines[path][3] != REF
                    if d_lines[path][3] == d_lines[path][4]:
                        try:
                            assert len(d_lines[path][-1])+1 == 2*d_lines[path][-1].count('0')
                        except:
                            print(REF, d_lines[path], path)
                            print(line_out)
                            stop
                    line_out += ' '+' '.join((str(1-i) for i in map(int, d_lines[path][-1].split())))
                else:
                    print(REF, A1, A2, ID, d_lines[path][:-1], path, len(d_lines[path][-1])+1, 2*d_lines[path][-1].count('0'))
                    stop
            fd_out.write(line_out+'\n')

    with open('.'.join(args.out.split('.')[:-2])+'.sample', 'w') as fout:
        for i, path in enumerate(l_paths):
            with open('.'.join(path.split('.')[:-2])+'.sample') as f:
                lines = f.readlines()
                if i > 0:
                    lines = lines[2:]
                fout.writelines(lines)

    return


def read_fai(path_fai):

    assert os.path.isfile(path_fai)

    d = {}
    with open(path_fai) as file_fai:
        for line_fai in file_fai:
            l = line_fai.rstrip().split()
            chrom = l[0]
            byte_length = int(l[1])
            byte_start = int(l[2])
            bytes_per_line_excl_line_break = int(l[3])
            bytes_per_line_incl_line_break = int(l[4])
            d[chrom] = {'length': byte_length, 'start': byte_start}

    d['23'] = d['X']
    d['24'] = d['Y']
    d['25'] = d['X']
    d['26'] = d['MT']

    return d


def parse_ref(d_fai, fd_ref, CHROM, POS, size=1):

    ## parse ref
    cnt = cnt_chars_per_line_excluding_newline = 60
    row1 = (POS - 1) // cnt
    row2 = (POS - 1 + size) // cnt
    size += row2 - row1
    col = (POS - 1) % cnt
    byte_init = d_fai[CHROM]['start']
    offset = byte_init + (cnt + 1) * row1 + col
    fd_ref.seek(offset)
    read = fd_ref.read(size).replace('\n', '')

    return read


def argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument('--i', '--in', nargs='+')
    parser.add_argument('--out', required=True)
    parser.add_argument('--ref', required=True)

    args = namespace_args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()
