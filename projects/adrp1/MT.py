import glob
import os


def main():

    with open('/lustre/scratch115/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa') as f:
        f.seek(3147273397)
        refseq = f.read(16805).replace('\n', '')

    for path in glob.glob('phylotree/*.txt'):
        haplogroup = os.path.splitext(os.path.basename(path))[0]
        with open(path) as f:
            for line in map(str.strip, f):
                ## Skip blank lines.
                if not line:
                    continue
                for char in '!()':
                    line = line.replace(char, '')
                if line[-1] == 'd':
                    pos = int(line[1:-1])
                    alt = '-'
                    ref = refseq[pos-1]
                elif line[-3] == '.' and line[-2] == '1' and line[-1] in 'ACGT':
                    pos = int(line[:-3])
                    ref = refseq[pos-1]
                    alt = ref+line[-1]
                else:
                    pos = int(line[1:-1])
                    ref = refseq[pos-1]
                    if ref == line[0].upper():
                        alt = line[-1].upper()
                    else:
                        alt = line[0].upper()
                        assert ref == line[-1].upper()
                print(line, ref, pos, alt, sep='\t')

    with open('MT-RNR1.txt') as f:
        for line in map(str.strip, f):
            ID, pos, ref, alt = line.split()
            print(ref+pos+alt, ref, pos, alt, sep='\t')

    return


def conversion(pos):

    ## https://www.snpedia.com/index.php/MtDNA_Position_Conversions
    if pos >= 1 and pos <= 309:
        return pos

    return pos


if __name__ == '__main__':
    main()
