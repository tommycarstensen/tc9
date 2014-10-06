import sys
import os


def main():

    path_in = sys.argv[-1]
    path_out = os.path.basename(path_in)+'.transposed'
    separator = ' '

    d_seek = {}
    with open(path_in) as fd_in:
        i = 0
        print('indexing')
        while True:
            tell = fd_in.tell()
            if fd_in.readline() == '':
                break
            d_seek[i] = tell
            i += 1
        print('indexed')
    cols2 = rows1 = i

    with open(path_in) as fd_in:
        line = fd_in.readline()
    rows2 = cols1 = len(line.split(separator))
    del line

    with open(path_in) as fd_in, open(path_out, 'w') as fd_out:
        print('transposing')
        for row2 in range(rows2):
            print('row', row2)
            for row1 in range(rows1):
                fd_in.seek(d_seek[row1])
                s = ''
                while True:
                    char = fd_in.read(1)
                    if char == separator or char == '\n':
                        break
                    s += char
                d_seek[row1] += len(s)+1
                if row1+1 < rows1:
                    fd_out.write('{} '.format(s))
                else:
                    fd_out.write('{}\n'.format(s))

    return


if __name__ == '__main__':
    main()
