#!/bin/python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, June 2013

from ftplib import FTP


def main():

    parse_ftp()

    return


def parse_ftp():

    path_ftp = 'pub/release-71/emf/ensembl-compara/epo_6_primate'    
    ftp = FTP('ftp.ensembl.org')
    ftp.login()
    print(ftp.dir(path))
    print(ftp.nlst(path))
    print(ftp.mlsd(path))
    stop
    for f in l_files:
        if f[-3:] != '.bb':
            continue
        path_out = os.path.join(k,os.path.basename(f))
        if os.path.isfile(path_out):
            continue
        ftp.retrbinary(
            'RETR %s' %(f),
            open(path_out,'wb').write,
            )
    ftp.quit()
    for dn in d_paths.keys():
        l_files = os.listdir(dn)
        for f in l_files:
            if f[-3:] != '.bb':
                continue
            path_in = os.path.join(dn,f)
            path_out = os.path.join(dn,f[:-3]+'.bed')
            if os.path.isfile(path_out):
                continue
            ## assume bigBedToBed binary is in cwd
            cmd = './bigBedToBed %s %s' %(path_in,path_out)
            os.system(cmd)

    return


if __name__ == '__main__':
    main()
