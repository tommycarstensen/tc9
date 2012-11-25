import os, time

def main():

    path = os.getcwd()
    l_fn = os.listdir(path)
    for fn in l_fn:
        if os.path.isdir(fn):
            continue
        extension = fn[fn.rindex('.')+1:]
        prefix = fn[:fn.index('.')]
        if prefix == '305_PD4777b':
            continue ## tmp
        ## file is not an encrypted file
        if extension != 'gz':
            continue
        ## file was already gunzipped and untarred
        if os.path.isfile('%s.tar' %(prefix)) and os.path.isdir(prefix):
            remove(prefix)
            continue
        ## file is not currently being decrypted?
        seconds = time.time()-os.path.getmtime('%s.tar.gz' %(prefix))
        if seconds < 5*60:
            continue
        print 'gunzip', prefix
        os.system('gunzip %s.tar.gz' %(prefix,))
        print 'tar', prefix
        os.system('tar -xvf %s.tar' %(prefix,))
        remove(prefix)

    return


def remove(prefix):

    if os.path.isfile('%s.gpgkey' %(prefix)):
        os.remove('%s.gpgkey' %(prefix))
    if os.path.isfile('%s.tar.gz.gpg' %(prefix)):
        os.remove('%s.tar.gz.gpg' %(prefix))

    return

if __name__ == '__main__':
    main()
