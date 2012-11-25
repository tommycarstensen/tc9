import os

path = os.getcwd()
l_fn = os.listdir(path)
for fn in l_fn:
    if os.path.isdir(fn):
        continue
    extension = fn[fn.rindex('.')+1:]
    prefix = fn[:fn.index('.')]
    ## file is not an encrypted file
    if extension != 'gpg':
        continue
    ## file was already decrypted
    if os.path.isfile('%s.tar.gz' %(prefix)):
        continue
    ## file was already decrypted
    if os.path.isdir('%s' %(prefix)):
        continue
    print prefix
##    if not os.path.isfile('%s.gpgkey' %(prefix)):
##        continue ## tmp
    fd = open('%s.gpgkey' %(prefix))
    s = fd.read()
    fd.close()
    passphrase = s.strip().split()[1]
    os.system('echo %s | gpg --passphrase-fd 0 --output %s.tar.gz --decrypt %s.tar.gz.gpg' %(passphrase,prefix,prefix,))

