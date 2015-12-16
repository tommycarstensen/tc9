import matplotlib as mpl
## http://matplotlib.org/faq/howto_faq.html#matplotlib-in-a-web-application-server
mpl.use('Agg')
import matplotlib.pyplot as plt

with open('tmp','r') as f:

    l_OR = []
    for line in f:

        l = line.rstrip().split()

        if l[0] == 'inf': continue

        try:
            OR = float(l[0])
            if OR > 200: continue
            l_OR += [OR]
        except: continue

    print(len(l_OR))
    plt.hist(l_OR,200)
    plt.savefig('hist_OR.png')

