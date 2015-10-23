#!/bin/env python3

import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import argparse
import matplotlib.patches as mpatches

colors = [
    '#a6cee3', '#1f78b4', '#b2df8a', '#33a02c',
    'red','blue','green','yellow','magenta','purple',
    ]

 
def draw_pie(ax, ratios, X, Y, size):

    xy = []
 
    start = 0.
    for ratio in ratios:
        c = np.linspace(2*math.pi*start,2*math.pi*(start+ratio), 30)
        x = [0] + np.cos(c).tolist()
        y = [0] + np.sin(c).tolist()
        xy.append(list(zip(x,y)))
        start += ratio
 
    for i, xyi in enumerate(xy):
        ax.scatter(
            [X],[Y] , marker=(xyi,0), s=2*size, facecolor=colors[i])


def main():

    args = parseargs()

    fig = plt.figure(figsize=(11.7,8.3))
    #Custom adjust of the subplots
    plt.subplots_adjust(
        left=0.05,right=0.95,top=0.90,bottom=0.05,wspace=0.15,hspace=0.05)
    ax = plt.subplot(111)
    #Let's create a basemap around Belgium
    m = Basemap(
        resolution='i',projection='merc',
        llcrnrlat=args.corners[0], urcrnrlat=args.corners[1],
        llcrnrlon=args.corners[2], urcrnrlon=args.corners[3],
        lat_ts=51.0)
    m.drawcountries(linewidth=0.5)
    m.drawcoastlines(linewidth=0.5)
     
##    m.drawparallels(
##        np.arange(args.corners[0],args.corners[2],1.),labels=[1,0,0,0],color='black',dashes=[1,0],
##        labelstyle='+/-',linewidth=0.2) # draw parallels
##    m.drawmeridians(
##        np.arange(args.corners[1],args.corners[3],1.),labels=[0,0,0,1],color='black',dashes=[1,0],
##        labelstyle='+/-',linewidth=0.2) # draw meridians

    with open(args.input) as f:
        for i, line in enumerate(f):
##    for i, (lat, lon, size, ratios) in enumerate(lats, lons, sizes, ratios):
            if line[0] == '#':
                continue
            if line == '\n':
                continue
            l = line.rstrip().split()
            text = l[0]
            lat = float(l[2])
            lon = float(l[1])
            lat += float(l[4])
            lon += float(l[3])
            size = float(l[5])
            ratios = list(map(float, l[6:]))
            assert len(ratios) == len(args.labels)
            X,Y = m(lat, lon)
            draw_pie(ax, ratios, X, Y, size)
            ## Do labels for legend.
            if i == 0:
                handles = []
                for j, label in enumerate(args.labels):
                    print(label, colors[j])
                    patch = mpatches.Patch(
                        color=colors[j], label=label)
                    handles.append(patch)
                    ax.scatter(
                        200, 200, marker='o', s=48, c=colors[j], label=label)
            ## Do text label.
            ax.annotate(
                text, xy=(X, Y),
                xytext=(0, 30), textcoords='offset points',
                horizontalalignment='center', verticalalignment='top')
                    
    plt.legend(numpoints=1, scatterpoints=1)

    plt.savefig(args.output+'.png')

    plt.show()


def parseargs():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input',
        help='Input file with columns label, lon, lat, offsetlon, offsetlat, circle size, slice sizes (%)',
        required=True,
        )
    parser.add_argument(
        '--labels', help='Pie chart titles', nargs='+')
    parser.add_argument(
        '--separator', help='File separator', choices=[',','\t'],
        default='\t')
    parser.add_argument(
        '--corners', help='llcrnrlat urcrnrlat llcrnrlon urcrnrlon',
        type=float, nargs=4)
    parser.add_argument('--output', required=True)

    args = parser.parse_args()

    assert args.corners[0] < args.corners[1]
    assert args.corners[2] < args.corners[3]

    return args


if __name__ == '__main__':
    main()
