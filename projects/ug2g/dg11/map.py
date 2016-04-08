import shapefile
import urllib.request
import os
import zipfile
import shapefile
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import cm

## Download the shp file:
url = 'http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/cultural/ne_50m_admin_0_countries.zip'

if not os.path.isfile(os.path.basename(url)):

    response = urllib.request.urlretrieve(url, os.path.basename(url))
    zip_ref = zipfile.ZipFile(os.path.basename(url), 'r')
    zip_ref.extractall('.')
    zip_ref.close()

ax = plt.subplot(111, frameon=False)

## Instanciate the Basemap
#m = Basemap(resolution='i',projection='merc', llcrnrlat=y1,urcrnrlat=y2,llcrnrlon=x1,urcrnrlon=x2,lat_ts=(x1+x2)/2)
m = Basemap(
    resolution='i',
    projection='merc',
    llcrnrlat=-37.5, urcrnrlat=40, llcrnrlon=-20, urcrnrlon=55,
    )
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
#m.drawparallels(np.arange(y1,y2,2.),labels=[1,0,0,0],color='black',dashes=[1,0],labelstyle='+/-',linewidth=0.2) # draw parallels
#m.drawmeridians(np.arange(x1,x2,2.),labels=[0,0,0,1],color='black',dashes=[1,0],labelstyle='+/-',linewidth=0.2) # draw meridians
m.fillcontinents(color='#a6bddb', lake_color='white')
 
## Then import and get the shapes (from the shp file) and records (from the dbf files):
r = shapefile.Reader(os.path.basename(url))
shapes = r.shapes()
records = r.records()
for index_sovereignt, field in enumerate(r.fields[1:]):
    if field[0] == 'sovereignt':
        break
for index_continent, field in enumerate(r.fields[1:]):
    if field[0] == 'continent':
        break

## then, we get the points we need to create the matplotlib polygons:
for record, shape in zip(records, shapes):
    if record[index_continent] != 'Africa':
        continue
    lons, lats = zip(*shape.points)
    data = np.array(m(lons, lats)).T
    sovereignt = record[index_sovereignt]
    print(sovereignt)


## within this loop, we will also check that the shape is not in multiple parts, and if yes, segment the points in different ensembles:

    if sovereignt == 'Uganda':
        color = '#2b8cbe'
    else:
        color = '#a6bddb'
        continue

    if len(shape.parts) == 1:
        segs = [data,]
    else:
        segs = []
        for i in range(1,len(shape.parts)):
            index = shape.parts[i-1]
            index2 = shape.parts[i]
            segs.append(data[index:index2])
        segs.append(data[index2:])
 
    lines = LineCollection(segs,antialiaseds=(1,))
    lines.set_facecolors(color)
    lines.set_edgecolors('k')
    lines.set_linewidth(0.1)
    ax.add_collection(lines)

m.drawcoastlines(linewidth=0.5)


plt.savefig('map.png', dpi=600)
plt.close()

