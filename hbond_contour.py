#!/usr/bin/env python

import os, sys
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import cm
nooa = 64

f = open('hbonds_contour.dat', 'r')

xcoords = []
ycoords = []
zcoords = []
hbonds = []

for line in f.readlines():
    x = float(line.split()[0])
    y = float(line.split()[1])
    #z = float(line.split()[2])
    hb = float(line.split()[2])

    xcoords.append(x)
    ycoords.append(y)
    #zcoords.append(z)
    hbonds.append(hb)	

plt.figure()

xlist = xcoords
ylist = ycoords               #np.linspace(-2., 1., 100)
#zlist = zcoords               #np.linspace(-1., 1., 100)
hblist = hbonds

xi = np.linspace(0., 12, 100)
yi = np.linspace(0., 12, 100)

hbi = griddata (xlist, ylist, hblist, xi, yi)


CS = plt.contour(xi,yi,hbi,15,linewidths=0.5,colors='k')
CS = plt.contourf(xi,yi,hbi,15,cmap=plt.cm.rainbow)
plt.colorbar()


plt.title('Number of hydrogen bonds per molecule wrt x and y axis')
plt.xlabel('x-axis [Angstroms]')
plt.ylabel('y-axis [Angstroms]')

#plt.clabel(CS, inline=1, fontsize=10)
plt.show()

