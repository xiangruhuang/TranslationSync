import sys
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os.path
from settings import *
import shapely
from shapely.geometry import Polygon
from descartes import PolygonPatch
from utils import *

fig, ax = plt.subplots()

def get_dist(filename):
    fin = open(filename, 'r')

    a = [float(line.strip()) for line in fin.readlines()]

    b = [a[0]]
    c = [1]
    density = 0

    for i in range(1, len(a)):
        if a[i] == b[-1]:
            c[-1] += 1
        else:
            b.append(a[i])
            c.append(c[-1]+1)

    c = [num*1.0/len(a) for num in c]
    return b, c

CD_x, CD_y = get_dist('CD.dist')
TL2_x, TL2_y = get_dist('TL2.dist')

ax.plot(CD_x, CD_y, color='k', label='CD')
ax.plot(TL2_x, TL2_y, color='b', label='TranSync')

legend = ax.legend(loc=(0.13, 0.03), shadow=True, fontsize=40)
plt.title('Cumulative Density Function', fontsize=30)
plt.xlabel('$\|t_{ij}-(x_i-x_j)\|$', fontsize=20)
plt.ylabel('Cumulative Density', fontsize=25)
plt.axis([0.0, 200, 0, 1])
plt.savefig('cdf.eps')

