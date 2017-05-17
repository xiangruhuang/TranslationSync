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

def plot_coords(ax, ob):
    x, y = ob.xy
    ax.plot(x, y, '+', color='grey')

folder = 'uniform_n1p2/graph1/'
print("drawing figure with data under %s" % folder)

min_TL2 = []
median_TL2 = []
max_TL2 = []

min_CD = []
median_CD = []
max_CD = []

ratios = []
for ratio in range(100):
    ml_TL2 = []
    time_TL2 = []
    ml_CD = []
    time_CD = []
    for eid in range(100):
        name_TL2 = folder+'/ratio'+str(ratio)+'_'+str(eid)+'.TL2'
        name_CD = folder+'/ratio'+str(ratio)+'_'+str(eid)+'.CD'
        if os.path.isfile(name_TL2) and os.path.isfile(name_CD):
            if len(ratios) == 0 or ratios[-1] < ratio/100.0-0.005:
                ratios.append(ratio/100.0)
            with open(name_TL2) as fin:
                lines = fin.readlines()
                min_loss=1e100
                time=0.0
                for line in lines:
                    ml = 0.0
                    t = 0.0
                    for token in line.strip().split(', '):
                        if token.startswith('min_loss'):
                            ml = float(token.split('=')[1])
                        if token.startswith('elapsed'):
                            t = float(token.split('=')[1])
                    if min_loss > ml:
                        min_loss = ml
                        time = t
                    if min_loss < 1e-5:
                        break
                ml_TL2.append(min_loss)
                time_TL2.append(time)

            with open(name_CD) as fin:
                lines = fin.readlines()
                min_loss=1e100
                time=0.0
                for line in lines:
                    ml = 0.0
                    t = 0.0
                    for token in line.strip().split(', '):
                        if token.startswith('min_loss'):
                            ml = float(token.split('=')[1])
                        if token.startswith('elapsed'):
                            t = float(token.split('=')[1])
                    if min_loss > ml:
                        min_loss = ml
                        time = t
                    if min_loss < 1e-5:
                        break
                ml_CD.append(min_loss)
                time_CD.append(time)
    if len(ml_TL2) > 0 and len(ml_CD) > 0:
        min_TL2.append(min(ml_TL2))
        median_TL2.append(numpy.median(ml_TL2))
        max_TL2.append(max(ml_TL2))
        min_CD.append(min(ml_CD))
        median_CD.append(numpy.median(ml_CD))
        max_CD.append(max(ml_CD)) 

ratios = [1.0-ratio for ratio in ratios]

plots = {'min_CD':min_CD, 'median_CD':median_CD, 'max_CD':max_CD,
        'min_TL2':min_TL2, 'median_TL2':median_TL2, 'max_TL2':max_TL2}

down_CD = [(x, y) for (x, y) in zip(ratios, min_CD)]
up_CD = [(x, y) for (x, y) in zip(ratios, max_CD)]
up_CD.reverse()

down_TL2 = [(x, y) for (x, y) in zip(ratios, min_TL2)]
up_TL2 = [(x, y) for (x, y) in zip(ratios, max_TL2)]
up_TL2.reverse()

area_1 = Polygon(down_CD+up_CD)

area_2 = Polygon(down_TL2+up_TL2)
    
solution = area_1.intersection(area_2)

fig, ax = plt.subplots()
#plot_coords(ax, area_1.exterior)
patch = PolygonPatch(area_1, facecolor=colors['min_CD'],
        edgecolor=colors['min_CD'], alpha=0.1)
ax.add_patch(patch)

#plot_coords(ax, area_2.exterior)
patch = PolygonPatch(area_2, facecolor=colors['min_TL2'], edgecolor=colors['min_TL2'], alpha=0.1)
ax.add_patch(patch)

#plot_coords(ax, solution.exterior)
patch = PolygonPatch(solution, facecolor='r', edgecolor='r')
ax.add_patch(patch)

for label in plots.keys():
    ax.plot(ratios, plots[label], color=colors[label], label=label,
            linestyle=linestyles[label], linewidth=linewidths[label],
            marker=markers[label])

#ax.fill_between(ratios, min_CD, max_CD, facecolor=colors['median_CD'], interpolate=True)
#ax.fill_between(ratios, min_TL2, max_TL2, facecolor=colors['median_TL2'], interpolate=True)

legend = ax.legend(loc=(0.65, 0.55), shadow=True, fontsize=15)
plt.title('Graph1', fontsize=40)
plt.xlabel('$p$', fontsize=25)
plt.ylabel('$\|x-x^{gt}\|_{\infty}$', fontsize=25)
plt.axis([0, 1.0, -0.1, 1])
plt.savefig('graph1.eps')

