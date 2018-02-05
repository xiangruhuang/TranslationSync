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

folder = 'uniform_n1p2/graph2_final/'
print("drawing figure with data under %s" % folder)

min_TL2, median_TL2, max_TL2, min_CD, median_CD, max_CD, tmean_TL2, tmean_CD, \
        ratios = process(folder)

if sys.argv[1].startswith('graph'):
    plots = {'min CD':min_CD, 'median CD':median_CD, 'max CD':max_CD,
            'min TranSync':min_TL2, 'median TranSync':median_TL2, 'max TranSync':max_TL2}
    
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
    patch = PolygonPatch(area_1, facecolor=colors['min CD'],
            edgecolor=colors['min CD'], alpha=1.0)
    ax.add_patch(patch)
    
    #plot_coords(ax, area_2.exterior)
    patch = PolygonPatch(area_2, facecolor=colors['min TranSync'],
            edgecolor=colors['min TranSync'], alpha=1.0)
    ax.add_patch(patch)
    
    #plot_coords(ax, solution.exterior)
    patch = PolygonPatch(solution, facecolor='r', edgecolor='r')
    ax.add_patch(patch)
    
    for label in ['min CD', 'median CD', 'max CD', 'min TranSync', 
            'median TranSync', 'max TranSync']:
        ax.plot(ratios, plots[label], color=colors[label], label=label,
                linestyle=linestyles[label], linewidth=linewidths[label],
                marker=markers[label])
    
    #ax.fill_between(ratios, min_CD, max_CD, facecolor=colors['median CD'], interpolate=True)
    #ax.fill_between(ratios, min_TL2, max_TL2, facecolor=colors['median TL2'], interpolate=True)
    
    legend = ax.legend(loc=(0.55, 0.55), shadow=True, fontsize=15)
    plt.title('Graph $G_{di}$', fontsize=40)
    plt.xlabel('$p$', fontsize=25)
    plt.ylabel('$\|x^*-x^{gt}\|_{\infty}$', fontsize=25)
    plt.axis([0.01, 1.0, 0, 1])
    plt.savefig('graph2.eps')
else:
    fig, ax = plt.subplots()
    plots = {'CD' : tmean_CD, 'TranSync':tmean_TL2}
    for label in ['CD', 'TranSync']:
        ax.plot(ratios, plots[label], color=colors[label], label=label,
                linestyle=linestyles[label], linewidth=linewidths[label],
                marker=markers[label])

    legend = ax.legend(loc=(0.55, 0.55), shadow=True, fontsize=15)
    plt.title('Graph $G_{dr}$', fontsize=40)
    plt.xlabel('$p$', fontsize=25)
    plt.ylabel('$Average Running Time$', fontsize=25)
    plt.axis([0.01, 1.0, 0, 3])
    plt.savefig('time1.eps')
