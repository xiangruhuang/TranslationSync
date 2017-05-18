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


folder = 'uniform_n1p2/graph'+sys.argv[1]+''
print("drawing figure with data under %s" % folder)

min_TL2, median_TL2, max_TL2, min_CD, median_CD, max_CD, tmean_TL2, tmean_CD, \
        zp_TL2, zp_CD, ratios = process(folder)

if sys.argv[2].startswith('graph'):
    plots = {'min CD':min_CD, 'median CD':median_CD, 'max CD':max_CD,
            'min TranSync':min_TL2, 'median TranSync':median_TL2, 
            'max TranSync':max_TL2}
    
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
    patch = PolygonPatch(area_1, facecolor=colors['min CD'],
            edgecolor=colors['min CD'], alpha=1.0)
    ax.add_patch(patch)
    
    patch = PolygonPatch(area_2, facecolor=colors['min TranSync'],
            edgecolor=colors['min TranSync'], alpha=1.0)
    ax.add_patch(patch)
    
    patch = PolygonPatch(solution, facecolor='r', edgecolor='r')
    ax.add_patch(patch)
    
    for label in ['min CD', 'median CD', 'max CD', 'min TranSync', 
            'median TranSync', 'max TranSync']:
        ax.plot(ratios, plots[label], color=colors[label], label=label,
                linestyle=linestyles[label], linewidth=linewidths[label],
                marker=markers[label])
    
    gs = setting_map[folder.split('/')[-1]]
    legend = ax.legend(loc=gs['loc_g'], shadow=True, fontsize=gs['fs_g'])
    plt.title(gs['title'], fontsize=40)
    plt.xlabel('$p$', fontsize=25)
    plt.ylabel('$\|x^*-x^{gt}\|_{\infty}$', fontsize=25)
    plt.axis([0.01, 1.0, 0, gs['max_diff']])
    plt.savefig('graph'+str(gs['id'])+'.eps')
elif sys.argv[2].startswith('time'):
    fig, ax = plt.subplots()
    plots = {'CD' : tmean_CD, 'TranSync':tmean_TL2}
    for label in ['CD', 'TranSync']:
        ax.plot(ratios, plots[label], color=colors[label], label=label,
                linestyle=linestyles[label], linewidth=linewidths[label],
                marker=markers[label])
    
    gs = setting_map[folder.split('/')[-1]]
    legend = ax.legend(loc=gs['loc_t'], shadow=True, fontsize=gs['fs_t'])
    plt.title(gs['title'], fontsize=40)
    plt.xlabel('$p$', fontsize=25)
    plt.ylabel('Average Running Time', fontsize=25)
    plt.axis([0.01, 1.0, 0, gs['max_time']])
    plt.savefig('time'+str(gs['id'])+'.eps')
else:
    fig, ax = plt.subplots()
    plots = {'CD' : zp_CD, 'TranSync':zp_TL2}
    for label in ['CD', 'TranSync']:
        ax.plot(ratios, plots[label], color=colors[label], label=label,
                linestyle=linestyles[label], linewidth=linewidths[label],
                marker=markers[label])
    
    gs = setting_map[folder.split('/')[-1]]
    legend = ax.legend(loc=gs['loc_t'], shadow=True, fontsize=gs['fs_t'])
    plt.title(gs['title'], fontsize=40)
    plt.xlabel('$p$', fontsize=25)
    plt.ylabel('Success Rate of Exact Recovery', fontsize=15)
    plt.axis([0.01, 1.0, 0, 1.0])
    plt.savefig('zp'+str(gs['id'])+'.eps')

