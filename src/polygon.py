import matplotlib.pyplot as plt
import shapely
from shapely.geometry import Polygon
from descartes import PolygonPatch
import numpy as np

def create_tube(a,height):
    x_tube_up = np.linspace(-4,4,300)
    y_tube_up = a*x_tube_up**2 + height
    x_tube_down = np.flipud(x_tube_up)          #flip for correct definition of polygon
    y_tube_down = np.flipud(y_tube_up - 2)

    points_x = list(x_tube_up) + list(x_tube_down)
    points_y = list(y_tube_up) + list(y_tube_down)

    return Polygon([(points_x[i], points_y[i]) for i in range(600)])

def plot_coords(ax, ob):
    x, y = ob.xy
    ax.plot(x, y, '+', color='grey')


area_1 = Polygon()          #First area, a MultiPolygon object
for h in [-5, 0, 5]:
    area_1 = area_1.union(create_tube(2, h))

area_2 = Polygon()
for h in [8, 13, 18]:
    area_2 = area_2.union(create_tube(-1, h))

solution = area_1.intersection(area_2)      #What I was looking for

##########  PLOT  ##########

fig = plt.figure()
ax = fig.add_subplot(111)

for tube in area_1:
    plot_coords(ax, tube.exterior)
    patch = PolygonPatch(tube, facecolor='g', edgecolor='g', alpha=0.1)
    ax.add_patch(patch)

for tube in area_2:
    plot_coords(ax, tube.exterior)
    patch = PolygonPatch(tube, facecolor='m', edgecolor='m', alpha=0.1)
    ax.add_patch(patch)

for sol in solution:
    plot_coords(ax, sol.exterior)
    patch = PolygonPatch(sol, facecolor='r', edgecolor='r')
    ax.add_patch(patch)

plt.show()
