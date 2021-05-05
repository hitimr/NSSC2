#include parent folder
import os, sys, inspect

currentdir = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

# replace with you module
from src.domain import *

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl
import numpy as np


# settings
rho = 0.8
N = 200
L = (float(N/rho)**(1/3))
speed = 0.01
particle_size = 0.2
orbit = 0.0
sigma = 1


# Window
app = pg.mkQApp("GLScatterPlotItem Example")
w = gl.GLViewWidget()
w.opts['distance'] = L*3
w.show()
w.setWindowTitle('pyqtgraph example: GLScatterPlotItem')

# Draw a box
# from https://stackoverflow.com/questions/53136301/plot-cube-using-pyqtgraph-in-python
vertexes = np.array([[1, 0, 0], #0
                     [0, 0, 0], #1
                     [0, 1, 0], #2
                     [0, 0, 1], #3
                     [1, 1, 0], #4
                     [1, 1, 1], #5
                     [0, 1, 1], #6
                     [1, 0, 1]])#7

vertexes = vertexes*L

faces = np.array([[1,0,7], [1,3,7],
                  [1,2,4], [1,0,4],
                  [1,2,6], [1,3,6],
                  [0,4,5], [0,7,5],
                  [2,4,5], [2,6,5],
                  [3,6,5], [3,7,5]])

colors = np.array([[1,0,0,1] for i in range(12)])
cube = gl.GLMeshItem(vertexes=vertexes, faces=faces, faceColors=colors,
                     drawEdges=True, edgeColor=(1,1,1,0.5), drawFaces=False)

cube.translate(-L/2, -L/2, -L/2)
w.addItem(cube)


# Replace with your initialized Domain
domain = Domain(Epot)
domain.fill(N, L, 0.01, sigma)
domain.minimizeEnergy()



pos = domain.pos # -> Replace with your initial position


# draw colored spheres
colors = np.random.rand(N,4)
for c in colors: c[3] = 1
sp = gl.GLScatterPlotItem(pos=pos, color=colors, size=particle_size, pxMode=False)
w.addItem(sp)


 
## update surface positions and colors
def update():  
    global sp
    
    # reeplace with you code
    domain.verlet_advance(speed)    # advance domain in time
    new_pos = np.array(domain.pos)  # get new position

    origin = new_pos[0]
    new_pos = new_pos - origin

    new_pos = new_pos - domain.length * np.round(new_pos/domain.length) # minimum image convention
    sp.setData(pos=new_pos)
    w.orbit(orbit,0)    # Rotate camera

    
t = QtCore.QTimer()
t.timeout.connect(update)
t.start(15) # update intervall in ms

if __name__ == '__main__':
    np.random.seed(1)
    pg.mkQApp().exec_()