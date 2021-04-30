#include parent folder
import os, sys, inspect

currentdir = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from src.domain import *

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl
import numpy as np

<<<<<<< HEAD
N = 50
L = 1
=======
N = 75
L = 10
>>>>>>> ca3cdee092a6bf5fe2afcc698686db8f8eebb09b
speed = 0.02

app = pg.mkQApp("GLScatterPlotItem Example")
w = gl.GLViewWidget()
w.opts['distance'] = L*1.
w.show()
w.setWindowTitle('pyqtgraph example: GLScatterPlotItem')

g = gl.GLGridItem()
w.addItem(g)


##
##  First example is a set of points with pxMode=False
##  These demonstrate the ability to have points with real size down to a very small scale 
## 

domain = Domain(Epot)
domain.fill(N, L, 1)
domain.minimizeEnergy()

pos3 = domain.pos
colors = np.random.rand(N,4)
for c in colors: c[3] = 1
sp3 = gl.GLScatterPlotItem(pos=pos3, color=colors, size=0.3, pxMode=False)
w.addItem(sp3)


 

def update():    
    ## update surface positions and colors
    global sp3
    
    domain.verlet_advance(speed)
    pos3 = np.array(domain.pos)
    #print(domain.Epot(domain.pos))
    sp3.setData(pos=pos3)
    w.orbit(0.2,0)

    
t = QtCore.QTimer()
t.timeout.connect(update)
t.start(15)

if __name__ == '__main__':
    np.random.seed(1)
    pg.mkQApp().exec_()
