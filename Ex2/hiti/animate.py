#include parent folder
import os, sys, inspect

currentdir = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from src.domain import Domain

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl
import numpy as np

app = pg.mkQApp("GLScatterPlotItem Example")
w = gl.GLViewWidget()
w.opts['distance'] = 20
w.show()
w.setWindowTitle('pyqtgraph example: GLScatterPlotItem')

g = gl.GLGridItem()
w.addItem(g)


##
##  First example is a set of points with pxMode=False
##  These demonstrate the ability to have points with real size down to a very small scale 
## 
pos = np.empty((53, 3))
size = np.empty((53))
color = np.empty((53, 4))
pos[0] = (1,0,0); size[0] = 0.5;   color[0] = (1.0, 0.0, 0.0, 0.5)
pos[1] = (0,1,0); size[1] = 0.2;   color[1] = (0.0, 0.0, 1.0, 0.5)
pos[2] = (0,0,1); size[2] = 2./3.; color[2] = (0.0, 1.0, 0.0, 0.5)

phase = 0.

#w.addItem(sp2)


##
##  Third example shows a grid of points with rapidly updating position
##  and pxMode = False
##
domain = Domain()
domain.fill(10, 20, 1)
#pos3 = np.zeros((100,100,3))
#pos3[:,:,:2] = np.mgrid[:100, :100].transpose(1,2,0) * [-0.1,0.1]
#pos3 = pos3.reshape(10000,3)
pos3 = domain.pos
d3 = (pos3**2).sum(axis=1)**0.5
sp3 = gl.GLScatterPlotItem(pos=pos3, color=(1,1,1,.3), size=0.1, pxMode=False)
w.addItem(sp3)



def update():
    ## update volume colors
    global phase

    phase -= 0.1
    
    ## update surface positions and colors
    global sp3
    
    domain.verlet_advance(0.1)
    pos3 = np.array(domain.pos)
    sp3.setData(pos=pos3)
    
t = QtCore.QTimer()
t.timeout.connect(update)
t.start(50)

if __name__ == '__main__':
    pg.mkQApp().exec_()
