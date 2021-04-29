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

domain = Domain(Epot)
domain.fill(40, 10, 1)
domain.minimizeEnergy()

pos3 = domain.pos
sp3 = gl.GLScatterPlotItem(pos=pos3, color=(1,1,1,10), size=1, pxMode=False)
w.addItem(sp3)

 

def update():    
    ## update surface positions and colors
    global sp3
    
    domain.verlet_advance(0.1)
    pos3 = np.array(domain.pos)
    #print(domain.Epot(domain.pos))
    sp3.setData(pos=pos3)
    
t = QtCore.QTimer()
t.timeout.connect(update)
t.start(50)

if __name__ == '__main__':
    np.random.seed(1)
    pg.mkQApp().exec_()
