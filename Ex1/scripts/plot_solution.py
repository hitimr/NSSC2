import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from itertools import chain
import os
import subprocess
import sys

from common import *


fileName = sys.argv[1]
data = np.genfromtxt(WORKSPACE_DIR + fileName)

x = np.linspace(0, 1, len(data[0]))
y = np.linspace(0, 1, len(data))

X, Y = np.meshgrid(x, y)


Z = list(chain.from_iterable(data))
Z = np.sinc(np.sqrt(X ** 2 + Y ** 2))

fig = plt.figure()
ax = fig.gca(projection = '3d')
ax.plot_surface(X, Y, data)
plt.title("parallel")
plt.show()