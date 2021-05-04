import numpy as np
import jax.numpy as jnp
import os
import platform
import inspect

# System specific directory separator
if platform.system() in ["Darwin", "Linux"]:
    SEP = "/"
else: # platform.system() == "Windows":
    SEP = "\\"

DIR_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))) + SEP
DIR_OUT = DIR_ROOT + "out" + SEP

def to3D(v):
    """Convert a 1D data to matrix (array of vector3)

    Args:
        v (array-like): vector to transform
    """    
    v = jnp.array(v) 
    if (len(v.shape) == 2) and (v.shape[1] == 3): return v  # data is already 3D

    assert(len(v.shape) == 1)
    assert(len(v) % 3 == 0)
    
    n = len(v) // 3
    return np.reshape(v, (n,3))
