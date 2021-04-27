import numpy as np


def to3D(v):
    """Convert a 1D data to matrix (array of vector3)

    Args:
        v (array-like): vector to transform
    """    
    v = np.array(v) 
    if (len(v.shape) == 2) and (v.shape[1] == 3): return v  # data is already 3D

    assert(len(v.shape) == 1)
    assert(len(v) % 3 == 0)
    
    n = len(v) // 3
    return np.reshape(v, (n,3))
