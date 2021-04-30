#include parent folder
import os, sys, inspect

currentdir = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)


import numpy as np
import jax
import scipy

import src.misc as misc
import src.domain as domain

def Epot(positions):
    """Lennard-Jones potential in reduced units.
    In this system of units, epsilon=1 and sigma=2**(-1. / 6.).
    """
    positions = misc.to3D(positions)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("positions must be an Mx3 array")
    # Compute all squared distances between pairs without iterating.
    delta = positions[:, np.newaxis, :] - positions
    r2 = (delta * delta).sum(axis=2)
    # Take only the upper triangle (combinations of two atoms).
    indices = np.triu_indices(r2.shape[0], k=1)
    rm2 = 1. / r2[indices]
    # Compute the potental energy recycling as many calculations as possible.
    rm6 = rm2 * rm2 * rm2
    rm12 = rm6 * rm6
    return (rm12 - 2. * rm6).sum()

Epot_jit = jax.jit(Epot)
grad_Epot_jit = jax.jit(jax.grad(Epot))

#pos = np.random.rand(10,3)
d = domain.Domain()
d.fill(10,1,1)
pos = d.pos
print(Epot(pos))

result = scipy.optimize.minimize(
    Epot, 
    pos.ravel(),
    jac=grad_Epot_jit,
    method='CG',
    options={"disp" : True, "maxiter" : 10}           
  )

print(result)