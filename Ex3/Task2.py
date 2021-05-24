import numpy as np
import matplotlib.pyplot as plt

from animate import animate

def func_gauss(x):
    return np.exp(-10.*np.power(4*x-1,2))

def func_square(x):
    if hasattr(x, "__len__"): return [square(xi) for xi in x]
    if(0.1 < x) and (x < 0.3): return 1.
    else: return 0.

def advance_upwind(C, Co):
    new_C = np.copy(C)
    for i in range(1, len(C)):
        new_C[i] = (1 - Co) * C[i] + Co * C[i-1]
    return new_C

def enforce_pb(x, L=1.):
    """Enforce Periodic Boundary Conditions using minimum image convention

    Args:
        x (double): coordinate
        L (double, optional): size of domain. Defaults to 1.

    Returns:
        (double): mapped coodrdinate
    """    
    return x - L * np.round(x/L) 

def advect_analytical(initFunc, Nx, steps, c=1.):
    x_vals = np.linspace(0,1, Nx)
    return [initFunc(enforce_pb(x_vals - c*t)) for t in np.linspace(0,1,steps)]


def advect_upwind(initFunc, Co, steps):
    x_vals = np.linspace(0,1, Nx)
    C = initFunc(x_vals)

    res = []
    for i in range(steps):
        C = advance_upwind(C, Co)
        res.append(C)

    res = np.array(res)
    return res


if __name__ == "__main__":
    dx = 0.01
    Nx = int(1/dx)
    x_vals = np.linspace(0, 1, Nx)

    #res = advect_analytical(func_gauss, Nx=200, steps=15, c=1)
    animate(res, "out/task2_gauss_co_1", title="Courant=1")