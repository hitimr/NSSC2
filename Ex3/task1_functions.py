import numpy as np
import matplotlib.pyplot as plt
import argparse
from misc import *
import animate



def analytical_sol(t, points=100, n_max=100):
    x_vals = np.linspace(1,0,points)
    n = np.arange(0, n_max+1, 1)

    # create vector with alternating signs so we dont need (-1)^n
    signs = np.empty((n_max+1,),int)
    signs[::2] = 1
    signs[1::2] = -1

    # pull out redundant calculations
    t = 10**(-6)*t
    pi2t = np.pi*np.pi*t
    

    # solution for a single point
    def analytical_sol_single(x):    
        return 1-2*np.sum(
           signs / ((n+0.5)*np.pi) * np.cos((n+0.5)*np.pi*x) * np.exp( -(n*n+n+0.25) * pi2t)
        )
    
    # vectorize over x_vals for parallel performance
    vectorized = np.vectorize(analytical_sol_single)
    return vectorized(x_vals)



def apply_boundaries(C,task):
    if task == 1:
        C[0] = 1
        C[-1] = C[-2]
    if task == 2:
        C[0] = 1
        C[len(C)-1] = 0
    return C


def advance(C,task,D,S):
    new_C = np.empty(len(C))
    for i in range(1, len(C)-1):
        new_C[i] = C[i] + D*S*(C[i-1] - 2*C[i] + C[i+1])
    apply_boundaries(new_C,task)
    return new_C



def run(Nx,dt,task,num_iter,printmod,filename):
    dx = (1/(Nx+1))
    S = dt/(dx*dx)
    D = 10**(-6)
    d = D*dt/dx**2

    C = np.zeros(Nx)
    C = apply_boundaries(C,task)
    res = [C]

    d=D*dt/dx**2

    for i in range(int(num_iter)):
        C = advance(C,task,D,S)
        if i % printmod == 0:
            res.append(C)
    np.savetxt(DIR_OUT + filename,res)

    return res





