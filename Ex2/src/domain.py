#include parent folder
import os, sys, inspect

currentdir = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import time
import numpy as np
import scipy
from scipy.optimize import minimize

import jax
import jax.numpy as jnp
from jax import grad, jit

from config import *
from src.misc import *

def Epot(x):
    x = to3D(x)
    n = len(x)
    broadcasted = x[:, jnp.newaxis, :] - x   # calculate the distance between all particles while utilizing jnp.broadcast
    broadcasted = broadcasted.reshape(n*n,3)
    r = jnp.linalg.norm(broadcasted, axis=1)  # calculate the euclidian norm of all distances
    r6 =  jnp.nan_to_num(jnp.power(r, -6), posinf=0.0)
    return 4*(2*jnp.dot(r6,r6) - r6.sum())



# Leonnard-Jones Potential in natural system of units
def VLJ(v, w):   
    r = jnp.linalg.norm(v-w) # euclidian distance between points
    #if(r == 0): return 0
    return 4 * (4*pow(1 / r, 12) - 2*pow(1 / r, 6)) # TODO: check if this is actually correct

# Leonard Jones with distance as argument
def VLJ_r(r):   
    if(r == 0): return 0
    return 4 * (4*pow(1 / r, 12) - 2*pow(1 / r, 6)) # TODO: check if this is actually correct



# https://jax.readthedocs.io/en/latest/jax.html#jax.grad
Epot_jit = jax.jit(Epot)
grad_Epot = grad(Epot)
grad_Epot_jit = jit(grad(Epot))  

class Domain:
    # Domain settings
    particle_count = int  # Number of particles in the box
    length = float  # length of the box
    std_dev = float  # standard deviation of the velocities
    pos = np.ndarray  # positional data of the particles
    vel = np.ndarray  # velocity data of the particle

    # Stats
    time_minimizeEnergy = float

    # Physics
    Epot = object
    grad_Epot = object

    def __init__(self, fEpot=Epot_jit, fGradEpot=grad_Epot_jit):
        self.particle_count = -1
        self.length = -1
        self.std_dev = -1

        # Benchmark stats
        self.time_minimizeEnergy = -1

        self.Epot = fEpot
        self.grad_Epot = fEpot

    def fill(self, particle_count, length, std_dev):
        # TODO: Move length arguemnt to constructor
        """Fill the domain with particles as specified in Task 2

        Args:
            particle_count (int): number of particles
            length (float): length of the simulation box
            std_dev (float): standard deviation for the velocity distribution
        """
        assert (isinstance(particle_count, int))
        assert (particle_count > 0)
        assert (length > 0)
        assert (std_dev >= 0)

        self.particle_count = particle_count
        self.length = length
        self.std_dev = std_dev

        self.initialize_pos()
        self.initialize_vel()

        assert (self.pos.shape == self.vel.shape)

    def verlet_advance(self, dt):
        if(dt == 0): return

        # TODO: optimize, maybe switch to jax
        f = -self.grad_Epot(self.pos)        

        new_pos = self.pos + (self.vel + 0.5 * f * dt) * dt
        new_f = -self.grad_Epot(new_pos)
        new_vel = self.vel + 0.5 * (f + new_f ) * dt

        self.pos = new_pos
        self.vel = new_vel

    def initialize_pos(self):
        # TODO, Reno: replace with custom distribution from Task 2.1
        self.pos =  np.random.rand(self.particle_count, 3) * self.length



    def initialize_vel(self):
        # TODO, Reno: replace with custom distriburtion from Task 2.3 and 2.4
        #self.vel = np.random.rand(self.particle_count, 3) - 0.5
        self.vel = np.zeros((self.particle_count, 3))

    def minimizeEnergy(self):
        start = time.time()
        result = minimize(
            self.Epot, 
            self.pos.ravel(),
            jac=self.grad_Epot,
            method='CG',
            options={"disp" : True, "maxiter" : 10}           
        )

        new_pos = to3D(result.x)
        self.pos = new_pos

        end = time.time()
        self.time_minimizeEnergy = end - start
        return new_pos  


    def read_from_file(self, fileName):
        """Fill the domain with data from a file
        #TODO, Hickel: read n-th block 

        Args:
            fileName (str): path to file
        """
        f = open(fileName+".txt", "r")
        self.particle_count=int(f.readline())
        comment=f.readline()
        self.length=float(f.readline())
        for a in range(0,self.particle_count):
            line_tmp=f.readline().split(" ")
            for b in range(0,2):
                self.pos[a][b]=float(line_tmp[b])
                self.vel[a][b]=float(line_tmp[b+2])
        f.close()
        assert (self.particle_count > 0)
        assert (self.length > 0)
        assert (self.pos.shape == self.vel.shape)

        return

    def write_to_file(self, fileName, comment=""):
        """Write domain data to a file

        #TODO, Hickel: Append Blocks

        Args:
            fileName (str): path to file
            comment (str): arbitrary comment/description for Line 2
        """
        f = open(fileName+".txt", "w")
        f.write(str(self.particle_count)+"\n")
        f.write(comment+"\n")
        f.write(str(self.length)+"\n")
        for a in range(0,self.particle_count):
            f.write(str(self.pos[a][0])+" "+str(self.pos[a][1])+" "+str(self.pos[a][2])+" "+str(self.vel[a][0])+" "+str(self.vel[a][1])+" "+str(self.vel[a][2])+"\n")
        f.close()		

        return

    def average_distance(self):
        distances = np.triu(scipy.spatial.distance.cdist(domain.pos, domain.pos)) # calculate distances but throws away lower triangle
        average = distances.sum() / (2*self.particle_count**2)
        return average


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    N = 20
    rho = 0.8

    L = ((float(N)/rho))**(1.0/3.0)
    np.random.seed(1)
    domain = Domain(Epot_jit, grad_Epot_jit)


    domain.fill(100, 20, 1)   

    print(domain.vel)
    pos = domain.pos.T
    vel = domain.vel.T


    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.quiver(pos[0], pos[1], pos[2], vel[0], vel[1], vel[2], length=3)
    plt.show()



