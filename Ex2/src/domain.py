#include parent folder
import os, sys, inspect

currentdir = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)


import numpy as np
import scipy
from scipy.optimize import minimize

import jax
import jax.numpy as jnp
from jax import grad, jit

from config import *
from misc import *


def Epot(pos):    
    # scipy optimize passes flattened array to E_pot. If pos is 1D transform
    # it to 3D and back again after the calculation has finished
    #if len(pos.shape) == 1:
    v = to3D(pos)

    #E = jnp.triu(scipy.spatial.distance.cdist(v, v, metric=VLJ)).sum() # LOL that actually works :D


    
    # Old way to calculate the potential. Still left in here to verify results
    # Calculate the potential energy for all pairs of molecules. Skip cases
    # where i=j and distances that have already been calculated
    E = 0
    for i in range(len(domain.pos)):
        for j in range(i + 1, len(domain.pos)):
            E += VLJ(domain.pos[i], domain.pos[j])
    
    return E   


# Leonnard-Jones Potential in natural system of units
def VLJ(v, w):   
    r = np.linalg.norm(v-w) # euclidian distance between points
    if(r == 0): return 0
    return 4 * (4*pow(1 / r, 12) - 2*pow(1 / r, 6)) # TODO: check if this is actually correct

# TODO remove unnecessry ones whern done
grad_Epot = jit(grad(Epot))  # https://jax.readthedocs.io/en/latest/jax.html#jax.grad
Epot_jit = jit(Epot)
VLJ_jit = jit(VLJ)

class Domain:
    particle_count = int  # Number of particles in the box
    length = float  # length of the box
    std_dev = float  # standard deviation of the velocities
    pos = np.ndarray  # positional data of the particles
    vel = np.ndarray  # velocity data of the particles

    def __init__(self):
        self.particle_count = -1
        self.length = -1
        self.std_dev = -1

    def fill(self, particle_count, length, std_dev):
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

    def initialize_pos(self):
        # TODO, Reno: replace with custom distribution from Task 2.1
        self.pos = np.ndarray((self.particle_count, 3))

        spacing = self.length / (self.particle_count)
        for i in range(self.particle_count):
            for j in range(3):
                self.pos[i][j] = spacing / 2 + i*spacing/10 # + 1/2 to avoid placing particles at (0|0)

    def initialize_vel(self):
        # TODO, Reno: replace with custom distriburtion from Task 2.3 and 2.4
        self.vel = np.random.rand(self.particle_count, 3) * 2.0 - 1

    def minimizeEnergy(self):
        result = minimize(
            Epot, 
            self.pos.ravel(),
            jac=grad_Epot,
            method='CG'            
        )

        new_pos = to3D(result.x)
        self.pos = new_pos
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


if __name__ == "__main__":
    domain = Domain()
    domain.fill(2, 1, 1)   
    print(domain.pos)

    old_energy = Epot(domain.pos))
    old_pos = 
    domain.minimizeEnergy()
    print(Epot(domain.pos))