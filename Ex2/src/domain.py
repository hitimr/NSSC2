#include parent folder
import os, sys, inspect

currentdir = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import numpy
import jax.numpy as np
import jax
from jax import grad, jit
from config import *
from scipy.optimize import minimize


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
        self.pos = numpy.ndarray((self.particle_count, 3))

        spacing = self.length / (self.particle_count)
        for i in range(self.particle_count):
            for j in range(3):
                self.pos[i][j] = spacing / 2 + i*spacing # + 1/2 to avoid placing particles at (0|0)

    def initialize_vel(self):
        # TODO, Reno: replace with custom distriburtion from Task 2.3 and 2.4
        self.vel = numpy.random.rand(self.particle_count, 3) * 2.0 - 1

    def minimizeEnergy(self):
        grad_Epot = jit(grad(self.E_pot), static_argnums=(1,1))
        result = minimize(
            self.E_pot, 
            self.pos,
            method='CG'
            )
        new_pos = result.x
        new_pos.shape = (int(len(new_pos) / 3), 3)
        self.pos = new_pos

        return new_pos

    def E_pot(self, pos):
        # TODO: optimize this routine
        E = 0

        # pos has MxN dimensions
        if len(pos.shape) == 2:
            for i in range(len(pos)):
                for j in range(i+1):
                    E += self.V_LJ(np.linalg.norm(pos[i] - pos[j], 2))
            return E

        # pos is 1D array
        if len(pos.shape) == 1:
            pos.shape = ( int(pos.shape[0] / 3), 3)
            for i in range(len(pos)):
                for j in range(i+1):
                    E += self.V_LJ(np.linalg.norm(pos[i] - pos[j], 2))
            pos.shape = (pos.shape[0]*3)
        return E      


    # Potential in natural system of units
    def V_LJ(self, r):      
        # With a uniform distribution its possible that r is almost 0
        # So for now we just ignore that case

        return 4 * (pow(0.25 / r, 12) - pow(1 / r, 6))



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
    domain.fill(5, 1, 1)
    domain.write_to_file("test","bla")
    domain.read_from_file("test")


    print("generated positions:")
    print(domain.pos)
    print(f"Position of first particle:")
    print(domain.pos[0])
    
    """avg_xvel = 0
    avg_yvel = 0
    avg_zvel = 0
    
    for i in range(len(domain.vel)):
        avg_xvel += domain.vel[i][0]
        avg_yvel += domain.vel[i][1]
        avg_zvel += domain.vel[i][2]

    avg_xvel = avg_xvel/len(domain.vel)
    avg_yvel = avg_yvel/len(domain.vel)
    avg_zvel = avg_zvel/len(domain.vel)

    print("-------\n","xvel:",avg_xvel,"yvel:",avg_yvel,"zvel:",avg_zvel)"""
