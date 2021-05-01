#!/usr/bin/env python3
#include parent folder
import os, sys, inspect
currentdir = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import time
import math
import numpy as np
numpy = np    # to be compatible with
import scipy
from scipy.optimize import minimize
import matplotlib.pyplot as plt


import jax
import jax.numpy as jnp
from jax import grad, jit


from config import *
from src.misc import *


def Epot(positions):
    """Lennard-Jones potential in reduced units.
    In this system of units, epsilon=1 and sigma=2**(-1. / 6.).
    """
    positions = to3D(positions)
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


class Domain:
    # Domain settings
    particle_count = int  # Number of particles in the box
    length = float  # length of the box
    spread = float  # spread of randomizing factor of initial position
    std_dev = float  # standard deviation of the velocities
    pos = np.ndarray  # positional data of the particles
    vel = np.ndarray  # velocity data of the particle

    # Stats
    time_minimizeEnergy = float

    # Physics
    fEpot = object
    grad_Epot = object

    def __init__(self, fEpot=Epot):
        # Domain parameters
        self.particle_count = -1
        self.length = -1
        self.std_dev = -1

        # Domain Physics
        self.fEpot = jax.jit(fEpot)
        self.grad_Epot = jax.jit(jax.grad(self.fEpot))

        # Benchmark stats
        self.time_minimizeEnergy = -1

    def fill(self, particle_count, length, spread, std_dev):
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
        self.spread = spread
        self.std_dev = std_dev

        self.initialize_pos()
        self.initialize_vel()

        assert (self.pos.shape == self.vel.shape)

    def Epot(self):
        """Returns the potential Energy of the whole domain

        Returns:
            float: Potential Energy
        """        
        return float(self.fEpot(self.pos))

    def Ekin(self):
        """Returns the kinetic energy of the domain

        Returns:
            [type]: [description]
        """        
        return float(0.5 * np.dot(self.vel, self.vel.T).sum())

    def verlet_advance(self, dt):
        if(dt == 0): return

        # TODO: optimize, maybe switch to jax
        f = -self.grad_Epot(self.pos)        
        new_pos = self.pos + (self.vel + 0.5 * f * dt) * dt
        
        new_f = -self.grad_Epot(new_pos)
        new_vel = self.vel + 0.5 * (f + new_f ) * dt

        #new_pos = new_pos - self.length * np.round(new_pos/self.length)
        
        self.pos = new_pos 
        self.vel = new_vel 

    def initialize_pos(self):
        #dimensions of domain
        self.pos = numpy.ndarray((self.particle_count, 3)) #store positions
        nr_edge = math.ceil(self.particle_count**(1/3)) #particle nr per dimension/axis
        spacing = self.length/(nr_edge+1) #distance between neighbours on same axis

        #fill cubic domain
        nr_particles = 0
        for i in range(nr_edge): #z coordinate in domain
            for j in range(nr_edge): #y coordinate in domain
                for k in range(nr_edge): #x coordinate in domain

                    for c in range(3): #2nd axis of self.pos array
                        if(nr_particles < self.particle_count): #check part count
                            
                            #place atom like this, e.g. for 3 atoms:
                            # | spacing atom spacing atom spacing atom spacing |
                            # here, each atom is spaced by 1/5 of self.length
                            
                            if (c==0): #place atom 
                                self.pos[nr_particles][0] = k*spacing + spacing
                            elif (c==1):
                                self.pos[nr_particles][1] = j*spacing + spacing
                            elif (c==2):
                                self.pos[nr_particles][2] = i*spacing + spacing
                        else: break #after all particles are set: break
                    
                    #this iteration over nr_particles 
                    #cant be replaced by for-loop, for some reason
                    nr_particles += 1

        #randomize positions
        mean = [0, 0, 0]
        matrix = [[1, 0, 0], [0, 1, 0],[0,0,1]]
        randomize = np.random.multivariate_normal(mean, matrix, self.particle_count)*self.spread

        #old code
        # randomize = ((numpy.random.rand(self.particle_count, 3)-0.5)*self.spread)
        # print(randomize)
        self.pos += randomize

    def initizalize_pos_old(self):
        #this is the old initializer
        #places particles on space diagonal
        spacing = self.length / (self.particle_count)
        for i in range(self.particle_count):
            for j in range(3):
                self.pos[i][j] = spacing / 2 + i*spacing
                # + 1/2 to avoid placing particles at (0|0)
    
    def initialize_vel(self):
        self.vel = numpy.ndarray((self.particle_count, 3))
        vel = numpy.random.normal(0, self.std_dev, self.particle_count*3)

        # multivariate_normal .. method in numpy?
        
        for i in range(self.particle_count):
            for c in range(3): #2nd axis of self.vel array
                self.vel[i][c] = vel[0]
                vel = numpy.delete(vel,0)

        #set velocity mean to 0
        xmean = numpy.average(numpy.transpose(self.vel)[0])
        ymean = numpy.average(numpy.transpose(self.vel)[1])
        zmean = numpy.average(numpy.transpose(self.vel)[2])

        for i in range(self.particle_count):
            for c in range(3): #2nd axis of self.vel array
                if (c==0):
                    self.vel[i][c] -= xmean
                elif (c==1):
                    #print(c)
                    self.vel[i][c] -= ymean
                elif (c==2):
                    #print(c)
                    self.vel[i][c] -= zmean

    def initialize_vel_old(self):
        self.vel = numpy.random.rand(self.particle_count, 3) * 2.0 - 1

    def visualize_pos(self, show=True, fileName=""):
        x=[];y=[];z=[]
        for i in range(len(self.pos)):
            x.append(self.pos[i][0])
            y.append(self.pos[i][1])
            z.append(self.pos[i][2])
        
        fig = plt.figure(figsize = (10, 7))
        ax = plt.axes(projection ="3d")
        ax.scatter3D(x, y, z, color = "blue")
        # TODO: set dynamically to box size
        ax.set_xlim(0,self.length)
        ax.set_ylim(0,self.length)
        ax.set_zlim(0,self.length)
        plt.title("domain")

        if(show == True):
            plt.show()

        if(fileName != ""):
            plt.savefig(fileName)

    def minimizeEnergy(self):
        start = time.time()
        result = minimize(
            self.fEpot, 
            self.pos.ravel(),
            jac=self.grad_Epot,
            method='CG',
            #options={"disp" : True, "maxiter" : 10}  # Uncomment for Debug stats
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


def playgroud_reno():
    domain = Domain()
    domain.fill(10**3, 1, 0.1, 0.1)

    print(domain.pos)
    # print(domain.vel)
    #domain.visualize_pos()
    # domain.integrate(1,10)
    #domain.visualize_pos()

    """
    def integrate(self,dt,N):
        # integrates over N steps within a time inverval of lenght dt
        # force = grad(E) is set to 0 at the moment

        force = 0
        for i in range(N):
            self.pos = self.pos + self.vel*(dt/N) + 0.5*force*(dt/N)**2
            self.pos = self.vel + (force+force)*dt/(2*N)
            self.write_to_file("trajectory","testcomment")"""

        # it would be nice, if read_to_file() would add the new positions,
        # insted of overwriting the old ones

def playground_hiti():

    N = 20
    rho = 0.8

    L = ((float(N)/rho))**(1.0/3.0)
    np.random.seed(1)

    domain = Domain(Epot)
    domain.fill(10, 10, 1)   
    #domain.visualize_pos(show=False, fileName="out/plot1.png")
    domain.minimizeEnergy()
    #domain.visualize_pos(show=True, fileName="out/plot2.png")
    domain.verlet_advance(0.1)

    print(domain.Epot())
    print(domain.Ekin())

    #print(domain.vel)
    pos = domain.pos.T
    vel = domain.vel.T

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.quiver(pos[0], pos[1], pos[2], vel[0], vel[1], vel[2], length=3)
    #plt.show()

def playground_hickel():
    print("Grillenzirpen...")
    

if __name__ == "__main__":
    playgroud_reno()
    #playground_hickel()
    # playground_hiti()
