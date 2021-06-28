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

    def __init__(self, deprecated_Epot=None):
        """Initialize the domain

        Args:
            deprecated_Epot (any, optional): old way of adding Epot to the
            domain. No longer has any effect but its still in there to be
            backwards compatible
        """
        # Domain parameters
        self.particle_count = -1
        self.length = -1
        self.std_dev = -1

        # Domain Physics
        self.fEpot = jax.jit(self.Epot_source)
        self.grad_Epot = jax.jit(jax.grad(self.fEpot))

        # Benchmark stats
        self.time_minimizeEnergy = -1


    def Epot_source(self, positions):
        """Lennard-Jones potential in reduced units.
        In this system of units, epsilon=1 and sigma=2**(-1. / 6.).
        """
        positions = to3D(positions)
        if positions.ndim != 2 or positions.shape[1] != 3:
            raise ValueError("positions must be an Mx3 array")
        # Compute all squared distances between pairs without iterating.

        delta = positions[:, np.newaxis, :] - positions
        delta = delta - self.length * np.round(delta/self.length)
        r2 = (delta * delta).sum(axis=2)
        # Take only the upper triangle (combinations of two atoms).
        indices = np.triu_indices(r2.shape[0], k=1)
        rm2 = 1. / r2[indices]
        # Compute the potental energy recycling as many calculations as possible.
        rm6 = rm2 * rm2 * rm2
        rm12 = rm6 * rm6
        return (rm12 - 2. * rm6).sum()


    def fill(self, particle_count, length, spread, std_dev, minimizeEpot=False, verbose=False):
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

        if(verbose): print(f"Filling domain with {self.particle_count} particles")
        self.initialize_pos()   # Placing M particles at random positions in the box
        if(verbose): print(f"Epot = {self.Epot()}")

        # Move particles to local energy minimum using CG Thisfunction initially
        # was called outside of fill(). The parameter minimizeEpot exists so old
        # code is still compatible
        if(minimizeEpot == True):
            if(verbose): print("Moving particles to local energy minimum")
            self.minimizeEnergy()
            if(verbose): print(f"Epot = {self.Epot()}")

        # Add velocities
        if(verbose): print("Adding velocities")
        self.initialize_vel()
        if(verbose): print(f"Epot = {self.Epot()}, Ekin = {self.Ekin()}")
        
        # Sanity checks
        # TODO: check if the sum of all forces = 0
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
        return float(0.5 * np.dot(self.vel.T, self.vel).sum())

    def verlet_advance(self, dt):
        if(dt == 0): return

        f = -self.grad_Epot(self.pos)       
        new_pos = self.pos + (self.vel + 0.5 * f * dt) * dt

        new_f = -self.grad_Epot(new_pos)
        new_vel = self.vel + 0.5 * (f + new_f ) * dt

        
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
        matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
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
            method='CG'
        )
        new_pos = to3D(result.x)
        self.pos = new_pos

        end = time.time()
        self.time_minimizeEnergy = end - start
        return new_pos  

    def read_from_file(self, fileName, block=0):
        """Fill the domain with data from a file
        #TODO, Hickel: read n-th block 

        Args:
            fileName (str): path to file
        """
        f = open(fileName, "r")
        self.particle_count=int(f.readline())
        comment=f.readline()
        self.length=float(f.readline())
        self.pos = numpy.ndarray((self.particle_count, 3))
        self.vel = numpy.ndarray((self.particle_count, 3))
        if block==0:
            for a in range(0,self.particle_count):
                line_tmp=f.readline().split(" ")
                for b in range(0,3):
                    self.pos[a][b]=float(line_tmp[b])
                    self.vel[a][b]=float(line_tmp[b+3])
        else:
            all_lines=f.readlines()
            for a in range(0,self.particle_count):
                line_tmp=all_lines[block*(self.particle_count+3)+a].split(" ")
                for b in range(0,3):
                    self.pos[a][b]=float(line_tmp[b])
                    self.vel[a][b]=float(line_tmp[b+3])
        f.close()
        assert (self.particle_count > 0)
        assert (self.length > 0)
        assert (self.pos.shape == self.vel.shape)

        return
#'/home/hiti/repos/NSSC2/Ex2/out/unitTest_testFile.txt'
    def write_to_file(self, fileName, comment="", append=0):
        """Write domain data to a file

        #TODO, Hickel: Append Blocks

        Args:
            fileName (str): path to file
            comment (str): arbitrary comment/description for Line 2
        """
        if append==1:
            f = open(fileName, "a")            
        else:
            f = open(fileName, "w")
        f.write(str(self.particle_count)+"\n")
        f.write(comment+"\n")
        f.write(str(self.length)+"\n")
        for a in range(0,self.particle_count):
            f.write(str(self.pos[a][0])+" "+str(self.pos[a][1])+" "+str(self.pos[a][2])+" "+str(self.vel[a][0])+" "+str(self.vel[a][1])+" "+str(self.vel[a][2]))
            f.write("\n")
        f.close()    


        return

    def writeenergy(self, fileName, comment="", append=0):
        """Write domain data to a file

        #TODO, Hickel: Append Blocks

        Args:
            fileName (str): path to file
            comment (str): arbitrary comment/description for Line 2
        """
        if append==1:
            f = open(fileName, "a")
        else:
            f = open(fileName, "w")
        f.write(str(self.particle_count)+"\n")
        f.write(comment+"\n")
        f.write(str(self.length)+"\n")
        for a in range(0,self.particle_count):
            f.write(str(self.pos[a][0])+" "+str(self.pos[a][1])+" "+str(self.pos[a][2])+" "+str(self.vel[a][0])+" "+str(self.vel[a][1])+" "+str(self.vel[a][2]))
            f.write("\n")
        f.close()

    def average_distance(self):
        distances = np.triu(scipy.spatial.distance.cdist(domain.pos, domain.pos)) # calculate distances but throws away lower triangle
        average = distances.sum() / (2*self.particle_count**2)
        return average


    def volumetric_density(self, sample_cnt=20):

        r = np.linspace(0, self.length/2., sample_cnt)

        # set 0th particle as the origin
        origin = self.pos[0]

        # calculate eudlicdian distances to that particle
        distances = np.array(np.linalg.norm(self.pos[1:] - origin, axis=1))

        # perform minimum image transformation
        distances = abs(MI_transform(distances, self.length))
        assert(len(distances) == self.particle_count-1)

        # ust histogram() to count particles in each shell
        hist, edges = np.histogram(distances, bins=r)
        
        densities = []
        for i in range(1,len(hist)):
            # calculate volume of the shell
            V = 4/3 * math.pi * (r[i]**3 - r[i-1]**3)

            # divide number of particles in that
            densities.append(hist[i]  / V)

        return densities



def playgroud_reno():
    domain = Domain()
    domain.fill(10**3, 1, 0.1, 0.1)

    import pandas as pd
    N = 20
    rho = 0.8

    L = ((float(N)/rho))**(1.0/3.0)
    np.random.seed(1)

    domain = Domain(Epot)
    domain.fill(10, 10, 1, 1)   
    #domain.visualize_pos(show=False, fileName="out/plot1.png")
    domain.minimizeEnergy()
    #domain.visualize_pos(show=True, fileName="out/plot2.png")
    for i in range(10):
                
        domain.verlet_advance(0.1)


        print(domain.Epot())
        print(domain.Ekin())



def playground_hiti():

    N = 200
    rho = 0.8

    L = ((float(N)/rho))**(1.0/3.0)
    np.random.seed(1)

    domain = Domain(Epot)
    domain.fill(N, L, 0, 1)  

    domain.volumetric_density()


def playground_hickel():
    domain = Domain()
    domain.fill(5, 9, 1, 0.1)
    domain.write_to_file("test", "test")
    domain.write_to_file("test", "test",1)
    domain_new = Domain()
    domain_new.read_from_file(DIR_OUT+"test",0)
    

if __name__ == "__main__":
    #playgroud_reno()

    playground_hickel()
    #playground_hiti()
