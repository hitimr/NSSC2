#!/usr/bin/env python3
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
import math
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt


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
        #dimensions of domain
        self.pos = numpy.ndarray((self.particle_count, 3))
        nr_edge = math.ceil(self.particle_count**(1/3)) #set cube size
        spacing = self.length/nr_edge #distance between neighbours on same axis
        
        #fill cubic domain
        nr_particles = 0
        for i in range(nr_edge): #z coordinate in domain
            for j in range(nr_edge): #y coordinate in domain
                for k in range(nr_edge): #x coordinate in domain

                    for c in range(3): #2nd axis of self.pos array
                        if(nr_particles < self.particle_count): #check part countprint(nr_particles,c,"----",i,j,k)
                            if (c==0): #multiply by self.length missing
                                self.pos[nr_particles][0] = k/nr_edge + spacing/2
                            elif (c==1):
                                self.pos[nr_particles][1] = j/nr_edge + spacing/2
                            elif (c==2):
                                self.pos[nr_particles][2] = i/nr_edge + spacing/2
                        else: break #after all particles are set: break
                    
                    #this iteration over nr_particles 
                    #cant be replaced by for-loop, for some reason
                    nr_particles += 1

        #randomize positions
        randomize = (numpy.random.rand(self.particle_count, 3)-0.5)/5
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

    def visualize_pos(self):
        x=[];y=[];z=[]
        for i in range(len(self.pos)):
            x.append(self.pos[i][0])
            y.append(self.pos[i][1])
            z.append(self.pos[i][2])
        
        fig = plt.figure(figsize = (10, 7))
        ax = plt.axes(projection ="3d")
        ax.scatter3D(x, y, z, color = "blue")
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        ax.set_zlim(0,1)
        plt.title("domain")
        plt.show()

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

    def V_LJ(self, r):      
        # Potential in natural system of units
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
    domain.fill(27, 1, 0.1)
    #domain.visualize_pos()
    #domain.write_to_file("test","bla")
    #domain.read_from_file("test")
    print(domain.pos)
    print(domain.vel)