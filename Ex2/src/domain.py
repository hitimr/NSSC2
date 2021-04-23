#include parent folder 
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) 
parentdir = os.path.dirname(currentdir) 
sys.path.insert(0,parentdir) 
 

import numpy as np
from config import *


class Domain:
    particle_count = int    # Number of particles in the box
    length = float          # length of the box
    std_dev = float         # standard deviation of the velocities
    pos = np.ndarray        # positional data of the particles
    vel = np.ndarray        # velocity data of the particles

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
        assert(isinstance(particle_count, int))
        assert(particle_count > 0)
        assert(length > 0)
        assert(std_dev > 0)

        self.particle_count = particle_count
        self.length = length
        self.std_dev = std_dev  

        self.pos = np.random.rand(self.particle_count, 3)   # TODO: replace with custom distribution. see description
        self.vel = np.random.rand(self.particle_count, 3)   # TODO: replace gaussian distribution


    def read_from_file(self, fileName):
        """Fill the domain with data from a file

        Args:
            fileName (str): path to file
        """

        return

    def write_to_file(self, fileName):
        """Write domain data to a file

        Args:
            fileName (str): path to file
        """       

        return


        
if __name__ == "__main__":
    domain = Domain()
    domain.fill(10, 1, 1)

    print("generated positions:")
    print(domain.pos)
    print(f"Position of first particle:")
    print(domain.pos[0])

    print("\ngenerated velocities:")
    print(domain.vel)
    print(f"Velocity of second particle:")
    print(domain.vel[1])