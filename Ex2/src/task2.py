#!/usr/bin/env python3
#include parent folder
import os, sys, inspect
currentdir = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from src.domain import Epot, Domain



if __name__ == "__main__":
    domain = Domain(Epot)   # Create domain. Epot (Argument) is the function that is used to calculate the potential
    domain.fill(8, 10, 1)  # Set size to 10 and fill it with 8 particles

    E = domain.Epot() + domain.Ekin()
    print(f"Total Energy = {E}")

    domain.minimizeEnergy() # Move particles to minimum
    
    E = domain.Epot() + domain.Ekin()
    print(f"Total Energy = {E}")

    domain.verlet_advance(0.01) # advance time by 0.01

    E = domain.Epot() + domain.Ekin()
    print(f"Total Energy = {E}")
