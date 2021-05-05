#!/usr/bin/env python3
#include parent folder
import os, sys, inspect
currentdir = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import argparse
from src.domain import Epot, Domain
from src.misc import *



def task2(M, L, sigma, fileName):

    spread = 1/(3*M**(1/3)) # Parameter for how much particles are distrubed from their ideal position


    print(f"Creating domain with M={M}, L={L}, Σ={sigma}")
    domain = Domain()
    domain.fill(M, L, spread, sigma, minimizeEpot=True, verbose=True)

    print(f"Domain Setup complete. Saving results to {fileName}")
    domain.write_to_file(fileName, comment="Initial conditions",append=False)




if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(description='Script to generate initiual conditions ')
        parser.add_argument('M', type=int, help='Number of aprticles')
        parser.add_argument('L', type=float, help='side lenght of the simulation box')
        parser.add_argument('sigma', type=float, help='standard deviation')
        parser.add_argument('fileName', type=str, help='soutput file')

        args = parser.parse_args()


        M = args.M
        L = args.L
        sigma = args.sigma
        fileName = args.fileName

    except:
        print("Error while parsing parameters. Using default parameters instead")
        M = 200
        L = 6.3
        sigma = 1
        fileName = DIR_OUT + "task2.txt"

    task2(M, L, sigma, fileName)
