#!/usr/bin/env python3
#include parent folder
import os, sys, inspect
currentdir = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import sys

import argparse
from src.domain import Epot, Domain


# see https://docs.python.org/3/library/argparse.html
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Script to generate trajectory')
    parser.add_argument('inputfile', type=str,
                    help='name of inputfile')
    parser.add_argument('deltat', type=float, 
                    help='timestep')
    parser.add_argument('numiter', type=int, 
                    help='number of iterations')
 
    args = parser.parse_args()

    domain = Domain(Epot)
    #domain.fill(10,1,1,0.1)
    domain.read_from_file(args.inputfile)
    domain.write_to_file("outfile", "comment",0)
    for a in range(0,args.numiter-1):
        domain.verlet_advance(args.deltat)
        domain.write_to_file("outfile", "comment",1)
    #domain.verlet_advance(0.1)  # TODO: repeat N-1
    domain.visualize_pos(True)
    #print(Epot(domain.pos))
