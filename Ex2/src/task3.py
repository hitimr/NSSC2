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
from misc import *


# see https://docs.python.org/3/library/argparse.html
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Script to generate trajectory')
    parser.add_argument('inputfile', type=str,help='name of inputfile', default=f'{DIR_OUT+"task2.txt"}', nargs='?', const=1)
    parser.add_argument('deltat', type=float, help='timestep', default=0.01, nargs='?', const=1)
    parser.add_argument('numiter', type=int, help='number of iterations', default=4, nargs='?', const=1)
 
    args = parser.parse_args()

    outFile = DIR_OUT + "task3.txt"
    domain = Domain(Epot)
    #domain.fill(10,1,1,0.1)
    domain.read_from_file(args.inputfile)
    domain.write_to_file(outFile, "comment",0)
    for a in range(0,args.numiter-1):
        domain.verlet_advance(args.deltat)
        domain.write_to_file(outFile, "comment",1)
    #domain.visualize_pos(True)
