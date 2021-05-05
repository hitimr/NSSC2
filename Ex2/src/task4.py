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


def part1(infile):
    # TODO: Hickel
    domain = Domain()
    f = open(infile, "r")
    all_lines=f.readlines()
    particle_count=float(all_lines[0])
    length=len(all_lines)
    dt_count=length/(particle_count+3)
    #print(dt_count)
    #print(len(all_lines))
    f.close()
    f = open("task4_out1.txt", "w")
    for a in range(0,int(dt_count)):
        domain.read_from_file(infile,a)
        #E_pot=domain.Epot()
        #E_kin=domain.Ekin()
        f.write(str(domain.Epot())+" "+str(domain.Ekin())+"\n")
    f.close()
    #pass


def part2():
    # TODO: Hiti
    # TODO: replace with trajectory file
    rho = 0.8
    N = 200
    L = (float(N/rho)**(1/3))
    dt = 0.01
    sigma = 1
    spread = 1/(3*N**(1/3))

    domain = Domain()
    domain.fill(N, L, spread, 1)
    domain.minimizeEnergy()

    




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to generate trajectory')
    parser.add_argument('inputfile', type=str,
                    help='name of inputfile')
    args = parser.parse_args()
    part1(args.inputfile)
    #part2()
