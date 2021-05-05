#!/usr/bin/env python3
#include parent folder
import os, sys, inspect
currentdir = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import sys
import matplotlib.pyplot as plt

import argparse
from src.domain import Epot, Domain
from src.misc import *


def part1(infile):
    domain = Domain()
    f = open(infile, "r")
    all_lines=f.readlines()
    particle_count=float(all_lines[0])
    length=len(all_lines)
    dt_count=length/(particle_count+3)
    f.close()
    f = open(DIR_OUT + "task4_out1.txt", "w")
    for a in range(0,int(dt_count)):
        domain.read_from_file(infile,a)
        f.write(str(domain.Epot())+" "+str(domain.Ekin())+"\n")
    f.close()


def part2(FileName):
    # TODO: replace with trajectory file
    rho = 0.8
    N = 200
    L = (float(N/rho)**(1/3))
    dt = 0.01
    sigma = 1
    spread = 1/(3*N**(1/3))*0

    domain = Domain()
    domain.fill(N, L, spread, 1)
    domain.minimizeEnergy()

    samples = 200
    time_steps = 2000
    densities = []


    discarded_steps = time_steps // 4
    tracked_steps = time_steps - discarded_steps

    # Discard the first 25% by advancing the domain but not tracking results
    for i in range(discarded_steps):
        domain.verlet_advance(dt)
    
    # calculate volumetric densities
    for i in range(tracked_steps):
        bins = domain.volumetric_density(samples)
        densities.append(bins)
        domain.verlet_advance(dt)

    densities = np.array(densities)
    average_densities = densities.mean(axis=0)
    r = np.linspace(0, domain.length/2, len(average_densities))
    plt.plot(r, average_densities)
    plt.title(f"N={N}, sigma={sigma}")
    plt.yscale("log")
    plt.xlabel("Distance r")
    plt.ylabel("Particle density [N/V]")
    plt.savefig(DIR_OUT+"volumetric_density.png")
    plt.show()

    




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to generate trajectory')
    parser.add_argument('inputfile', type=str, help='name of inputfile', default=f'{DIR_OUT+"task3.txt"}', nargs='?', const=1)
    args = parser.parse_args()
    part1(args.inputfile)
    part2(args.inputfile)
