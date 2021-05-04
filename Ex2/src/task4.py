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


def part1():
    # TODO: Hickel
    pass


def part2():
    # TODO: Hiti
    # settings
    rho = 0.8
    N = 200
    L = (float(N/rho)**(1/3))
    dt = 0.01
    sigma = 1
    spread = 1/(3*N**(1/3))

    domain = Domain()
    domain.fill(N, L, spread, 1)
    print(domain.Epot())
    domain.minimizeEnergy()
    print(domain.Epot())




if __name__ == "__main__":
    part1()
    part2()
