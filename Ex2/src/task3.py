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
    """
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('integers', metavar='N', type=int, nargs='+',
                    help='an integer for the accumulator')
    parser.add_argument('--sum', dest='accumulate', action='store_const',
                    const=sum, default=max,
                    help='sum the integers (default: find the max)')
    """
    #print(sys.argv)
    # oder sys.arg
    #args = parser.parse_args()
    #print(args.accumulate(args.integers))



    domain = Domain(Epot)
    domain.fill(10,1,1)
    print(Epot(domain.pos))
    domain.verlet_advance(0.1)  # TODO: repeat N-1

    print(Epot(domain.pos))