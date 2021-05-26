import numpy as np
import matplotlib.pyplot as plt
import argparse
from misc import *
import animate
import task1_functions as t1f
import plot



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to generate trajectory')
    parser.add_argument('Nx', type=int,help='Nx', default=100, nargs='?', const=1)
    parser.add_argument('deltat', type=float, help='timestep', default=10, nargs='?', const=1)
    parser.add_argument('num_iter', type=int, help='number of timesteps', default=10000, nargs='?', const=1)
    parser.add_argument('printmod', type=int, help='modulo for printing', default=1, nargs='?', const=1)
    parser.add_argument('filename', type=str, help='enter filename', default="task1_1", nargs='?', const=1)

    args = parser.parse_args()

    C=t1f.run(args.Nx,args.deltat,1,args.num_iter,args.printmod,args.filename)
    plot.setup()
    plt.title("Task 1_1")
    plt.plot(np.linspace(0,1, 100),C[-1])
    plt.savefig(DIR_OUT +args.filename+".png")
    #plt.show()
