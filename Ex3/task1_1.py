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
    
    num_iter=[100,500,1000,5000,10000,20000]
    C=t1f.run(args.Nx,args.deltat,1,args.num_iter,args.printmod,args.filename)
    plot.setup()
    plt.title("Task 1_1")
    plt.plot(np.linspace(0,1, 100),C[-1])
    plt.savefig(DIR_OUT +args.filename+".png")
    #plt.show()
    if args.filename=="task1_1":
        for b in num_iter:
            C=t1f.run(100,20,1,b,1,"Nx_100_dt_20_numIter_"+str(b)+"_task_1.txt")
        C=t1f.run(100,20,1,20000,1,"Nx_100_dt_20_numIter_20000_task_1.txt")
        C=t1f.run(100,20,2,20000,1,"Nx_100_dt_20_numIter_20000_task_2.txt")
        C=t1f.run(100,20,2,20000,1,"Nx_100_dt_20_numIter_20000_task_2.txt")
        C=t1f.run(100,1,1,10000,1,"Nx_100_dt_1_numIter_10000_task_1.txt")
        C=t1f.run(100,51,1,5000,1,"Nx_100_dt_51_numIter_5000_task_1.txt")
    #C=t1f.run(100,1,5,5000,1,"Nx_100_dt_51_numIter_5000_task_1.txt")
        plot.plot_num_iter_1_1()
        plot.plot_error()
        plot.plot_unstable()
        plot.plot_num_vs_analytical()
