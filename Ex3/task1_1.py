import numpy as np
import matplotlib.pyplot as plt
import argparse
from misc import *
import animate
import task1_functions as t1f



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to generate trajectory')
    parser.add_argument('Nx', type=int,help='Nx', default=100, nargs='?', const=1)
    parser.add_argument('deltat', type=float, help='timestep', default=0.01, nargs='?', const=1)
    parser.add_argument('num_iter', type=int, help='number of timesteps', default=1000, nargs='?', const=1)
    parser.add_argument('printmod', type=int, help='modulo for printing', default=10, nargs='?', const=1)
    #parser.add_argument('task', type=int, help='task', default=1, nargs='?', const=1)
    parser.add_argument('filename', type=str, help='enter filename', default="default_filename.txt", nargs='?', const=1)
    #parser.add_argument('plot', type=int, help='plot or not to plot', default=0, nargs='?', const=1)

    args = parser.parse_args()

    #print(args.task)
    C=t1f.run(args.Nx,args.deltat,1,args.num_iter,args.printmod,args.filename)
    #plotsingle(C,"bla",args.filename)
    #t=10000
    #result_analitic = analytical_sol(t, args.Nx, 100)
    #plt.plot(np.linspace(0,1, args.Nx), result_analitic, label="analytic")
    #plt.plot(np.linspace(0,1, args.Nx), C[t],marker='o',linestyle = 'None', label="numeric")
    #plt.plot( result_analitic, label="analytic")
    #plt.plot(np.abs(result_analitic-C[t]))
    #plt.title(f"Numerical Error at t={t}, d=D*dt/dx² = {d}")
    #plt.legend()
    #plt.show()

    #fgr,axs=plt.subplots(2)
    #C1=run(100,2,1,1000)
    #C2=run(100,2,2,1000)
    #C3=run(100,2,1,100000)
    #C4=run(100,2,2,100000)
    #plt.yscale("log")
    #plt.plot(C3)
    #axs[0].plot(C1)
    #axs[0].plot(C2)
    #axs[1].plot(C3)
    #axs[1].plot(C4)
