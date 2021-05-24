import numpy as np
import matplotlib.pyplot as plt
import argparse
from misc import *


def apply_boundaries(C,task):
    if task == 1:
        C[0] = 1
        C[-1] = C[-2]
    if task == 2:
        C[0] = 1
        C[len(C)-1] = 0
    return C


def advance(C,task,D,S):
    new_C = np.empty(len(C))
    for i in range(1, len(C)-1):
        new_C[i] = C[i] + D*S*(C[i-1] - 2*C[i] + C[i+1])
    apply_boundaries(new_C,task)
    return new_C



def run(Nx,dt,task,num_iter,printmod,filename):
    dx = (1/(Nx+1))
    S = dt/(dx*dx)
    D = 10**(-6)
    d = D*dt/dx**2

    C = np.zeros(Nx)
    C = apply_boundaries(C,task)
    res = [C]

    d=D*dt/dx**2
    #print(d)

    for i in range(int(num_iter)):
        C = advance(C,task,D,S)
        #res.append(C)
        if i % printmod == 0:
            res.append(C)
    np.savetxt(DIR_OUT + filename,res)

    #print(res)
    return C

def plotsingle(C,desc,filename):
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('x')
    ax.set_ylabel('C')
    ax.legend(loc='best')
    plt.plot(C,label=str(desc))
    plt.legend(loc="upper right")
    plt.savefig(DIR_OUT + filename + ".png")
    #plt.show()



parser = argparse.ArgumentParser(description='Script to generate trajectory')
parser.add_argument('Nx', type=int,help='Nx', default=100, nargs='?', const=1)
parser.add_argument('deltat', type=float, help='timestep', default=0.01, nargs='?', const=1)
parser.add_argument('num_iter', type=int, help='number of timesteps', default=1000, nargs='?', const=1)
parser.add_argument('printmod', type=int, help='modulo for printing', default=10, nargs='?', const=1)
parser.add_argument('task', type=int, help='task', default=1, nargs='?', const=1)
parser.add_argument('filename', type=str, help='enter filename', default="default_filename.txt", nargs='?', const=1)
#parser.add_argument('plot', type=int, help='plot or not to plot', default=0, nargs='?', const=1)

args = parser.parse_args()

#print(args.task)
C=run(args.Nx,args.deltat,args.task,args.num_iter,args.printmod,args.filename)
plotsingle(C,"bla",args.filename)


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
