import numpy as np
import matplotlib.pyplot as plt
import argparse
import task1_functions as tf1
from misc import *

def setup():
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('x')
    ax.set_ylabel('C')
    ax.legend(loc='best')
    plt.legend(loc="upper right")


def plot_num_iter():
    num_iter=[100,500,1000,5000,10000,20000]
    for a in range(1,3):
        fgr,axs=plt.subplots(1)
        axs.set_xlabel('x')
        axs.set_ylabel('C')
        axs.legend(loc='best')
        for b in num_iter:
            C_tmp=np.loadtxt(DIR_OUT + "Nx_100_dt_20_numIter_"+str(b)+"_task_"+str(a)+".txt")
            axs.plot(np.linspace(0,1, 100),C_tmp[-1],label="num_iter = "+str(b))
            plt.title("C for different number of iterations task1."+str(a))
        plt.legend(loc="upper right")
        plt.savefig(DIR_OUT +"task1_"+str(a)+"numIter_task.png")
        #plt.show()

def plot_num_vs_analytical():
    setup()
    t=10000
    result_analitic = tf1.analytical_sol(t, 100, 100)
    C=np.loadtxt(DIR_OUT + "Nx_100_dt_1_numIter_10000_task_1.txt")
    plt.title("Numerical vs analytical solution")
    plt.plot(np.linspace(0,1, 100), result_analitic, label="analytic")
    plt.plot(np.linspace(0,1, 100), C[-1],marker='o',linestyle = 'None', label="numeric")
    plt.legend()
    plt.savefig(DIR_OUT + "task_1_1_num_vs_analytical.png")
    #plt.show()

def plot_error():
    setup()
    t=10000
    C=np.loadtxt(DIR_OUT + "Nx_100_dt_1_numIter_10000_task_1.txt")
    plt.title("Absolute error numerical and analytical solution")
    result_analitic = tf1.analytical_sol(t, 100, 100)
    plt.plot(np.linspace(0,1, 100),np.abs(result_analitic-C[-1]))
    plt.savefig(DIR_OUT + "task_1_1_error.png")
    #plt.show()

def plot_t_infinity():
    setup()
    C1=np.loadtxt(DIR_OUT + "Nx_100_dt_20_numIter_20000_task_1.txt")
    C2=np.loadtxt(DIR_OUT + "Nx_100_dt_20_numIter_20000_task_2.txt")
    plt.title("Behavior for large t")
    plt.plot(np.linspace(0,1, 100), C1[-1], label="analytic")
    plt.plot(np.linspace(0,1, 100), C2[-1], label="numeric")
    plt.legend()
    plt.savefig(DIR_OUT + "task_1_2_t_infinity.png")
    #plt.show()

def plot_unstable():
    setup()
    C=np.loadtxt(DIR_OUT + "Nx_100_dt_51_numIter_5000_task_1.txt")
    plt.title("Unstable solution for explicit solution")
    plt.plot(np.linspace(0,1, 100), C[-1], label="numeric unstable")
    plt.legend()
    plt.savefig(DIR_OUT + "task_1_1_unstable.png")
    #plt.show()


if __name__ == "__main__":

    plot_num_iter()
    plot_num_vs_analytical()
    plot_error()
    plot_t_infinity()
    plot_unstable()
