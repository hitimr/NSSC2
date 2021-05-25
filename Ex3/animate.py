import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import argparse


def animate(yNumeric, fileName, title="no title", yAnalytic=[]):
    fig, ax = plt.subplots()
    

    xNumeric =  np.linspace(0,1,len(yNumeric[0]))
    lineNumeric, = ax.plot(xNumeric, yNumeric[0], label="numeric evolution")

    # Add analytic data if passed as argument
    if len(yAnalytic) != 0:
        if(len(yNumeric) != len(yAnalytic)): raise ValueError("numeric solution and analytic solution must have the same number of time steps")

        xAnalytic = np.linspace(0,1,len(yAnalytic[0])) 
        lineAnalytic, = ax.plot(xAnalytic, yAnalytic[0], "x", label="analytic evolution")
    
    global i
    i = 0

    initial_line = ax.plot(xNumeric, yNumeric[0], "--", label="initial condition")
    ax.legend(loc="upper right")
    ax.set_ylim([np.min(yNumeric), 1.1*np.max(yNumeric)])
    ax.set_title(f"{title}, Frame: {i}")

    def update(i):
        i += 1
        # Update Title
        ax.set_title(f"{title}, Frame: {i}")

        # Update Data
        lineNumeric.set_ydata(yNumeric[i-1])

        if len(yAnalytic) != 0:
            lineAnalytic.set_ydata(yAnalytic[i-1])


    ani = animation.FuncAnimation(fig, update, frames=len(yNumeric))
    ani.save(fileName+".gif", fps=30, writer='imagemagick')

    fig.savefig(fileName+".png")

    return ani



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Animate a data File')
    parser.add_argument("inFile", type=str)
    parser.add_argument("outFile", type=str)
    args = parser.parse_args()

    data = np.loadtxt(args.inFile)
    animate(data, args.outFile, title="Task2")