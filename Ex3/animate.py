import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import argparse


def animate(ydata, fileName, title=""):
    xdata = np.linspace(0,1,len(ydata[0]))
    fig, ax = plt.subplots()
    line, = ax.plot(xdata, ydata[0], label="time evolution")
    initial_line = ax.plot(xdata, ydata[0], "--", label="initial condition")
    ax.legend(loc="upper right")
    ax.set_ylim([np.min(ydata), 1.1*np.max(ydata)])
    ax.set_title(title)
    global i
    i = 0

    def animate(i):
        line.set_ydata(ydata[i])
        i = i+1
        return line,


    ani = animation.FuncAnimation(fig, animate, frames=len(ydata))
    ani.save(fileName+".gif", fps=30, writer='imagemagick')

    fig.savefig(fileName+".png")

    return ani



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Animate a data File')
    parser.add_argument("inFile", type=str)
    parser.add_argument("outFile", type=str)
    args = parser.parse_args()

    data = np.loadtxt(args.inFile)
    animate(data, args.outFile)