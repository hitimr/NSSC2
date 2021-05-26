import numpy as np
import matplotlib.pyplot as plt

from animate import animate

def func_gauss(x):
    return np.exp(-10.*np.power(4*x-1,2))

def func_square(x):
    if hasattr(x, "__len__"): return [func_square(xi) for xi in x]
    if(0.1 < x) and (x < 0.3): return 1.
    else: return 0.



def enforce_pb(x, L=1.):
    """Enforce Periodic Boundary Conditions using minimum image convention

    Args:
        x (double): coordinate
        L (double, optional): size of domain. Defaults to 1.

    Returns:
        (double): mapped coodrdinate
    """    
    return x - L * np.round(x/L)*1.01

def advect_analytical(initFunc, Nx, steps, c=1.):
    x_vals = np.linspace(0,1, Nx)
    return np.array([initFunc(enforce_pb(x_vals - c*t)) for t in np.linspace(0,1,steps)])

def advance_upwind(C, Co):
    new_C = np.copy(C)
    for i in range(len(C)):
        new_C[i] = (1 - Co) * C[i] + Co * C[i-1]
    return new_C

def advect_upwind(initFunc, Co, steps, Nx):
    x_vals = np.linspace(0,1, Nx)
    C = initFunc(x_vals)

    res = [C]
    for i in range(steps-1):
        C = advance_upwind(C, Co)
        res.append(C)

    res = np.array(res)
    return res


def animations():
    dx = 0.01
    Nx = int(1/(dx+0))
    x_vals = np.linspace(0, 1, Nx)
    steps = 130

    for co in [0.8, 1.0, 1.05]:
        print(f"Animating for Co={co}")
        res_analytic = advect_analytical(func_gauss, Nx=Nx, steps=steps, c=dx*steps*co)
        res_numeric =  advect_upwind(func_gauss, co, steps, Nx)
        animate(res_numeric, f"animations/task2_gauss_co_{str(co)}", title=f"Courant = {str(co)}", yAnalytic=res_analytic)


    #plt.show()


def courantPlots(courant):
    dx = 0.01
    Nx = int(1/(dx))-1
    x_vals = np.linspace(0, 1, Nx)
    steps = 100


    rowCnt = 3
    colCnt = 2
    f, (ax1, ax2, ax3) = plt.subplots(rowCnt, colCnt, sharey=True, figsize=(7, 5))
    f.subplots_adjust(hspace=0.5)
    f.suptitle(f"Co = {courant}")
    
    initFunc = [func_gauss, func_square]
    #initFunc = [func_gauss]

    for i in range(len(initFunc)):
        res_analytic = advect_analytical(initFunc[i], Nx=Nx, steps=steps, c=dx*steps*courant)
        res_numeric =  advect_upwind(initFunc[i], courant, steps, Nx)

        markersize = 3.5
        frames = [0, 25, 80]
        # (1,1)
        ax1[i].plot(x_vals, res_numeric[frames[0]], label="numeric")
        ax1[i].plot(x_vals, res_analytic[frames[0]], "x", label="analytic", markersize=markersize)
        ax1[i].set_title(f'Frame: {frames[0]}')
        ax1[i].set_ylim([0, 1.2])
        ax1[i].legend()

        ax2[i].plot(x_vals, res_numeric[frames[1]], label="numeric")
        ax2[i].plot(x_vals, res_analytic[frames[1]], "x", label="analytic", markersize=markersize)
        ax2[i].set_title(f'Frame: {frames[1]}')
        

        ax3[i].plot(x_vals, res_numeric[frames[2]], label="numeric")
        ax3[i].plot(x_vals, res_analytic[frames[2]], "x", label="analytic", markersize=markersize)
        ax3[i].set_title(f'Frame: {frames[2]}')

    plt.savefig(f"out/courant_{str(courant)}.png")
    #plt.show()


def accuracyPlot():
    dx = 0.01
    Nx = int(1/(dx))
    x_vals = np.linspace(0, 1, Nx)
    steps = 100


    res_analytic_gauss = advect_analytical(func_gauss, Nx=Nx, steps=steps, c=dx*steps)
    res_numeric_gauss =  advect_upwind(func_gauss, 1, steps, Nx)


    res_analytic_square = advect_analytical(func_square, Nx=Nx, steps=steps, c=dx*steps)
    res_numeric_square =  advect_upwind(func_square, 1, steps, Nx)

    # Single frame Error
    t = 25
    plt.plot(x_vals, np.abs(res_analytic_gauss[t] - res_numeric_gauss[t]), label="Error gauss")
    plt.plot(x_vals, np.abs(res_analytic_square[t] - res_numeric_square[t]), label="Error square")
    #plt.yscale("log")
    plt.title(f"Absolute Error for Frame {t}")
    plt.xlabel("x (admim)")
    plt.ylabel("Absolute Error (admin)")
    plt.legend()
    plt.savefig("out/courant_1_error.png")


if __name__ == "__main__":
    accuracyPlot()
    animations()
    courantPlots(0.8)
    courantPlots(1.0)
    courantPlots(1.2)
    