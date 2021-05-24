import numpy
import numpy as np
import matplotlib.pyplot as plt

#https://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
def TDMAsolve(A, b):
    Ac = A.copy()
    bc = b.copy()
    n = A.shape[0]
    x = numpy.zeros(n)

    #eliminate lower diagonal
    for k in range(1, n):
        m = Ac[k,k-1]/Ac[k-1,k-1]
        Ac[k,k] = Ac[k,k] - m*Ac[k-1,k]
        bc[k] = bc[k] - m*bc[k-1]
        #Ac[k,k-1] = 0 #not necessary

    #reverse substitution
    x[n-1] = bc[n-1]/Ac[n-1,n-1]
    k=n-2
    for i in range(0, n-1):
        x[k] = (bc[k]-Ac[k,k+1]*x[k+1])/Ac[k,k]
        k=k-1

    return x


#fill discretisation matrix
def fill_matrix(matrix):
    for i in range(len(matrix[0])):
        matrix[i,i] = 1+2*S*D
    for i in range(len(matrix[0])-1):
        matrix[i+1,i] = -S*D
    for i in range(len(matrix[0])-1):
        matrix[i,i+1] = -S*D


#apply BC to concentration vector
def apply_BC(C_vector):
    C_vector[0] = C_vector[0] + S*D*dirichlet
    #C_vector[-1] = C_vector[-1] + S*D*(2*dx*neuman + C_vector[-2])
    #C_vector[-1] = C_vector[-2] #+ S*D*C_vector[-2]
    return C_vector





if __name__ == "__main__":
    #settings of time and space
    xmax = 100
    tmax = 1000
    dx = 1/xmax
    dt = 60
    D = 10**(-6)
    #D = 0.1
    #################################################
    dirichlet = 1 #at x=0
    neuman = 0 #at x=h

    #initialization
    space = np.linspace(1,xmax,xmax)
    time = np.linspace(0,tmax,tmax)
    S = dt/(dx*dx)
    C = [] #concentration values
    A = numpy.zeros([xmax,xmax]) #discretisation matrix
    C_0 = numpy.zeros([xmax]) #concentration vector at given time step

    #solve
    fill_matrix(A)
    #C_0 = apply_BC(C_0)
    C_old = C_0
    for i in range(tmax):
        C_old = apply_BC(C_old)
        #C_new = np.linalg.solve(A,C_old) #compute next C vector
        C_new = TDMAsolve(A,C_old)
        C.append(C_new)
        C_old = C_new #set old C vector to new for next iteration

    C = np.array(C)


    #plots
    plottimesteps = [0,round((tmax-1)*0.01),round((tmax-1)*0.1),round((tmax-1)*0.5),tmax-1]
    fig = plt.figure(figsize = (7,5))
    for i in plottimesteps:
        #plt.plot(space,C[i,:],label="C")
        plt.plot(space,C[i,:],label="C at t="+str(i))

    plt.legend()
    plt.xlabel("x")
    plt.ylabel("C")
    plt.grid()
    plt.show()
    #plt.savefig('plots/asdfasdfasdfasdf.pdf')