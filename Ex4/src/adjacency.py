import math
import numpy


#position of finite element in 9x9 mes
def posFE(FE):
    posFEx = math.ceil(((FE-2)%18)/2)
    posFEy = math.ceil((FE-2)/18-1)
    posrel = "lower"
    if FE%2 == 0: posrel = "upper"
    return [posFEx, posFEy, posrel]


#returns adjacency matrix for regular mesh in Figure 1 (assignment)
def generate_adjacency(nnodes):

    if int(nnodes**0.5)*int(nnodes**0.5)!=nnodes: return "[ENTER NODE NUMBER OF SQUARE-SHAPED DOMAIN]"

    n = int(nnodes**0.5)
    A = numpy.zeros([nnodes,nnodes])

    #fill matrix neglecting boundary nodes
    hshift = 1 #horizontal shift
    vshift = n #vertical shift
    dshift = n - 1 #diagonal shift from bottom right to upper left
    dshift2 = n + 1 #diagonal shift other direction (not used!)
    for i in range(nnodes-hshift): #horizontal
        A[i,i+hshift] = 1
        A[i+hshift,i] = 1
    for i in range(nnodes-vshift): #vertical
        A[i,i+vshift] = 1
        A[i+vshift,i] = 1
    for i in range(nnodes-dshift): #diagonal
        A[i,i+dshift] = 1
        A[i+dshift,i] = 1
    
    #correct boundary nodes
    hskip = n
    dskip = n
    for i in range(n,nnodes-hshift,hskip): #horizontal
        A[i-hshift,i] = 0
        A[i,i-hshift] = 0
    #vertikal muss nicht korrigert werden
    for i in range(n-1,nnodes,dskip): #diagonal
        A[i-dshift,i] = 0
        A[i,i-dshift] = 0

    return A


if __name__ == "__main__":
    print(posFE(FE=40))
    print(generate_adjacency(nnodes=16))