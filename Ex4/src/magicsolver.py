import numpy

#for the equation A*T = P (A matrix, T, P vectors)
#we know some values from T and some from P, s.t. there exists a unique sol
def magicsolver(A,T,P):
    n = len(A)
    T_known = []
    T_known_i = []
    T_unknown = []
    T_unknown_i = []
    P_known = []
    P_known_i = []
    P_unknown = []
    P_unknown_i = []
    
    #find unkown values
    for i in range(n):
        if T[i] == None:
            T_unknown.append(T[i])
            T_unknown_i.append(i)
        else:
            T_known.append(T[i])
            T_known_i.append(i)
    for i in range(n):
        if P[i] == None:
            P_unknown.append(P[i])
            P_unknown_i.append(i)
        else:
            P_known.append(P[i])
            P_known_i.append(i)
    
    #compute T_unknown
    A_sub = A[P_known_i[0]:P_known_i[-1]+1, P_known_i[0]:P_known_i[-1]+1]

    for Pi in range(len(P_known)):
        for Ti in range(len(T_known)):
            P_known[Pi] = P_known[Pi] - A[P_known_i[Pi],T_known_i[Ti]] * T_known[Ti]
    
    T_unknown = numpy.linalg.solve(A_sub, P_known)

    j = 0
    for i in range(n):
        if T[i] == None:
            T[i] = T_unknown[j]
            j+=1

    #compute P_unknown
    for i in range(len(P_unknown)):
        P_unknown[i] = numpy.dot(A[P_unknown_i[i]],T)

    j = 0
    for i in range(n):
        if P[i] == None:
            P[i] = P_unknown[j]
            j+=1

    return T,P


if __name__ == "__main__":
    A = numpy.array([[9,10,11,12],[1,2,3,4],[5,6,7,8],[13,14,15,16]])
    T = [100,None,None,400]
    P = [None,2000,3000,None]

    T,P = magicsolver(A,T,P)
    print("T =",T)
    print("P =",P)