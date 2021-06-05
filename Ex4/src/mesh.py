import numpy as np
import matplotlib.pyplot as plt


class Mesh:
    variation = "noVariation"
    size = -1
    adj_mat = np.array

    nodal_temps = np.array
    nodal_forces = np.array
    nodal_coords_x = np.array
    nodal_coords_y = np.array

    

    def __init__(self, variation):
        assert(type(variation) == str)

        if variation == "V0": 
            self.init_V0()
        elif variation == "V1":
            self.init_V1()
        elif variation == "V3":
            self.init_V3()
        else: 
            raise ValueError(f"Unknown argumen '{variation}' for variation")

        self.variation = variation

    def init_V0(self):
        """Initialize using the simplest case:

        2 - 3
        | \ |
        0 - 1

        adjacency matrix is hardcoded
            0   1   2   3   
        0   1   1   1   
        1   1   1        1
        2   1       1    1
        3       1   1    1

        """
        self.adj_mat = np.array([
            [1,1,1,0],
            [1,1,0,1],
            [1,0,1,1],
            [0,1,1,1],
        ])

        # Temperature distribution
        # 2.0   1.0
        # 2.0   1.0
        self.nodal_temps = [4.0, 2.0, 4.0, 2.0]
        
        pass

    def init_V1(self):
        self.nx = 5
        self.ny = self.nx
        n =  self.nx * self.nx
        k = -1

        self.adj_mat = self.generate_adj_mat(n)    
        upper = np.triu(self.adj_mat) 

        stiff_mat = np.zeros((n,n))
        for y in range(n):
            for x in range(n):
                if(upper[x][y] == 1):
                    sub_mat =  np.zeros((n,n))
                    sub_mat[x][x] = k
                    sub_mat[y][y] = k
                    #sub_mat[x][-y] = -k
                    #sub_mat[-x][y] = -k

                    stiff_mat += sub_mat

        for y in range(n):
            for x in range(n):
                sub_mat[x][y] = 1

        self.stiff_mat = stiff_mat
        self.nodal_forces = np.zeros(n)
        self.nodal_forces[0:self.nx] = 1

        x = np.linspace(0, 1, self.nx)
        y = np.linspace(0, 1, self.ny)
        self.xv, self.yv = np.meshgrid(x,y)

    def init_V3(self):
        self.nx = 10
        self.ny = 10
        self.L = 10

        # values in polar coordinates
        r_vals = np.linspace(self.L, 2*self.L, self.nx)
        phi_vals = np.linspace(0, 1/4*np.pi, self.ny) + np.pi

        # generate polar mesh
        r, phi = np.meshgrid(r_vals, phi_vals)

        # transform to cartesian
        x,y = [], []
        for i in range(len(r)):
            x.append(r[i]*np.cos(phi[i]))
            y.append(r[i]*np.sin(phi[i]))

        # store in mesh  
        nodal_coords_x = np.array(x)
        nodal_coords_y = np.array(y)

    #returns adjacency matrix for regular mesh in Figure 1 (assignment)
    def generate_adj_mat(self, nnodes):

        if int(nnodes**0.5)*int(nnodes**0.5)!=nnodes: return "[ENTER NODE NUMBER OF SQUARE-SHAPED DOMAIN]"

        n = int(nnodes**0.5)
        A = np.zeros([nnodes,nnodes])

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


    def solve(self):
        # solve Ax = b
        self.nodal_temps = np.linalg.solve(self.stiff_mat, self.nodal_forces)
        self.nodal_temps.shape = (self.nx, self.ny)

    def plot(self):
        plt.scatter(self.nodal_coords_x, self.nodal_coords_y)

if __name__ == "__main__":
    mesh = Mesh("V3")

    mesh.plot()