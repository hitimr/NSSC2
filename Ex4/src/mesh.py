#include parent folder
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import numpy as np
import matplotlib.pyplot as plt

from src.magicsolver import *
from src.misc import *

class Mesh:
    variation = "noVariation"
    size = -1

    index_mat = np.array    # martirx that represents the index position of every node
    adj_mat = np.array  # adjacency matrix
    stiff_mat = np.array # global stiffness matrix

    nodal_temps = []
    nodal_forces = []
    nodal_coords_x = np.array
    nodal_coords_y = np.array


    

    def __init__(self, variation=None):

        # common settings for all variations
        file_params = read_from_file(DIR_SRC + "inputfile_group_5.txt")

        self.nx = 10
        self.ny = 10
        self.n = self.nx*self.ny
        self.L = float(file_params["L"])
        self.hz = float(file_params["hz"])
        self.k = float(file_params["k"])
        self.k_mod = float(file_params["c"])
        self.num_faces = (self.nx - 1) * (self.ny - 1) * 2
        


            #self.nodal_forces.append(float(file_params["q_y_L"]))
        #for i in range(90,self.n): self.nodal_forces.append(10)

        # initialize depending on the variation
        if variation == "V0": 
            self.init_V0()
        elif variation == "V1":
            self.init_V1()
        elif variation == "V3":
            self.init_V3()
        elif variation == None:
            pass
        elif variation == 'debug':
            self.init_V0()
        else: 
            raise ValueError(f"Unknown argumen '{variation}' for variation")

        self.variation = variation


        # apply boundary conditions
        # append is not ideal but size is only 100
        self.nodal_temps = []
        for i in range(0, 10): self.nodal_temps.append(float(file_params["T_y_0"]))     # • N1 – N10: Dirichlet BCs (depending on the “load case”) T = ...
        for i in range(10, 100): self.nodal_temps.append(None)

        self.nodal_forces = []
        for i in range(0,10): self.nodal_forces.append(None)  
        for i in range(10, 90): self.nodal_forces.append(0)     # • N1 – N10: Dirichlet BCs (depending on the “load case”) T = ...
        # cross section are recangles. lenghts may vary. height is constant
        for i in range(90, 100): 
            if(i == 90): 
                 # left most node
                dx = (self.nodal_coords_x[91] - self.nodal_coords_x[90]) / 2
            elif(i == 99):
                # right most node  
                dx = (self.nodal_coords_x[99] - self.nodal_coords_x[98]) / 2 
            else: 
                dx = (self.nodal_coords_x[i+1] - self.nodal_coords_x[i]) / 2 + (self.nodal_coords_x[i] - self.nodal_coords_x[i-1]) / 2

            # TODO: replace with Q from book (pdf p. 164)
            area = dx*self.hz
            nodal_power = -float(file_params["q_y_L"])*area
            self.nodal_forces.append(nodal_power)
            

        # generate data from init values
        self.index_mat = self.generate_nodal_indices()
        self.adj_mat = self.generate_adj_mat(self.nx*self.ny)
        self.stiff_mat = self.generate_global_stiffness_mat()
        return

    def get_face_area(self, face):
        nodes = self.get_face_nodes(face)

        # get absolute coordinates of every node
        y = [self.nodal_coords_x[node] for node in nodes]
        x = [self.nodal_coords_y[node] for node in nodes]

        b = [
            y[1] - y[2],
            y[2] - y[0],
            y[0] - y[1]
        ]

        c = [
            x[2] - x[1],
            x[0] - x[2],
            x[1] - x[0]
        ]

        area = (x[0] * b[0] + x[1]*b[1] + x[2]*b[2]) / 2
        return area

        

    def generate_nodal_indices(self):
        """Generate a matrix where every entry corresponds to the index of a givern node

        N91 N92 N93 ... N99
        ...
        N10 N11 N12 ... N19
        N00 N01 N02 ... N09

        Returns:
            [type]: [description]
        """        
        indices = np.arange(0, self.nx*self.ny, step=1)
        indices.shape = (self.nx, self.ny)
        return indices.T

    def get_face_nodes(self, face):
        # TODO: generalize for arbitrary nx, ny if there is time left

        # treat mesh as 2 triangles forming a square
        row = int((face-1)/18)
        col = int((face-1) % 18 / 2)

        #node1 = self.index_mat[col+1][row]
        #node2 = self.index_mat[col][row+1]

        #diagonal nodes
        if (face % 2 == 0):
            # upper triangle
            node1 = self.index_mat[col+1][row]
            node2 = self.index_mat[col+1][row+1]
            node3 = self.index_mat[col][row+1]
        else:
            node1 = self.index_mat[col][row]
            node2 = self.index_mat[col+1][row]
            node3 = self.index_mat[col][row+1]

        # return as list in ascending order
        nodes = [node1, node2, node3]
        #nodes.sort()
        return np.array(nodes)

    def generate_global_stiffness_mat(self):
        global_stiff_mat = np.zeros((self.n, self.n))

        for face in range(1, self.num_faces+1):
            nodes = self.get_face_nodes(face)

            element_stiffness_mat = self.generate_element_stiffness_mat(face)            
            # connectivity matrix
            # TODO replace witzh hand code 
            connectivity_mat = np.zeros((3,self.n))
            connectivity_mat[0][nodes[0]] = 1
            connectivity_mat[1][nodes[1]] = 1
            connectivity_mat[2][nodes[2]] = 1


            global_stiff_mat += connectivity_mat.T @ element_stiffness_mat.T @ connectivity_mat

            # for i in range(3):
            #     global_stiff_mat[nodes[i]][nodes[0]] += element_stiffness_mat[i][0]
            #     global_stiff_mat[nodes[i]][nodes[1]] += element_stiffness_mat[i][1]
            #     global_stiff_mat[nodes[i]][nodes[2]] += element_stiffness_mat[i][2]

        return global_stiff_mat  
        
    def generate_element_stiffness_mat(self, face):   # TODO: repalce with nodes to save redundant function call
        assert(face > 0)
        assert(face < (self.num_faces + 1))

        # TODO: check for kmod
        k = self.k
        nodes = self.get_face_nodes(face)

        # get absolute coordinates of every node
        x = [self.nodal_coords_x[node] for node in nodes]
        y = [self.nodal_coords_y[node] for node in nodes]

        # Zienkiewicz p. 120
        """
        a = [
            x[1]*y[2] - x[2]*y[1], 
            x[2]*y[0] - x[0]*y[2],
            x[0]*y[1] - x[1]*y[0]
        ]"""

        b = [
            y[1] - y[2],
            y[2] - y[0],
            y[0] - y[1]
        ]

        c = [
            x[2] - x[1],
            x[0] - x[2],
            x[1] - x[0]
        ]

        area = (x[0] * b[0] + x[1]*b[1] + x[2]*b[2]) / 2

        # calculate sub matrix
        H_xx = np.array([
            [b[0]*b[0], b[0]*b[1], b[0]*b[2]],
            [b[1]*b[0], b[1]*b[1], b[1]*b[2]],
            [b[2]*b[0], b[2]*b[1], b[2]*b[2]],
        ])

        H_yy = np.array([
            [c[0]*c[0], c[0]*c[1], c[0]*c[2]],
            [c[1]*c[0], c[1]*c[1], c[1]*c[2]],
            [c[2]*c[0], c[2]*c[1], c[2]*c[2]],
        ])
        H = np.zeros((3,3))
        H += k*self.hz / (4*area) * (H_xx + H_yy)
        return np.array(H)
        

    def init_V0(self):
        x = np.linspace(0.0, self.L, self.nx)
        y = np.linspace(0.0, self.L, self.ny)

        nodal_coords_x, nodal_coords_y = np.meshgrid(x,y)
        self.nodal_coords_x =  nodal_coords_x.ravel()
        self.nodal_coords_y =  nodal_coords_y.ravel()

    def init_V1(self):
        self.nx = 10
        self.ny = self.nx
        n =  self.nx * self.ny


        self.nodal_forces = np.zeros(n)
        self.nodal_forces[0:self.nx] = 1

        x = np.linspace(0, 1, self.nx)
        y = np.linspace(0, 1, self.ny)
        self.xv, self.yv = np.meshgrid(x,y)



    def init_V3(self):
        # values in polar coordinateL
        a, b = 1.25, 1
        r_vals = np.linspace(2*self.L, self.L, self.nx)
        phi_vals = np.linspace(a*np.pi, b*np.pi, self.ny) 

        # generate polar mesh
        r, phi = np.meshgrid(r_vals, phi_vals)


        x,y = [], []
        for i in range(len(r)):
            x.append(r[i]*np.cos(phi[i]))
            y.append(r[i]*np.sin(phi[i]))

        x = np.array(x).ravel() + 2*self.L
        y = np.array(y).ravel()

        # store in mesh  
        self.nodal_coords_x = np.array(x)
        self.nodal_coords_y = np.array(y)

    #returns adjacency matrix for regular mesh in Figure 1 (assignment)
    def generate_adj_mat(self, nnodes):

        if int(nnodes**0.5)*int(nnodes**0.5)!=nnodes: return "[ENTER NODE NUMBER OF SQUARE-SHAPED DOMAIN]"


        #mod_faces = [9,10,11,12,13,14,27,28,29,29,30,31,32]
        mod_faces = []

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

        #for i in range(nnodes):
            #A[i][i] = 1

        return A


    def solve(self):
        # solve Ax = b
        self.nodal_temps = np.linalg.solve(self.stiff_mat, self.nodal_forces)
        self.nodal_temps.shape = (self.nx, self.ny)

    def plot(self):
        T = np.reshape(T, (-1, len(x)))
        fig, ax = plt.subplots()
        plt.xlim(x[0], x[-1])
        plt.ylim(y[0], y[-1])
        ax.set_title('Temperature plot')
        ax.pcolormesh(T)
        plt.show()

    def plot_test(self, T):
        elements=[]
        for a in range(0,mesh.num_faces+1):
            elements.append(mesh.get_face_nodes(a))
        triangulation = tri.Triangulation(mesh.nodal_coords_x, mesh.nodal_coords_y,elements)
        plt.triplot(triangulation, '-k')
        plt.tricontourf(triangulation, T)
        plt.colorbar()
        plt.show()



if __name__ == "__main__":
    mesh = Mesh("V0")

    plt.matshow(mesh.stiff_mat)
    plt.show()
    T,P = magicsolver(mesh.stiff_mat, mesh.nodal_temps, mesh.nodal_forces)
    pass

