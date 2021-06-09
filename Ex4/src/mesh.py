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

    face_ids = np.array # number of every face [1,2,3,...] as stated in the assignment
    face_k = np.array   # k value for a given face_id
    face_center_x = np.array    # x coordinate of the center of every face
    face_center_y = np.array    # y coordinate of the center of every face
    face_flux_x = np.array  # x component of the flux
    face_flux_y = np.array  # y component of the flu
    face_gradient_x = np.array
    face_gradient_y = np.array

    def __init__(self, variation=None):

        # common settings for all variations
        self.file_params = read_from_file(DIR_SRC + "inputfile_group_5.txt")

        self.nx = 10
        self.ny = 10
        self.n = self.nx*self.ny
        self.L = float(self.file_params["L"])
        self.hz = float(self.file_params["hz"])
        self.k = float(self.file_params["k"])
        self.k_mod = float(self.file_params["c"])
        self.num_faces = (self.nx - 1) * (self.ny - 1) * 2
        #self.mod_faces = file_params["all_elements"]
        self.mod_faces = [-1]   # no faces to be modiefied
        

        # initialize depending on the variation
        if variation == "V0": 
            self.init_V0()
        elif variation == "V1":
            self.init_V1()
        elif variation == "V2":
            self.init_V2()
        elif variation == "V3":
            self.init_V3()
        elif variation == "V4a":
            self.init_V4a()
        elif variation == "V4b":
            self.init_V4b()           
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
        for i in range(0, 10): self.nodal_temps.append(float(self.file_params["T_y_0"]))     # • N1 – N10: Dirichlet BCs (depending on the “load case”) T = ...
        for i in range(10, 100): self.nodal_temps.append(None)

        # define heat sink
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
            self.nodal_power = -float(self.file_params["q_y_L"])*area
            self.nodal_forces.append(self.nodal_power)
            

        # generate matrices for simulation
        self.index_mat = self.generate_nodal_indices()
        self.adj_mat = self.generate_adj_mat(self.nx*self.ny)
        self.stiff_mat = self.generate_global_stiffness_mat()

        # calculate stuff about faces
        self.face_ids = np.arange(start=1, stop=self.num_faces+1, step=1)
        self.face_center_x, self.face_center_y = [], []
        for face_id in self.face_ids:
            x,y = self.get_face_center(face_id)
            self.face_center_x.append(x)
            self.face_center_y.append(y)

        self.face_center_x = np.array(self.face_center_x)
        self.face_center_y = np.array(self.face_center_y)

        return

    def solve(self):
        T,P = magicsolver(self.stiff_mat, self.nodal_temps, self.nodal_forces)

        self.nodal_temps = np.array(T)
        self.nodal_forces = np.array(P)

        flux = np.array([self.flux(face) for face in self.face_ids]).T
        gradient = np.array([self.gradient(face) for face in self.face_ids]).T

        self.face_flux_x = flux[0]
        self.face_flux_y = flux[1]
        self.face_gradient_x = -gradient[0]
        self.face_gradient_y = -gradient[1]

        pass

    def get_face_center(self, face):
        nodes = self.get_face_nodes(face)

        center_x = (self.nodal_coords_x[nodes[0]] + self.nodal_coords_x[nodes[1]] + self.nodal_coords_x[nodes[2]]) / 3.
        center_y = (self.nodal_coords_y[nodes[0]] + self.nodal_coords_y[nodes[1]] + self.nodal_coords_y[nodes[2]]) / 3.

        return center_x, center_y

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
        #assert(face > 0)
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
            # connectivity_mat = np.zeros((3,self.n))
            # connectivity_mat[0][nodes[0]] = 1
            # connectivity_mat[1][nodes[1]] = 1
            # connectivity_mat[2][nodes[2]] = 1


            # global_stiff_mat += connectivity_mat.T @ element_stiffness_mat.T @ connectivity_mat

            for i in range(3):
                global_stiff_mat[nodes[i]][nodes[0]] += element_stiffness_mat[i][0]
                global_stiff_mat[nodes[i]][nodes[1]] += element_stiffness_mat[i][1]
                global_stiff_mat[nodes[i]][nodes[2]] += element_stiffness_mat[i][2]

        return global_stiff_mat 

    def get_face_k(self, face):
        if face in self.mod_faces:
            k = self.k_mod
        else:
            k = self.k

        return k

        
    def generate_element_stiffness_mat(self, face): 
        assert(face > 0)
        assert(face < (self.num_faces + 1))

        k = self.get_face_k(face)

        nodes = self.get_face_nodes(face)

        # get absolute coordinates of every node
        x = [self.nodal_coords_x[node] for node in nodes]
        y = [self.nodal_coords_y[node] for node in nodes]

        # Zienkiewicz p. 120
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

        self.mod_faces = [-1] # remove modded faces

    def init_V1(self):
        # trapezoidal shape
        x_lengths = np.linspace(self.L, self.L/2, self.ny)
       
        nodal_coords_x = []
        for x_len in x_lengths:
            nodal_coords_x.append(np.linspace(0, x_len, self.nx))
        nodal_coords_x = np.array(nodal_coords_x).ravel()
        
        nodal_coords_y = []
        for y in np.linspace(0, self.L, self.ny):
            nodal_coords_y.append([y]*self.ny)
        nodal_coords_y = np.array(nodal_coords_y).ravel()

        self.nodal_coords_x = nodal_coords_x
        self.nodal_coords_y = nodal_coords_y
        return

    def init_V2(self):
        # generate square grid coords
        x = np.linspace(0.0, self.L, self.nx)
        y = np.linspace(0.0, self.L, self.ny)

        nodal_coords_x, nodal_coords_y = np.meshgrid(x,y)
        self.nodal_coords_x =  nodal_coords_x.ravel()
        self.nodal_coords_y =  nodal_coords_y.ravel()

        for i in range(len(self.nodal_coords_x)):
            B = 1/(2*self.L) * (self.L-self.nodal_coords_y[i])
            self.nodal_coords_x[i] =  self.nodal_coords_x[i] * (B/self.L*self.nodal_coords_x[i] - B + 1)


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

    def init_V4a(self):
        self.init_V0() # start with regular grid from V0

        self.mod_faces = self.file_params["all_elements"]
        self.k = float(self.file_params["k"])
        self.k_mod = float(self.file_params["k"])*float(self.file_params["c"])

    def init_V4b(self):
        self.init_V0() # start with regular grid from V0

        self.mod_faces = self.file_params["all_elements"]
        self.k = float(self.file_params["k"])
        self.k_mod = float(self.file_params["k"])/float(self.file_params["c"])

    def generate_adj_mat(self, nnodes):
        #returns adjacency matrix for regular mesh in Figure 1 (assignment)

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

    def perform_sanity_checks(self):
        # TODO: finish and integrate in V0
        # The temperature gradient must be constant in the entire model and equal to the overallgradient ∆T/∆y
        delta_T = self.nodal_temps.max() - self.nodal_temps.min()
        delta_y = self.L
        global_gradient = delta_T / delta_y
        grad = np.linalg.norm([self.face_gradient_x, self.face_gradient_y], axis=0)
        #res = np.allclose(, global_gradient)

        print("Local gradients constant and equal to the overall gradient?",
            np.allclose(grad_t[:, 1], global_gradient))
        print("Local fluxes constant and equal to the applied flux?", np.allclose(heat_flux[:, 1], q_0))
        print("Heat flux and temperature gradient yield k on each element?", np.allclose(-heat_flux[:, 1] / grad_t[:, 1], k))

    def flux(self, face):
        # calculate flux according to eq. (2) in the assignment
        gradient = self.gradient(face)
        k_ij = self.get_face_k(face)
        flux = k_ij * self.gradient(face)

        return flux

    def gradient(self, face):
        # calculate the gradient on a given face
        area = self.get_face_area(face)
        nodes = self.get_face_nodes(face)

        x = [self.nodal_coords_x[node] for node in nodes]
        y = [self.nodal_coords_y[node] for node in nodes]

        # I know this aprt is redundant but its not erally worth to make that efficient
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

        T = np.array([self.nodal_temps[nodes[0]], self.nodal_temps[nodes[1]], self.nodal_temps[nodes[2]]])
        M = np.array([
            b,
            c
            ])
        gradient = 1./(2*area) * M @ T # eq. (5) from assignment
        return gradient

    def plot_flux(self, filename='',title='Flux Plot'):
        elements=[]
        for a in range(0,self.num_faces+1):
            elements.append(self.get_face_nodes(a))
        triangulation = tri.Triangulation(self.nodal_coords_x, self.nodal_coords_y,elements)
        x,y=self.face_center_x, self.face_center_y
        plt.triplot(triangulation, '-k')
        plt.tricontourf(triangulation, self.nodal_temps)
        cbar = plt.colorbar()
        cbar.set_label('T [K]')
        plt.quiver(x,y, self.face_flux_x, self.face_flux_y)
        if filename != '':
            plt.savefig(filename)
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.title(title)
        plt.show()

    def plot_gradient(self, filename='',title='Gradient Plot'):
        elements=[]
        for a in range(0,self.num_faces+1):
            elements.append(self.get_face_nodes(a))
        triangulation = tri.Triangulation(self.nodal_coords_x, self.nodal_coords_y,elements)
        x,y=self.face_center_x, self.face_center_y
        plt.triplot(triangulation, '-k')
        plt.tricontourf(triangulation, self.nodal_temps)
        cbar = plt.colorbar()
        cbar.set_label('T [K]')
        plt.quiver(x,y, self.face_gradient_x, self.face_gradient_y)
        if filename != '':
            plt.savefig(filename)
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.title(title)
        plt.show()

    def plot_conservation_of_flow(self, filename='', title="no title",):
        first_face = 144
        x_vals = self.face_center_x[first_face:]

        # plot calculated values
        fluxes_out = np.linalg.norm([self.face_flux_x[first_face:], self.face_flux_y[first_face:]], axis=0)*10**(-6.)
        plt.plot(x_vals, fluxes_out, "x-", label="Calculated output flow")

        # plot input (it's just a straight line)
        fluxes_in = np.array([float(self.file_params["q_y_L"])]*len(fluxes_out))*10**(-6.)
        plt.plot(x_vals, fluxes_in, label="defined input flow")

        error = (fluxes_out - fluxes_in).sum() / fluxes_in.sum()

        plt.ylim([fluxes_in[0]*0.8, fluxes_in[0]*1.2])
        plt.title(f"{title}\nRelative Error = {error*100:.4}%")    # TODO: 2 kommastellen
        plt.grid()
        plt.xlabel("x [m]")
        plt.ylabel("flux [MW/m]")
        plt.legend()
        if filename != '':
            plt.savefig(filename)
        plt.show()
        pass



if __name__ == "__main__":
    mesh = Mesh("V4a")
    mesh.solve()
    mesh.plot_flux()
    mesh.plot_conservation_of_flow()

