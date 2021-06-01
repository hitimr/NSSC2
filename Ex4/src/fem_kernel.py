import numpy as np

from mesh import Mesh


class FEM_Kernel():
    def __init__(self):
        pass

    def set_mesh(self, mesh):
        assert isinstance(mesh, Mesh)

        self.mesh = mesh

    def solve(self):
        # solve Ax = b
        A = np.linalg.inv(self.mesh.adj_mat)
        b = self.mesh.nodal_temps

        x = np.linalg.solve(A, b)
        print(x)



if __name__ == "__main__":
    mesh = Mesh("V0", 4)

    kernel = FEM_Kernel()
    kernel.set_mesh(mesh)
    kernel.solve()
