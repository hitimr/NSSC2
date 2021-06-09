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
        A = self.mesh.stiff_mat
        b = self.mesh.nodal_forces

        x = np.linalg.solve(A, b)
        print(x)



if __name__ == "__main__":
    mesh = Mesh("V1")

    kernel = FEM_Kernel()
    kernel.set_mesh(mesh)
    kernel.solve()
