import numpy as np


class Mesh:
    variation = "noVariation"
    size = -1
    adj_mat = np.array

    def __init__(self, variation, size):
        assert(type(variation) == str)
        assert(type(size) == int and size > 2)

        if variation == "V0": self.init_V0()
        else: raise ValueError(f"Unknown argumen '{variation}' for variation")

        self.variation = variation

    def init_V0(self):
        """Initialize using the simplest case:

        2 - 3
        | \ |
        0 - 1

        adjacency matrix is hardcoded
            0   1   2   3   
        0       1   1   
        1   1           1
        2   1           1
        3       1   1

        """
        self.adj_mat = np.array([
            [0,1,1,0],
            [1,0,0,1],
            [1,0,0,1],
            [0,1,1,0],
        ])
        print(self.adj_mat)
        pass



if __name__ == "__main__":
    mesh = Mesh("V0", 4)