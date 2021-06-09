#include parent folder
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import pytest
import numpy as np
from src.mesh import *


def test_init():
    mesh = Mesh('debug')
    assert(mesh.num_faces == 162)
    assert(len(mesh.face_ids) == 162)
    assert(mesh.face_ids[0] == 1)
    assert(mesh.face_ids[1] == 2)




def test_generate_nodal_indices():
    mesh = Mesh('debug')
    
    mat = mesh.index_mat

    assert(mat[0][0] == 0)
    assert(mat[1][0] == 1)
    assert(mat[0][1] == 10)
    assert(mat[9][9] == 99)
    assert(mat[1][1] == 11)

def test_get_face_nodes():
    mesh = Mesh('debug')
    a = mesh.get_face_nodes(0)
    pass
    return # TODO: fix tests if there is time 
    assert(np.allclose(mesh.get_face_nodes(1).sort(), [0,1,10]))
    assert(np.allclose(mesh.get_face_nodes(2).sort(), [1,10,11]))
    assert(np.allclose(mesh.get_face_nodes(5).sort(), [2,3,12]))
    assert(np.allclose(mesh.get_face_nodes(6).sort(), [3,12,13]))
    assert(np.allclose(mesh.get_face_nodes(7).sort(), [3,4,13]))
    assert(np.allclose(mesh.get_face_nodes(8).sort(), [4,13,14]))
    assert(np.allclose(mesh.get_face_nodes(9).sort(), [4,5,14]))
    assert(np.allclose(mesh.get_face_nodes(10).sort(), [5,14,15]))
    assert(np.allclose(mesh.get_face_nodes(11).sort(), [5,6,15]))
    assert(np.allclose(mesh.get_face_nodes(12).sort(), [6,15,16]))
    assert(np.allclose(mesh.get_face_nodes(13).sort(), [6,7,16]))
    assert(np.allclose(mesh.get_face_nodes(16).sort(), [8, 17,18]))
    assert(np.allclose(mesh.get_face_nodes(17).sort(), [8,9,18]))
    assert(np.allclose(mesh.get_face_nodes(18).sort(), [9,18,19]))
    assert(np.allclose(mesh.get_face_nodes(29).sort(), [15,16,25]))


def test_adj_mat():
    mesh = Mesh('debug')
    mat = mesh.adj_mat

    #face 1
    assert(mat[0][1] == 1)
    assert(mat[0][11] == 0)
    assert(mat[0][10] == 1)
    assert(mat[1][10] == 1)

    #face 2
    assert(mat[1][11] == 1)
    assert(mat[1][10] == 1)
    assert(mat[10][11] == 1)

    # face 5
    assert(mat[2][3] == 1)
    assert(mat[3][2] == 1)
    assert(mat[2][12] == 1)
    assert(mat[12][2] == 1)
    assert(mat[3][12] == 1)
    assert(mat[12][3] == 1)

    #face 6
    assert(mat[3][4] == 1)
    assert(mat[3][13] == 1)
    assert(mat[4][13] == 1)

    # node 11
    assert(mat[11][2] == 1)
    assert(mat[11][10] == 1)
    assert(mat[11][21] == 1)
    assert(mat[11][12] == 1)
    assert(mat[11][2] == 1)


if __name__ == "__main__":
    test_get_face_nodes()
    test_generate_adj_mat