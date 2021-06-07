#include parent folder
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import pytest
import numpy as np
from src.mesh import *



def test_generate_adj_mat():
    mesh = Mesh('debug')
    mat = mesh.generate_adj_mat(10*10)


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
    
    assert(mesh.get_face_nodes(1).all() == np.array([0,1,10]).all())
    assert(mesh.get_face_nodes(2).all() == np.array([1,10,11]).all())
    assert(mesh.get_face_nodes(162).all() == np.array([89,98,99]).all())

def test_generate_adj_mat():
    mesh = Mesh('V1')

    mat = mesh.adj_mat
    # TODO unit test

if __name__ == "__main__":
    test_get_face_nodes()