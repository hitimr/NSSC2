#include parent folder
import os, sys, inspect

currentdir = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import numpy as np
import pytest
from src.misc import *
from src.domain import *

def test_to3d():
    res = np.array([
        [1,2,3],
        [4,5,6],
        [7,8,9]
    ])

    assert np.isclose( to3D([1,2,3,4,5,6,7,8,9]), res).all
    assert np.isclose( to3D(np.array([1,2,3,4,5,6,7,8,9])), res).all

    for i in range(1000):
        res = np.random.rand(i,3)
        assert np.isclose( to3D(res.ravel()), res).all
        assert np.isclose( to3D(res), res).all 

    


if __name__ == "__main__":
    test_to3d()