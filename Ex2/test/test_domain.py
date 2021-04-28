#include parent folder
import os, sys, inspect

currentdir = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import numpy as np
import pytest

from src.domain import *


def test_Epot():
    pos = np.array([
        [1,1,1],
        [2,2,2],
        [3,3,3]
    ])
    
    E = VLJ(pos[0], pos[1]) + VLJ(pos[0], pos[2]) + VLJ(pos[1], pos[2])

    assert(Epot(pos) == E)




if __name__ == "__main__":
    test_Epot()