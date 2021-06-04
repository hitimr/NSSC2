#include parent folder
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import pytest
from src.magicsolver import *



def test_magicseolver():
    assert(True)

    with pytest.raises(Exception):
        assert(False)




if __name__ == "__main__":
    test_magicseolver()