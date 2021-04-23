#include parent folder 
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) 
parentdir = os.path.dirname(currentdir) 
sys.path.insert(0,parentdir)  

import pytest
import pathlib
from config import *


def test_dirs():
    assert os.path.exists(WORKPSACE_DIR)