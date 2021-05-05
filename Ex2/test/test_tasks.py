#include parent folder
import os, sys, inspect

currentdir = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import numpy as np
import pytest

from src.task2 import *
from src.misc import *


def test_task2():
    fileName = DIR_OUT + "unitTestOut_task2.txt"
    task2(200, 6.3, 1, fileName)
    assert os.path.exists(fileName)
