#include parent folder
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 


import pytest

from config import *

def test_task2():
    # change dir to task2
    os.system(f"cd {TASK2_DIR}")
    os.system(f"rm -rf out")




# for debugging
if __name__ == "__main__":
    test_task2()