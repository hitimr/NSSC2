import os 

WORKPSACE_DIR = os.path.dirname(os.path.realpath(__file__)) + "/"

COORD_X = 0
COORD_Y = 1
COORD_Z = 2

EPS = 0.0001

# enable 64 bit numbering for higher precision
#import jax.config
#jax.config.update("jax_enable_x64", True)

