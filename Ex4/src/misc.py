#!/bin/python3
import os
import platform
import inspect
import argparse
import configparser
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from matplotlib.collections import LineCollection
import matplotlib.tri as tri



# System specific directory separator
if platform.system() in ["Darwin", "Linux"]:
    SEP = "/"
else: # platform.system() == "Windows":
    SEP = "\\"

DIR_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))) + SEP
DIR_OUT = DIR_ROOT + "out" + SEP
DIR_SRC = DIR_ROOT + "src" + SEP



UNKNOWN_NODAL_TEMP = None
UNKNOWN_NODAL_FORCE = None


def read_from_file(filename):
    f = open(filename, "r")
    all_lines=f.readlines()
    input_dict={}
    i=0
    for line in all_lines:
        line_tmp=line.split(" = ")
        if line_tmp[0]=="groupnr":
            groupnr=line_tmp[1]
            input_dict.update({"groupnr": line_tmp[1].replace('\n','')})
        if line_tmp[0]=="L":
            L=line_tmp[1]
            input_dict.update({"L": line_tmp[1].replace('\n','')})
        if line_tmp[0]=="hz":
            hz=line_tmp[1]
            input_dict.update({"hz": line_tmp[1].replace('\n','')})
        if line_tmp[0]=="k":
            k=line_tmp[1].replace('.','')
            input_dict.update({"k": line_tmp[1].replace('.\n','')})
        if line_tmp[0]=="c":
            c=line_tmp[1].replace('.','')
            input_dict.update({"c": line_tmp[1].replace('.\n','')})
        if line_tmp[0]=="q(y=L)":
            q_y_L=line_tmp[1].replace('.','')
            input_dict.update({"q_y_L": line_tmp[1].replace('.\n','')})
        if line_tmp[0]=="T(y=0)":
            T_y_0=line_tmp[1].replace('.','')
            input_dict.update({"T_y_0": line_tmp[1].replace('.\n','')})
        if line_tmp[0]=="elements_to_be_modified":
            all_elements=[]
            for a in range(1,8):
                elements=all_lines[i+a].split("-")
                num_elements=int(elements[1])-int(elements[0])
                for b in range(0,num_elements+1):
                    all_elements.append(int(elements[0])+b)
            input_dict.update({"all_elements": all_elements})
        i=i+1
    #return L,hz,k,c,q_y_L,T_y_0,all_elements
    return input_dict




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to generate trajectory')
    parser.add_argument('inputfile', type=str,help='inputfile', default="inputfile_group_5.txt", nargs='?', const=1)

    args = parser.parse_args()
    #L,hz,k,c,q_y_L,T_y_0,all_elements=read_from_file(args.inputfile)
    input_dict=read_from_file(args.inputfile)
    print(input_dict.get("k"))
