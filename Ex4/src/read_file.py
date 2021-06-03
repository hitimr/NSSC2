#!/bin/python3
import argparse
import os, sys, inspect
import configparser

def read_from_file(filename):
    f = open(filename, "r")
    all_lines=f.readlines()
    i=0
    for line in all_lines:
        line_tmp=line.split(" = ")
        if line_tmp[0]=="groupnr":
            groupnr=line_tmp[1]
            #print(groupnr)
        if line_tmp[0]=="L":
            L=line_tmp[1]
            #print(L)
        if line_tmp[0]=="hz":
            hz=line_tmp[1]
            #print(hz)
        if line_tmp[0]=="k":
            k=line_tmp[1].replace('.','')
            #print(k)
        if line_tmp[0]=="c":
            c=line_tmp[1].replace('.','')
            #print(c)
        if line_tmp[0]=="q(y=L)":
            q_y_L=line_tmp[1].replace('.','')
            #print(q_y_L)
        if line_tmp[0]=="T(y=0)":
            T_y_0=line_tmp[1].replace('.','')
            #print(T_y_0)
        if line_tmp[0]=="elements_to_be_modified":
            all_elements=[]
            for a in range(1,8):
                elements=all_lines[i+a].split("-")
                num_elements=int(elements[1])-int(elements[0])
                for b in range(0,num_elements+1):
                    all_elements.append(int(elements[0])+b)
        i=i+1
    #print(all_elements)
    return L,hz,k,c,q_y_L,T_y_0,all_elements





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to generate trajectory')
    parser.add_argument('inputfile', type=str,help='inputfile', default="inputfile_group_5.txt", nargs='?', const=1)

    args = parser.parse_args()
    L,hz,k,c,q_y_L,T_y_0,all_elements=read_from_file(args.inputfile)
    print("L = "+L+"\nhz = "+hz+"\nk = "+k+"\nc = "+c+"\nq_y_L = "+q_y_L+"\nT_y_0 = "+T_y_0+"\nelements_to_be_modified = "+str(all_elements))
