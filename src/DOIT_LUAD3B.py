#!/usr/bin/env python3

import sys
import subprocess

def main():

    a0 = sys.argv[0]
    pp = pathlib.Path(a0)
    SRCDIR = str(pp.parent)

    WORKDIR = "."

    command = SRCDIR + "/DOIT_GetDATA_LUAD3B.py" 
    print(command)
    p = subprocess.run(command, shell=True, check=True)

    command = SRCDIR + "/DOIT_Visium_LUAD3B.py" 
    print(command)
    p = subprocess.run(command, shell=True, check=True)

    command = SRCDIR + "/DOIT_PhenoCycler_LUAD3B.py" 
    print(command)
    p = subprocess.run(command, shell=True, check=True)

    command = SRCDIR + "/DOIT_Normalize_LUAD3B.py" 
    print(command)
    p = subprocess.run(command, shell=True, check=True)

    command = SRCDIR + "/DOIT_Align_LUAD3B.py" 
    print(command)
    p = subprocess.run(command, shell=True, check=True)

    command = SRCDIR + "/DOIT_Eval_LUAD3B.py" 
    print(command)
    p = subprocess.run(command, shell=True, check=True)
    
    command = SRCDIR + "/DOIT_Clustering_LUAD3B.py" 
    print(command)
    p = subprocess.run(command, shell=True, check=True)
    
    command = SRCDIR + "/DOIT_GetROI_LUAD3B.py" 
    print(command)
    p = subprocess.run(command, shell=True, check=True)
    
#####

if __name__ == "__main__":
    main()

    
