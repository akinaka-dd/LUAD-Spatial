#!/usr/bin/env python3

import sys
import subprocess

def main():

    a0 = sys.argv[0]
    pp = pathlib.Path(a0)
    SRCDIR = str(pp.parent)

    WORKDIR = "."

    sid = "FFPE_LUAD_3_B"    

    # Clustering
    command = SRCDIR + "/Clustering_Visium.py" + " " + sid
    print(command)
    p = subprocess.run(command, shell=True, check=True)
    

#####

if __name__ == "__main__":
    main()

    
