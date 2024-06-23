#!/usr/bin/env python3

import sys
import subprocess


def main():

    a0 = sys.argv[0]
    pp = pathlib.Path(a0)
    SRCDIR = str(pp.parent)

    WORKDIR = "."

    sid = "FFPE_LUAD_3_B"    
    filename_qptiff = "./dat/PhenoCycler_LUAD_3_B/Experiment_LUAD_3_B.qptiff"
    
    # Parse PhenoCycler
    command = SRCDIR + "/Parse_PhenoCycler_QP.py" + " " + filename_qptiff  + " " + sid
    print(command)

    p = subprocess.run(command, shell=True, check=True)

#####

if __name__ == "__main__":
    main()

    
