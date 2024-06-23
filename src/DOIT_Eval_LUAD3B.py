#!/usr/bin/env python3

import sys
import subprocess

def main():
    a0 = sys.argv[0]
    pp = pathlib.Path(a0)
    SRCDIR = str(pp.parent)

    WORKDIR = "."

    sid = "FFPE_LUAD_3_B"

    filename_seg = WORKDIR + "/dat/PhenoCycler_LUAD_3_B/segmentation_LUAD_3_B.tsv"

    # Eval
    command = SRCDIR + "/Eval_Visium-PhenoCycler.py" + " " + sid + " " + filename_seg
    print(command)
    p = subprocess.run(command, shell=True, check=True)
    

#####

if __name__ == "__main__":
    main()

    
