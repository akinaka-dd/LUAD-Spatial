#!/usr/bin/env python3

import sys
import subprocess

def main():

    a0 = sys.argv[0]
    pp = pathlib.Path(a0)
    SRCDIR = str(pp.parent)

    WORKDIR = "."

    sid = "FFPE_LUAD_3_B"
    dirname_spaceranger = WORKDIR + "/dat/Visium_FFPE_LUAD_3_B"
    
    # Normalize
    command = SRCDIR + "/Normalize_Visium.py" + " " + dirname_spaceranger + " " + sid
    print(command)
    p = subprocess.run(command, shell=True, check=True)

    
#####

if __name__ == "__main__":
    main()

    
