#!/usr/bin/env python3

import sys
import subprocess
import pathlib

def main():

    a0 = sys.argv[0]
    pp = pathlib.Path(a0)
    SRCDIR = str(pp.parent)

    WORKDIR = "."

    sid = "FFPE_LUAD_3_B"

    filename_dat = "https://gds.k.u-tokyo.ac.jp/pub/LUAD/LUAD3B_dat.tgz"    

    # GetDATA
    command = SRCDIR + "/GetDATA.py" + " " + sid + " " + filename_dat
    print(command)
    p = subprocess.run(command, shell=True, check=True)
    

#####

if __name__ == "__main__":
    main()
