#!/usr/bin/env python3

import sys
import subprocess
import pathlib

def main():

    if len(sys.argv) < 3:
        print("sys.argv")
        sys.exit()

    sid = sys.argv[1]
    filename_dat = sys.argv[2]

    WORKDIR = "."

    filename_dat = "https://gds.k.u-tokyo.ac.jp/pub/LUAD/LUAD3B_dat.tgz"

    pp = pathlib.Path(filename_dat)
    filename_tgz = str(pp.name)

    # wget
    command = "curl -OL" + " " + filename_dat
    print(command)
    p = subprocess.run(command, shell=True, check=True)

    # tar+gzip
    command = "tar -zxvf" + " " + filename_tgz
    print(command)
    p = subprocess.run(command, shell=True, check=True)

    

#####

if __name__ == "__main__":
    main()
