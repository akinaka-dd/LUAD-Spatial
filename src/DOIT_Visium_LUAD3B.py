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
    filename_hires_png = dirname_spaceranger + "/spatial/tissue_hires_image.png"
    filename_fulres_png = WORKDIR + "/dat/Visium_FFPE_LUAD_3_B/FFPE_LUAD_3_B_x10.tif"
    
    # Parse Visium
    command = SRCDIR + "/Parse_Visium_HE.py" + " " + dirname_spaceranger  + " " + filename_fulres_png + " " + sid 
    print(command)
    p = subprocess.run(command, shell=True, check=True)

#####

if __name__ == "__main__":
    main()

    
