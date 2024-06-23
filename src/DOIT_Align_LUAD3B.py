#!/usr/bin/env python3

import sys
import subprocess

def main():

    a0 = sys.argv[0]
    pp = pathlib.Path(a0)
    SRCDIR = str(pp.parent)

    sid = "FFPE_LUAD_3_B"    

    WORKDIR = "."
    RESULTDIR = WORKDIR + "/" + sid    
    
    PREFIX_v = sid + "-" + "Visium"    
    MASKDIR_v = RESULTDIR + "/" + PREFIX_v + "-" + "MASK"
    MASKDIR_v_DAT = MASKDIR_v + "/" + "dat"
    filename_visium_mask = MASKDIR_v_DAT + "/" + PREFIX_v + "_tr_gray_without_circle_nega_mask.png"        

    PREFIX_p = sid + "-" + "PhenoCycler"
    MASKDIR_p = RESULTDIR + "/" + PREFIX_p + "-" + "MASK"
    MASKDIR_p_DAT = MASKDIR_p + "/" + "dat"
    MASKDIR_p_DAT_MULTI = MASKDIR_p_DAT + "/" + "MULTI"

    filename_phenocycler_mask_01 = MASKDIR_p_DAT_MULTI + "/" + PREFIX_p + "-" + "MULTI_gray_nega_mask_01.png"

    # Align

    command = SRCDIR + "/Align_Visium-PhenoCycler.py" + " " + filename_visium_mask + " " + filename_phenocycler_mask_01 + " " + sid
    print(command)
    p = subprocess.run(command, shell=True, check=True)
    

#####

if __name__ == "__main__":
    main()

    
