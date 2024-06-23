#!/usr/bin/env python3

import sys
import re
import pathlib
import shutil
import numpy as np
import math
import cv2

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from util.image import *
from util.html import *

binary_th = 10

def main():

    color_and = lightgray

    if len(sys.argv) < 2:
        print("sys.argv")
        sys.exit()

    filename_visium_mask = sys.argv[1]
    filename_phenocycler_mask_01 = sys.argv[2]
    
    sid = sys.argv[3]

    print(filename_visium_mask)
    print(filename_phenocycler_mask_01)
    print("SID", sid)

    # directories
    
    WORKDIR = "."
    RESULTDIR = WORKDIR + "/" + sid

    PREFIX = sid    
    
    # Visium
    PREFIX_v = sid + "-" + "Visium"
    MASKDIR_v = RESULTDIR + "/" + PREFIX_v + "-" + "MASK"
    MASKDIR_v_DAT = MASKDIR_v + "/" + "dat"    
    print("MASKDIR_v", MASKDIR_v)
    
    # PhenoCycler
    PREFIX_p = sid + "-" + "PhenoCycler"
    MASKDIR_p = RESULTDIR + "/" + PREFIX_p + "-" + "MASK"
    MASKDIR_p_DAT = MASKDIR_p + "/" + "dat"
    MASKDIR_p_DAT_DAPI = MASKDIR_p_DAT + "/" + "DAPI"
    MASKDIR_p_DAT_MULTI = MASKDIR_p_DAT + "/" + "MULTI"
    print("MASKDIR_p_DAT_DAPI", MASKDIR_p_DAT_DAPI)
    print("MASKDIR_p_DAT_MULTI", MASKDIR_p_DAT_MULTI)

    # ALIGNDIR 
    ALIGNDIR = RESULTDIR + "/" + PREFIX + "-" + "ALIGN"
    pp = pathlib.Path(ALIGNDIR)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)
    ALIGNDIR_DAT = ALIGNDIR + "/" + "dat"        
    pp = pathlib.Path(ALIGNDIR_DAT)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)

    ALIGNDIR_DAT_PPMNG = ALIGNDIR_DAT + "/" + "ppmng"
    pp = pathlib.Path(ALIGNDIR_DAT_PPMNG)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)

    # Masks
    # Visium
    m = re.match("(.+?)\.png", filename_visium_mask)
    if m:
        filename_body = m.group(1)
        filename_visium_mask_wh = filename_body + "_with_hole.png"
    else:
        print("Error:VISIUM MASK")
        sys.exit()
        
    # PhenoCycler
    m = re.match("(.+?)_01\.png", filename_phenocycler_mask_01)
    if m:
        filename_body = m.group(1)
        filename_phenocycler_mask_wh_01 = filename_body + "_with_hole_01.png"
    else:
        print("Error:PHENOCYCLER MASK")
        sys.exit()

    # Original images (nega)
    
    # without hole
    img_visium_mask_00 = cv2.imread(filename_visium_mask, cv2.IMREAD_GRAYSCALE)
    img_phenocycler_mask_01 = cv2.imread(filename_phenocycler_mask_01, cv2.IMREAD_GRAYSCALE)
    
    # with hole
    img_visium_mask_wh_00 = cv2.imread(filename_visium_mask_wh, cv2.IMREAD_GRAYSCALE)
    img_phenocycler_mask_wh_01 = cv2.imread(filename_phenocycler_mask_wh_01, cv2.IMREAD_GRAYSCALE)

    # Expanded images ... add margins for transformation
    dm0 = 0
    dm1 = 0
    if sid == "FFPE_LUAD_2_C":
        dm0 = 300
        dm1 = 0
    elif sid == "FFPE_LUAD_4_C":
        dm0 = 300
        dm1 = 100
    else: # "FFPE_LUAD_3_B"
        dm0 = 0
        dm1 = 0
        
    vh, vw = img_visium_mask_00.shape
    img_visium_mask = np.zeros((vh+dm0+dm1, vw+dm0+dm1), dtype=np.uint8)
    img_visium_mask_wh = np.zeros((vh+dm0+dm1, vw+dm0+dm1), dtype=np.uint8)

    if dm1 < 0:
        vh += dm1
        vw += dm1

    img_visium_mask[dm0:dm0+vh, dm0:dm0+vw] = img_visium_mask_00[:vh, :vw]
    img_visium_mask_wh[dm0:dm0+vh, dm0:dm0+vw] = img_visium_mask_wh_00[:vh, :vw]

    ph, pw = img_phenocycler_mask_01.shape
    img_phenocycler_mask = np.zeros((ph+dm0+dm1, pw+dm0+dm1), dtype=np.uint8)
    img_phenocycler_mask_wh = np.zeros((ph+dm0+dm1, pw+dm0+dm1), dtype=np.uint8)
    
    if dm1 < 0:
        ph += dm1
        pw += dm1
    
    img_phenocycler_mask[dm0:dm0+ph, dm0:dm0+pw] = img_phenocycler_mask_01[:ph, :pw]
    img_phenocycler_mask_wh[dm0:dm0+ph, dm0:dm0+pw] = img_phenocycler_mask_wh_01[:ph, :pw]
    
    #####

    filename_html = ALIGNDIR + "/" + PREFIX + "-Visium" + "-" "PhenoCycler" + "-ALIGN" + ".html"
    print(filename_html)
    fp_html = open(filename_html, "wt")

    ww = "200px"
    ww_mp4 = "300px"
    ww_mp4_s = "200px"
    w_s = 200 # width (pixels) of the small mp4    

    print_html_header(fp_html)    

    print("<h1>", file=fp_html)
    print("Alignment of Visium and PhenoCycler data" + " &mdash; " + sid, sid, file=fp_html)
    print("</h1>", file=fp_html)

    ######################

    print("<h2>", file=fp_html)
    print("Image alignment", file=fp_html)
    print("</h2>", file=fp_html)

    filename_txt = ALIGNDIR_DAT + "/" + sid + "_matrix" + ".txt"
    filename_txt_rel = "dat" + "/" + sid + "_matrix" + ".txt"    

    print("<a href=\"" + filename_txt_rel + "\">" + "Affine trasformation matrices" + "</a>", file=fp_html)
    print("<br><br>", file=fp_html)

    # MP4
    
    filename_mp4_rel = "dat/ppmng" + "/" + sid + "_transition" + ".mp4"
    filename_mp4_wh_rel = "dat/ppmng" + "/" + sid + "_transition_wh" + ".mp4"

    filename_mp4_s_rel = "dat/ppmng" + "/" + sid + "_transition_{:d}".format(w_s) + ".mp4"
    filename_mp4_wh_s_rel = "dat/ppmng" + "/" + sid + "_transition_wh_{:d}".format(w_s) + ".mp4"

    print("<h3>", file=fp_html)
    print("Alignment animation", file=fp_html)
    print("</h3>", file=fp_html)


    print("<video src=" + "\"" + filename_mp4_rel + "\"" + " width=\"" + ww_mp4 + "\"" + " controls playsinline" + "></video>", file=fp_html)
    print("&emsp;", file=fp_html)    
    print("<video src=" + "\"" + filename_mp4_wh_rel + "\"" + " width=\"" + ww_mp4 + "\"" + " controls playsinline" + "></video>", file=fp_html)

    print("&emsp;", file=fp_html)
    print("&emsp;", file=fp_html)        

    print("<video src=" + "\"" + filename_mp4_s_rel + "\"" + " width=\"" + ww_mp4_s + "\"" + " controls playsinline" + "></video>", file=fp_html)
    print("&emsp;", file=fp_html)    
    print("<video src=" + "\"" + filename_mp4_wh_s_rel + "\"" + " width=\"" + ww_mp4_s + "\"" + " controls playsinline" + "></video>", file=fp_html)
    
    #####

    h_v, w_v = img_visium_mask.shape
    h_p, w_p = img_phenocycler_mask.shape
    h_max = max(h_v, h_p)
    w_max = max(w_v, w_p)

    # INIT

    # phenocycler (codex)
    contours_cod, hierachy_cod = cv2.findContours(img_phenocycler_mask, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    contour_cod_max_area = max(contours_cod, key=lambda x: cv2.contourArea(x))
    contour_cod_max = cv2.contourArea(contour_cod_max_area)
    img_phenocycler_mask_resized = np.zeros((h_max, w_max), dtype=np.uint8)
    img_phenocycler_mask_resized[:h_p, :w_p] = img_phenocycler_mask[:h_p, :w_p]
    
    # visium
    contours_vis, hierachy_vis = cv2.findContours(img_visium_mask, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    contour_vis_max_area = max(contours_vis, key=lambda x: cv2.contourArea(x))
    contour_vis_max = cv2.contourArea(contour_vis_max_area)
    img_visium_mask_resized = np.zeros((h_max, w_max), dtype=np.uint8)
    img_visium_mask_resized[:h_v, :w_v] = img_visium_mask[:h_v, :w_v]
    
    img_and = cv2.bitwise_and(img_phenocycler_mask_resized, img_visium_mask_resized)    
    
    img_init = np.ones((h_max, w_max)) * 255
    img_init = cv2.cvtColor(img_init.astype(np.uint8), cv2.COLOR_GRAY2BGR)
    img_init[np.where(img_visium_mask_resized == 255)] = orangered    
    img_init[np.where(img_phenocycler_mask_resized == 255)] = royalblue
    img_init[np.where(img_and == 255)] = color_and 
    img_init = bg_transp(img_init)

    filename_init = ALIGNDIR_DAT + "/" + sid + "_init" + ".png"
    filename_init_rel = "dat" + "/" + sid + "_init" + ".png"
    cv2.imwrite(filename_init, img_init)

    # with hole

    img_phenocycler_mask_resized_wh = img_phenocycler_mask_resized.copy()
    img_phenocycler_mask_resized_wh[np.where(img_phenocycler_mask_wh != 255)] = 0

    img_visium_mask_resized_wh = img_visium_mask_resized.copy()
    img_visium_mask_resized_wh[np.where(img_visium_mask_wh != 255)] = 0

    img_and_wh = cv2.bitwise_and(img_phenocycler_mask_resized_wh, img_visium_mask_resized_wh)
    
    img_init_wh = np.ones((h_max, w_max)) * 255
    img_init_wh = cv2.cvtColor(img_init_wh.astype(np.uint8), cv2.COLOR_GRAY2BGR)
    img_init_wh[np.where(img_visium_mask_resized_wh == 255)] = orangered    
    img_init_wh[np.where(img_phenocycler_mask_resized_wh == 255)] = royalblue
    img_init_wh[np.where(img_and_wh == 255)] = color_and 
    img_init_wh = bg_transp(img_init_wh)

    filename_init_wh = ALIGNDIR_DAT + "/" + sid + "_init_wh" + ".png"
    filename_init_wh_rel = "dat" + "/" + sid + "_init_wh" + ".png"
    cv2.imwrite(filename_init_wh, img_init_wh)

    count_and, count_or, count_vis, count_cod = COUNT(img_visium_mask_resized_wh, img_phenocycler_mask_resized_wh)
    cr = count_and / count_or
    xor = (count_or - count_and)

    print("<h3>", file=fp_html)
    print("Original images", file=fp_html)
    print("</h3>", file=fp_html)

    print("IoU = {0:.4f} = {1:d} / {2:d} = {1:d} / ({3:d} + {4:d} - {1:d}) = |V &cap; P| / (|V| + |P| - |V &cap; P|) <br>"\
          .format(cr, count_and, count_or, count_vis, count_cod), file=fp_html)    
    print("XOR = {0:d} = {1:d} - {2:d} = |V &cup; P| - |V &cap; P| <br>".format(xor, count_or, count_and), file=fp_html)

    print("<img src=" + "\"" + filename_init_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    print("<img src=" + "\"" + filename_init_wh_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)    

    fp_html.flush()
    
    #####

    ## ZOOM
    
    ratio = math.sqrt(contour_vis_max / contour_cod_max)
    print("RATIO", ratio)

    img_phenocycler_mask_resized_zoomed, img_phenocycler_mask_resized_zoomed_wh = \
        ZOOM(img_phenocycler_mask_resized, img_phenocycler_mask_resized_wh, h_max, w_max, ratio)

    # without hole
    img_and = cv2.bitwise_and(img_phenocycler_mask_resized_zoomed, img_visium_mask_resized)

    img_zoom = np.ones((h_max, w_max)) * 255
    img_zoom = cv2.cvtColor(img_zoom.astype(np.uint8), cv2.COLOR_GRAY2BGR)
    img_zoom[np.where(img_visium_mask_resized == 255)] = orangered    
    img_zoom[np.where(img_phenocycler_mask_resized_zoomed == 255)] = royalblue
    img_zoom[np.where(img_and == 255)] = color_and 
    img_zoom = bg_transp(img_zoom)
    
    filename_zoom = ALIGNDIR_DAT + "/" + sid + "_zoom" + ".png"
    filename_zoom_rel = "dat" + "/" + sid + "_zoom" + ".png"
    cv2.imwrite(filename_zoom, img_zoom)

    # with hole
    img_and_wh = cv2.bitwise_and(img_phenocycler_mask_resized_zoomed_wh, img_visium_mask_resized_wh)

    img_zoom_wh = np.ones((h_max, w_max)) * 255
    img_zoom_wh = cv2.cvtColor(img_zoom_wh.astype(np.uint8), cv2.COLOR_GRAY2BGR)
    img_zoom_wh[np.where(img_visium_mask_resized_wh == 255)] = orangered    
    img_zoom_wh[np.where(img_phenocycler_mask_resized_zoomed_wh == 255)] = royalblue
    img_zoom_wh[np.where(img_and_wh == 255)] = color_and 
    img_zoom_wh = bg_transp(img_zoom_wh)

    filename_zoom_wh = ALIGNDIR_DAT + "/" + sid + "_zoom_wh" + ".png"
    filename_zoom_wh_rel = "dat" + "/" + sid + "_zoom_wh" + ".png"
    cv2.imwrite(filename_zoom_wh, img_zoom_wh)

    count_and, count_or, count_vis, count_cod = COUNT(img_visium_mask_resized_wh, img_phenocycler_mask_resized_zoomed_wh)
    cr = count_and / count_or
    xor = (count_or - count_and)

    print("<h3>", file=fp_html)
    print("Zooming in/out", file=fp_html)
    print("</h3>", file=fp_html)

    print("PhenoCycler &times; {:.4f}<br>".format(ratio), file=fp_html)
    print("IoU = {0:.4f} = {1:d} / {2:d} = {1:d} / ({3:d} + {4:d} - {1:d}) = |V &cap; P| / (|V| + |P| - |V &cap; P|) <br>"\
          .format(cr, count_and, count_or, count_vis, count_cod), file=fp_html)    
    print("XOR = {0:d} = {1:d} - {2:d} = |V &cup; P| - |V &cap; P| <br>".format(xor, count_or, count_and), file=fp_html)
    
    print("<img src=" + "\"" + filename_zoom_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    print("<img src=" + "\"" + filename_zoom_wh_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    fp_html.flush()

    
    ## SHIFT/TRANSLATE

    contours_cod_resized, hierachy_cod_resized = cv2.findContours(img_phenocycler_mask_resized_zoomed, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)    
    contour_cod_resized_max_area = max(contours_cod_resized, key=lambda x: cv2.contourArea(x))
    contour_cod_resized_max = cv2.contourArea(contour_cod_resized_max_area)

    mom_vis = cv2.moments(contour_vis_max_area)
    mom_cod = cv2.moments(contour_cod_resized_max_area)

    xc_vis = int(mom_vis["m10"]/mom_vis["m00"])
    yc_vis = int(mom_vis["m01"]/mom_vis["m00"])    

    xc_cod = int(mom_cod["m10"]/mom_cod["m00"])
    yc_cod = int(mom_cod["m01"]/mom_cod["m00"])

    dx = xc_vis - xc_cod
    dy = yc_vis - yc_cod

    img_phenocycler_mask_resized_zoomed_shifted, img_phenocycler_mask_resized_zoomed_shifted_wh = \
        SHIFT(img_phenocycler_mask_resized_zoomed, img_phenocycler_mask_resized_zoomed_wh, dx, dy)

    img_and = cv2.bitwise_and(img_phenocycler_mask_resized_zoomed_shifted, img_visium_mask_resized)

    img_shift = np.ones((h_max, w_max)) * 255    
    img_shift = cv2.cvtColor(img_shift.astype(np.uint8), cv2.COLOR_GRAY2BGR)
    img_shift[np.where(img_visium_mask_resized == 255)] = orangered
    img_shift[np.where(img_phenocycler_mask_resized_zoomed_shifted == 255)] = royalblue    
    img_shift[np.where(img_and == 255)] = color_and 
    img_shift = bg_transp(img_shift)
    
    filename_shift = ALIGNDIR_DAT + "/" + sid + "_shift" + ".png"
    filename_shift_rel = "dat" + "/" + sid + "_shift" + ".png"
    cv2.imwrite(filename_shift, img_shift)

    # with hole
    
    img_and_wh = cv2.bitwise_and(img_phenocycler_mask_resized_zoomed_shifted_wh, img_visium_mask_resized_wh)

    img_shift_wh = np.ones((h_max, w_max)) * 255    
    img_shift_wh = cv2.cvtColor(img_shift_wh.astype(np.uint8), cv2.COLOR_GRAY2BGR)
    img_shift_wh[np.where(img_visium_mask_resized_wh == 255)] = orangered    
    img_shift_wh[np.where(img_phenocycler_mask_resized_zoomed_shifted_wh == 255)] = royalblue
    img_shift_wh[np.where(img_and_wh == 255)] = color_and 
    img_shift_wh = bg_transp(img_shift_wh)
    
    filename_shift_wh = ALIGNDIR_DAT + "/" + sid + "_shift_wh" + ".png"
    filename_shift_wh_rel = "dat" + "/" + sid + "_shift_wh" + ".png"
    cv2.imwrite(filename_shift_wh, img_shift_wh)
    
    count_and, count_or, count_vis, count_cod = COUNT(img_visium_mask_wh, img_phenocycler_mask_resized_zoomed_shifted_wh)
    cr = count_and / count_or
    xor = (count_or - count_and)

    print("<h3>", file=fp_html)
    print("Translation", file=fp_html)
    print("</h3>", file=fp_html)

    print("(dx, dy) = ({:d}, {:d})<br>".format(dx, dy), file=fp_html)
    print("IoU = {0:.4f} = {1:d} / {2:d} = {1:d} / ({3:d} + {4:d} - {1:d}) = |V &cap; P| / (|V| + |P| - |V &cap; P|) <br>"\
          .format(cr, count_and, count_or, count_vis, count_cod), file=fp_html)
    print("XOR = {0:d} = {1:d} - {2:d} = |V &cup; P| - |V &cap; P| <br>".format(xor, count_or, count_and), file=fp_html)    

    print("<img src=" + "\"" + filename_shift_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    print("<img src=" + "\"" + filename_shift_wh_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    fp_html.flush()


    ## ROTATE

    rr = 1.0
    rn = int(180 / rr)

    max_cr = -999999.0
    max_count_and = 0
    max_count_or = 0        
    max_count_vis = 0
    max_count_cod = 0
    opt_dg = 0
    min_xor = 999999

    max_img = img_phenocycler_mask_resized_zoomed_shifted
    max_img_wh = img_phenocycler_mask_resized_zoomed_shifted_wh

    for kr in range(-rn, rn+1):

        dg = kr * rr

        img_phenocycler_opt_rr, img_phenocycler_opt_rr_wh = \
            ROTATE(img_phenocycler_mask_resized_zoomed_shifted, img_phenocycler_mask_resized_zoomed_shifted_wh, xc_vis, yc_vis, h_max, w_max, dg)

        count_and, count_or, count_vis, count_cod = COUNT(img_visium_mask_resized_wh, img_phenocycler_opt_rr_wh)        

        cr = count_and / count_or
        xor = (count_or - count_and)

        if cr > max_cr:
            min_xor = xor
            max_cr = cr
            max_count_and = count_and
            max_count_or = count_or
            max_count_vis = count_vis
            max_count_cod = count_cod
            opt_dg = dg

            max_img = img_phenocycler_opt_rr
            max_img_wh = img_phenocycler_opt_rr_wh         

            print("ROTATE", "{:.4f}".format(cr), "{:.3f}".format(dg))            

    degree = opt_dg
    
    img_phenocycler_mask_resized_zoomed_shifted_rotated = max_img
    img_phenocycler_mask_resized_zoomed_shifted_rotated_wh = max_img_wh

    img_and = cv2.bitwise_and(img_phenocycler_mask_resized_zoomed_shifted_rotated, img_visium_mask_resized)

    img_rotate = np.ones((h_max, w_max)) * 255
    img_rotate = cv2.cvtColor(img_rotate.astype(np.uint8), cv2.COLOR_GRAY2BGR)
    img_rotate[np.where(img_visium_mask_resized == 255)] = orangered    
    img_rotate[np.where(img_phenocycler_mask_resized_zoomed_shifted_rotated == 255)] = royalblue
    img_rotate[np.where(img_and == 255)] = color_and
    img_rotate = bg_transp(img_rotate)
    
    filename_rotate = ALIGNDIR_DAT + "/" + sid + "_rotate" + ".png"
    filename_rotate_rel = "dat" + "/" + sid + "_rotate" + ".png"
    cv2.imwrite(filename_rotate, img_rotate)

    # with hole

    img_and_wh = cv2.bitwise_and(img_phenocycler_mask_resized_zoomed_shifted_rotated_wh, img_visium_mask_resized_wh)    

    img_rotate_wh = np.ones((h_max, w_max)) * 255
    img_rotate_wh = cv2.cvtColor(img_rotate_wh.astype(np.uint8), cv2.COLOR_GRAY2BGR)
    img_rotate_wh[np.where(img_visium_mask_resized_wh == 255)] = orangered    
    img_rotate_wh[np.where(img_phenocycler_mask_resized_zoomed_shifted_rotated_wh == 255)] = royalblue
    img_rotate_wh[np.where(img_and_wh == 255)] = color_and 
    img_rotate_wh = bg_transp(img_rotate_wh)
    
    filename_rotate_wh = ALIGNDIR_DAT + "/" + sid + "_rotate_wh" + ".png"
    filename_rotate_wh_rel = "dat" + "/" + sid + "_rotate_wh" + ".png"
    cv2.imwrite(filename_rotate_wh, img_rotate_wh)

    print("<h3>", file=fp_html)
    print("Rotation", file=fp_html)
    print("</h3>", file=fp_html)

    print("Degree(anticlockwise) = {:.1f}@({:d}, {:d})<br>".format(opt_dg, xc_vis, yc_vis), file=fp_html)

    print("IoU = {0:.4f} = {1:d} / {2:d} = {1:d} / ({3:d} + {4:d} - {1:d}) = |V &cap; P| / (|V| + |P| - |V &cap; P|) <br>"\
          .format(max_cr, max_count_and, max_count_or, max_count_vis, max_count_cod), file=fp_html)
    print("XOR = {0:d} = {1:d} - {2:d} = |V &cup; P| - |V &cap; P| <br>".format(min_xor, max_count_or, max_count_and), file=fp_html)

    print("<img src=" + "\"" + filename_rotate_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    print("<img src=" + "\"" + filename_rotate_wh_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    fp_html.flush()
   
    img_phenocycler_opt = img_phenocycler_mask_resized_zoomed_shifted_rotated    
    img_phenocycler_opt_wh = img_phenocycler_mask_resized_zoomed_shifted_rotated_wh
    
    count_and, count_or, count_vis, count_cod = COUNT(img_visium_mask_wh, img_phenocycler_opt_wh)

    
    ######################
    ######################

    # OPTIMIZE (exhaustive grid search w/o pruning)

    print("===== OPTIMIZATION =====")    

    max_count_and = 0
    max_count_or = 0
    max_count_vis = 0
    max_count_cod = 0

    max_pr = 0
    max_px = 0
    max_py = 0
    max_pd = 0
    max_dxc = 0
    max_dyc = 0
    max_img = img_phenocycler_opt_wh
    
    if sid == "FFPE_LUAD_2_C":
        zz = 0.010
        zw = 0.080
        zn = int(zw / zz) 
        sn = 8
        ss = 2
        rn = 5
        rr = 0.5
        
    elif sid == "FFPE_LUAD_4_C":
        zz = 0.010
        zw = 0.050
        zn = int(zw / zz) 
        sn = 5
        ss = 10
        rn = 5
        rr = 1

    else: # "FFPE_LUAD_3_B"
        zz = 0.001
        zw = 0.012
        zn = int(zw / zz)
        sn = 5
        ss = 1
        rn = 10
        rr = 0.05

    zz_list = [0]
    for z in range(1, zn+1):
        zz_list.append(z)
        zz_list.append(-1*z)

    xx_list = [0]        
    for x in range(1, sn+1):
        xx_list.append(x)
        xx_list.append(-1*x)

    rr_list = [0]
    for r in range(1, rn+1):
        rr_list.append(r)
        rr_list.append(-1*r)

        
    DO_OPTIMIZE = True
    if DO_OPTIMIZE:

        ### ZOOM
        for kz in zz_list:

            pr = 1.0 + kz * zz

            print("ZOOM", "{:.3f}".format(pr))

            img_phenocycler_mask_resized_tmp, img_phenocycler_opt_wh_tmp = \
                ZOOM(img_phenocycler_opt, img_phenocycler_opt_wh, h_max, w_max, pr)

            img_phenocycler_opt_zz, img_phenocycler_opt_zz_wh, dxc, dyc = \
                CENTERING(img_phenocycler_mask_resized_tmp, img_phenocycler_opt_wh_tmp, xc_vis, yc_vis)

            ### SHIFT            
            for kx in xx_list:
                for ky in xx_list:

                    px = kx * ss
                    py = ky * ss

                    img_phenocycler_opt_ss, img_phenocycler_opt_ss_wh = \
                        SHIFT(img_phenocycler_opt_zz, img_phenocycler_opt_zz_wh, px, py)

                    ### ROTATE
                    for kr in rr_list:

                        dg = kr * rr

                        img_phenocycler_opt_rr, img_phenocycler_opt_rr_wh = \
                            ROTATE(img_phenocycler_opt_ss, img_phenocycler_opt_ss_wh, xc_vis+px, yc_vis+py, h_max, w_max, dg)

                        count_and, count_or, count_vis, count_cod = COUNT(img_visium_mask_wh, img_phenocycler_opt_rr_wh)
                        cr = count_and / count_or
                        xor = (count_or - count_and)
                        
                        if cr > max_cr:

                            (max_cr, min_xor, max_pr, max_px, max_py, max_pd, max_dxc, max_dyc) = (cr, xor, pr, px, py, dg, dxc, dyc)
                            (max_count_and, max_count_or, max_count_vis, max_count_cod) = (count_and, count_or, count_vis, count_cod)

                            print("  IMPROVED", "{:.4f}".format(cr), "{:.3f}".format(pr), px, py, "{:.3f}".format(dg))
                            
                            max_img = img_phenocycler_opt_rr                        
                            max_img_wh = img_phenocycler_opt_rr_wh
                            

        #####

        img_phenocycler_mask_resized_shifted_rotated_opt = max_img
        img_phenocycler_mask_resized_shifted_rotated_opt_wh = max_img_wh

        img_and = cv2.bitwise_and(img_phenocycler_mask_resized_shifted_rotated_opt, img_visium_mask_resized)

        img_opt = np.ones((h_max, w_max)) * 255
        img_opt = cv2.cvtColor(img_opt.astype(np.uint8), cv2.COLOR_GRAY2BGR)
        img_opt[np.where(img_visium_mask_resized == 255)] = orangered    
        img_opt[np.where(img_phenocycler_mask_resized_shifted_rotated_opt == 255)] = royalblue
        img_opt[np.where(img_and == 255)] = color_and 
        img_opt = bg_transp(img_opt)

        filename_opt = ALIGNDIR_DAT + "/" + sid + "_opt" + ".png"
        filename_opt_rel = "dat" + "/" + sid + "_opt" + ".png"
        cv2.imwrite(filename_opt, img_opt)

        # with hole

        img_and_wh = cv2.bitwise_and(img_phenocycler_mask_resized_shifted_rotated_opt_wh, img_visium_mask_resized_wh)    

        img_opt_wh = np.ones((h_max, w_max)) * 255
        img_opt_wh = cv2.cvtColor(img_opt_wh.astype(np.uint8), cv2.COLOR_GRAY2BGR)
        img_opt_wh[np.where(img_visium_mask_wh == 255)] = orangered    
        img_opt_wh[np.where(img_phenocycler_mask_resized_shifted_rotated_opt_wh == 255)] = royalblue
        img_opt_wh[np.where(img_and_wh == 255)] = color_and #darkmagenta #darkgreen
        img_opt_wh = bg_transp(img_opt_wh)

        filename_opt_wh = ALIGNDIR_DAT + "/" + sid + "_opt_wh" + ".png"
        filename_opt_wh_rel = "dat" + "/" + sid + "_opt_wh" + ".png"
        cv2.imwrite(filename_opt_wh, img_opt_wh)
        
        # cache to skip optimization

        filename_opt_params_txt = ALIGNDIR_DAT + "/" + sid + "_opt_params" + ".txt"
        fp_opt = open(filename_opt_params_txt, "wt")

        print("{}\t{}".format("max_cr", max_cr), file=fp_opt)
        print("{}\t{}".format("max_pr", max_pr), file=fp_opt)
        print("{}\t{}".format("max_px", max_px), file=fp_opt)
        print("{}\t{}".format("max_py", max_py), file=fp_opt)
        print("{}\t{}".format("max_pd", max_pd), file=fp_opt)
        print("{}\t{}".format("max_dxc", max_dxc), file=fp_opt)
        print("{}\t{}".format("max_dyc", max_dyc), file=fp_opt)        
        
        fp_opt.close()
                        
    else: # LOAD

        filename_opt_params_txt = ALIGNDIR_DAT + "/" + sid + "_opt_params" + ".txt"
        fp_opt = open(filename_opt_params_txt, "rt")

        params = []
        count = 0
        for line in fp_opt:
            line = line.rstrip()
            vals = line.split("\t")

            print(line)

            kk = vals[0]
            vv = vals[1]
            
            params.append(vv)

        max_cr = float(params[0])
        max_pr = float(params[1])
        max_px = int(params[2])
        max_py = int(params[3])
        max_pd = float(params[4])
        max_dxc = int(params[5])
        max_dyc = int(params[6])

        print(max_cr, max_pr, max_px, max_py, max_pd, max_dxc, max_dyc)

        filename_opt = ALIGNDIR_DAT + "/" + sid + "_opt" + ".png"
        filename_opt_rel = "dat" + "/" + sid + "_opt" + ".png"
        
        filename_opt_wh = ALIGNDIR_DAT + "/" + sid + "_opt_wh" + ".png"
        filename_opt_wh_rel = "dat" + "/" + sid + "_opt_wh" + ".png"

        img_opt = cv2.imread(filename_opt, cv2.IMREAD_UNCHANGED)
        img_opt_wh = cv2.imread(filename_opt_wh, cv2.IMREAD_UNCHANGED)

    #####
    
    print("<h3>", file=fp_html)
    print("Optimization", file=fp_html)
    print("</h3>", file=fp_html)

    print("PhenoCycler &times; {:.4f}<br>".format(max_pr), file=fp_html)
    print("(dx, dy) = ({:d}, {:d})<br>".format(max_dxc, max_dyc), file=fp_html)    
    print("(dx, dy) = ({:d}, {:d})<br>".format(max_px, max_py), file=fp_html)
    print("Degree(anticlockwise) = {:.1f}@({:d}, {:d})<br>".format(max_pd, xc_cod+max_px, yc_cod+max_py), file=fp_html)

    print("IoU = {0:.4f} = {1:d} / {2:d} = {1:d} / ({3:d} + {4:d} - {1:d}) = |V &cap; P| / (|V| + |P| - |V &cap; P|) <br>".\
          format(max_cr, max_count_and, max_count_or, max_count_vis, max_count_cod), file=fp_html)
    print("XOR = {0:d} = {1:d} - {2:d} = |V &cup; P| - |V &cap; P| <br>".format(min_xor, max_count_or, max_count_and), file=fp_html)    
    
    print("<img src=" + "\"" + filename_opt_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    print("<img src=" + "\"" + filename_opt_wh_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    fp_html.flush()

    #####
    
    print("<h3>", file=fp_html)
    print("Frames for animation", file=fp_html)
    print("</h3>", file=fp_html)

    html_ss = ""
    html_ss_wh = ""

    FILES = []
    FILES_wh = []
    
    ####
    # Transition
    # from the initial to the optimum

    uh = 2.0 # height of the image : h_max * uh
    uw = 1.1
    h_max_tr = int(h_max * uh)
    w_max_tr = int(w_max * uw)    

    ww_u = "20px"
    
    img_v = img_visium_mask_resized
    img_p = img_phenocycler_mask_resized

    img_v_wh = img_visium_mask_resized_wh
    img_p_wh = img_phenocycler_mask_resized_wh

    def plot(img_transition, img_xv, img_xp, img_and):

        img_transition[np.where(img_xv == 255)] = orangered            
        img_transition[np.where(img_xp == 255)] = royalblue
        img_transition[np.where(img_and == 255)] = color_and 
        img_transition = bg_transp(img_transition)

    # and
    img_pp = np.zeros((h_max_tr, w_max_tr), dtype=np.uint8) 
    img_pp[np.where(img_p == 255)] = 255

    img_pp_wh = np.zeros((h_max_tr, w_max_tr), dtype=np.uint8) 
    img_pp_wh[np.where(img_p_wh == 255)] = 255
    
    M = np.float32([[1,0,0],[0,1,int(0.0*h_max)]])
    img_pp_x = cv2.warpAffine(img_pp, M, (w_max_tr, h_max_tr))

    # 
    M = np.float32([[1,0,0],[0,1,int(1.0*h_max)]])
    img_vv = cv2.warpAffine(img_v, M, (w_max_tr, h_max_tr))
    img_vv_wh = cv2.warpAffine(img_v_wh, M, (w_max_tr, h_max_tr))
    
    img_and = cv2.bitwise_and(img_pp, img_vv)
    img_and_wh = cv2.bitwise_and(img_pp_wh, img_vv_wh)
    
    # 
    M = np.float32([[1,0,0],[0,1,int(0.0*h_max)]])
    img_vv_0 = cv2.warpAffine(img_v, M, (w_max_tr, h_max_tr))
    img_vv_0_wh = cv2.warpAffine(img_v_wh, M, (w_max_tr, h_max_tr))
    

    #####
    # INIT
    count = 0
    # transition : background/base image
    img_transition = np.ones((h_max_tr, w_max_tr), dtype=np.uint8) * 255 
    img_transition = cv2.cvtColor(img_transition, cv2.COLOR_GRAY2BGR)
    plot(img_transition, img_vv, img_pp, img_and)

    # wh
    img_transition_wh = np.ones((h_max_tr, w_max_tr), dtype=np.uint8) * 255 
    img_transition_wh = cv2.cvtColor(img_transition_wh, cv2.COLOR_GRAY2BGR)
    plot(img_transition_wh, img_vv_wh, img_pp_wh, img_and_wh)

    filename_transition = ALIGNDIR_DAT_PPMNG + "/" + sid + "_transition_{:02d}".format(count) + ".png"
    filename_transition_rel = "dat/ppmng" + "/" + sid + "_transition_{:02d}".format(count) + ".png"

    cv2.imwrite(filename_transition, img_transition)

    html_ss += "<img src=" + "\"" + filename_transition_rel + "\"" + " width=\"" + ww_u + "\"" + " border=\"1\"" + ">" + "\n"

    # wh
    filename_transition_wh = ALIGNDIR_DAT_PPMNG + "/" + sid + "_transition_wh_{:02d}".format(count) + ".png"
    filename_transition_wh_rel = "dat/ppmng" + "/" + sid + "_transition_wh_{:02d}".format(count) + ".png"
    
    cv2.imwrite(filename_transition_wh, img_transition_wh)
    html_ss_wh += "<img src=" + "\"" + filename_transition_wh_rel + "\"" + " width=\"" + ww_u + "\"" + " border=\"1\"" + ">" + "\n"

    FILES.append(filename_transition)
    FILES_wh.append(filename_transition_wh)    

    # ZOOM
    num_zoom = 10
    for k in range(1, num_zoom):
        ratio_dd = 1.0 - k * ((1.0 - ratio) / num_zoom)

        img_pp_x, img_pp_x_wh = ZOOM(img_pp, img_pp_wh, h_max_tr, w_max_tr, ratio_dd)
        
        img_and = cv2.bitwise_and(img_pp_x, img_vv)
        img_and = cv2.bitwise_and(img_pp_x_wh, img_vv_wh)        
        
        count += 1
        # transition : background/base image    
        img_transition = np.ones((h_max_tr, w_max_tr), dtype=np.uint8) * 255 # base
        img_transition = cv2.cvtColor(img_transition, cv2.COLOR_GRAY2BGR)
        plot(img_transition, img_vv, img_pp_x, img_and)

        # wh
        img_transition_wh = np.ones((h_max_tr, w_max_tr), dtype=np.uint8) * 255 # base
        img_transition_wh = cv2.cvtColor(img_transition_wh, cv2.COLOR_GRAY2BGR)
        plot(img_transition_wh, img_vv_wh, img_pp_x_wh, img_and_wh)
        

        filename_transition = ALIGNDIR_DAT_PPMNG + "/" + sid + "_transition_{:02d}".format(count) + ".png"
        filename_transition_rel = "dat/ppmng" + "/" + sid + "_transition_{:02d}".format(count) + ".png"
        
        cv2.imwrite(filename_transition, img_transition)
        html_ss += "<img src=" + "\"" + filename_transition_rel + "\"" + " width=\"" + ww_u + "\"" + " border=\"1\"" + ">" + "\n"

        # wh
        filename_transition_wh = ALIGNDIR_DAT_PPMNG + "/" + sid + "_transition_wh_{:02d}".format(count) + ".png"
        filename_transition_wh_rel = "dat/ppmng" + "/" + sid + "_transition_wh_{:02d}".format(count) + ".png"
        
        cv2.imwrite(filename_transition_wh, img_transition_wh)
        html_ss_wh += "<img src=" + "\"" + filename_transition_wh_rel + "\"" + " width=\"" + ww_u + "\"" + " border=\"1\"" + ">" + "\n"

        FILES.append(filename_transition)
        FILES_wh.append(filename_transition_wh)    

        
    img_pp, img_pp_wh = ZOOM(img_pp, img_pp_wh, h_max_tr, w_max_tr, ratio)
    img_and = cv2.bitwise_and(img_pp, img_vv)

    img_and_wh = cv2.bitwise_and(img_pp_wh, img_vv_wh)    

    count += 1
    # transition : background/base image    
    img_transition = np.ones((h_max_tr, w_max_tr), dtype=np.uint8) * 255 # base
    img_transition = cv2.cvtColor(img_transition, cv2.COLOR_GRAY2BGR)
    plot(img_transition, img_vv, img_pp, img_and)

    # wh
    img_transition_wh = np.ones((h_max_tr, w_max_tr), dtype=np.uint8) * 255 # base
    img_transition_wh = cv2.cvtColor(img_transition_wh, cv2.COLOR_GRAY2BGR)
    plot(img_transition_wh, img_vv_wh, img_pp_wh, img_and_wh)

    # solid
    filename_transition = ALIGNDIR_DAT_PPMNG + "/" + sid + "_transition_{:02d}".format(count) + ".png"
    filename_transition_rel = "dat/ppmng" + "/" + sid + "_transition_{:02d}".format(count) + ".png"
    
    cv2.imwrite(filename_transition, img_transition)
    html_ss += "<img src=" + "\"" + filename_transition_rel + "\"" + " width=\"" + ww_u + "\"" + " border=\"1\"" + ">" + "\n"

    # wh
    filename_transition_wh = ALIGNDIR_DAT_PPMNG + "/" + sid + "_transition_wh_{:02d}".format(count) + ".png"
    filename_transition_wh_rel = "dat/ppmng" + "/" + sid + "_transition_wh_{:02d}".format(count) + ".png"
    
    cv2.imwrite(filename_transition_wh, img_transition_wh)
    html_ss_wh += "<img src=" + "\"" + filename_transition_wh_rel + "\"" + " width=\"" + ww_u + "\"" + " border=\"1\"" + ">" + "\n"

    FILES.append(filename_transition)
    FILES_wh.append(filename_transition_wh)    
    
    img_pp, img_pp_wh = SHIFT(img_pp, img_pp_wh, dx, dy)

    # ROTATE
    num_rotate = 20
    degree_dd = degree / num_rotate
    center_dd = 0.5 * h_max / num_rotate    
    for k in range(1, num_rotate+1):
        degree_k = k * degree_dd
        img_pp_x, img_pp_x_wh = ROTATE(img_pp, img_pp_wh, xc_vis, yc_vis, h_max_tr, w_max_tr, degree_k)

        center_k = k * center_dd

        img_pp_x, img_pp_x_wh = SHIFT(img_pp_x, img_pp_x_wh, 0, center_k)            
        img_vv_x, img_vv_x_wh = SHIFT(img_vv, img_vv_wh, 0,  -center_k)
        
        img_and = cv2.bitwise_and(img_pp_x, img_vv_x)
        img_and_wh = cv2.bitwise_and(img_pp_x_wh, img_vv_x_wh)        

        count += 1
        # transition : background/base image    
        img_transition = np.ones((h_max_tr, w_max_tr), dtype=np.uint8) * 255 # base
        img_transition = cv2.cvtColor(img_transition, cv2.COLOR_GRAY2BGR)
        plot(img_transition, img_vv_x, img_pp_x, img_and)

        # wh
        img_transition_wh = np.ones((h_max_tr, w_max_tr), dtype=np.uint8) * 255 # base
        img_transition_wh = cv2.cvtColor(img_transition_wh, cv2.COLOR_GRAY2BGR)
        plot(img_transition_wh, img_vv_x_wh, img_pp_x_wh, img_and_wh)

        # solid
        filename_transition = ALIGNDIR_DAT_PPMNG + "/" + sid + "_transition_{:02d}".format(count) + ".png"
        filename_transition_rel = "dat/ppmng" + "/" + sid + "_transition_{:02d}".format(count) + ".png"
        
        cv2.imwrite(filename_transition, img_transition)
        html_ss += "<img src=" + "\"" + filename_transition_rel + "\"" + " width=\"" + ww_u + "\"" + " border=\"1\"" + ">" + "\n"

        # wh
        filename_transition_wh = ALIGNDIR_DAT_PPMNG + "/" + sid + "_transition_wh_{:02d}".format(count) + ".png"
        filename_transition_wh_rel = "dat/ppmng" + "/" + sid + "_transition_wh_{:02d}".format(count) + ".png"
        
        cv2.imwrite(filename_transition_wh, img_transition_wh)
        html_ss_wh += "<img src=" + "\"" + filename_transition_wh_rel + "\"" + " width=\"" + ww_u + "\"" + " border=\"1\"" + ">" + "\n"

        FILES.append(filename_transition)
        FILES_wh.append(filename_transition_wh)    
        
    img_pp, img_pp_wh = ROTATE(img_pp, img_pp_wh, xc_vis, yc_vis, h_max_tr, w_max_tr, degree)    

    # Calibration

    img_pp, img_pp_wh = ZOOM(img_pp, img_pp_wh, h_max_tr, w_max_tr, max_pr)
    img_pp, img_pp_wh = SHIFT(img_pp, img_pp_wh, max_dxc, max_dyc)    
    img_pp, img_pp_wh = SHIFT(img_pp, img_pp_wh, max_px, max_py)
    img_pp, img_pp_wh = ROTATE(img_pp, img_pp_wh, xc_vis+max_px, yc_vis+max_py, h_max_tr, w_max_tr, max_pd)


    #####

    img_and = cv2.bitwise_and(img_pp, img_vv_0)
    img_and_wh = cv2.bitwise_and(img_pp_wh, img_vv_0_wh)    
    
    img_pp_x, img_pp_x_wh = SHIFT(img_pp, img_pp_wh, 0, 0.5*h_max)    
    img_vv_x, img_vv_x_wh = SHIFT(img_vv_0, img_vv_0_wh, 0, 0.5*h_max) 
    img_and_x, _ = SHIFT(img_and, img_and, 0, 0.5*h_max)

    img_and_x_wh, _ = SHIFT(img_and_wh, img_and_wh, 0, 0.5*h_max)    

    count += 1
    # transition : background/base image    
    img_transition = np.ones((h_max_tr, w_max_tr), dtype=np.uint8) * 255 
    img_transition = cv2.cvtColor(img_transition, cv2.COLOR_GRAY2BGR)
    plot(img_transition, img_vv_x, img_pp_x, img_and_x)

    # wh
    img_transition_wh = np.ones((h_max_tr, w_max_tr), dtype=np.uint8) * 255 
    img_transition_wh = cv2.cvtColor(img_transition_wh, cv2.COLOR_GRAY2BGR)
    plot(img_transition_wh, img_vv_x_wh, img_pp_x_wh, img_and_x_wh)
    
    # solid
    filename_transition = ALIGNDIR_DAT_PPMNG + "/" + sid + "_transition_{:02d}".format(count) + ".png"
    filename_transition_rel = "dat/ppmng" + "/" + sid + "_transition_{:02d}".format(count) + ".png"
    
    cv2.imwrite(filename_transition, img_transition)
    html_ss += "<img src=" + "\"" + filename_transition_rel + "\"" + " width=\"" + ww_u + "\"" + " border=\"1\"" + ">" + "\n"

    # wh
    filename_transition_wh = ALIGNDIR_DAT_PPMNG + "/" + sid + "_transition_wh_{:02d}".format(count) + ".png"
    filename_transition_wh_rel = "dat/ppmng" + "/" + sid + "_transition_wh_{:02d}".format(count) + ".png"
    
    cv2.imwrite(filename_transition_wh, img_transition_wh)
    html_ss_wh += "<img src=" + "\"" + filename_transition_wh_rel + "\"" + " width=\"" + ww_u + "\"" + " border=\"1\"" + ">" + "\n"

    FILES.append(filename_transition)
    FILES_wh.append(filename_transition_wh)    


    #####

    print(html_ss, file=fp_html)
    print("<br>", file=fp_html)    
    print(html_ss_wh, file=fp_html)

    # solid

    IMGS = []
    for filename in FILES:
        #print(filename)

        img = cv2.imread(filename)
        IMGS.append(img)
        
    h, w = IMGS[0].shape[:2]

    h_s = int(h * (w_s / w))

    filename_mp4 = ALIGNDIR_DAT_PPMNG + "/" + sid + "_transition" + ".mp4"    
    filename_mp4_s = ALIGNDIR_DAT_PPMNG + "/" + sid + "_transition_{:d}".format(w_s) + ".mp4"    

    video = cv2.VideoWriter(filename_mp4, cv2.VideoWriter_fourcc('m', 'p', '4', 'v'), 10.0, (w, h))
    video_s = cv2.VideoWriter(filename_mp4_s, cv2.VideoWriter_fourcc('m', 'p', '4', 'v'), 10.0, (w_s, h_s))
    
    for img in IMGS:

        img_s = cv2.resize(img, dsize=(w_s, h_s), interpolation=cv2.INTER_NEAREST)

        video.write(img)
        video_s.write(img_s)        

    video.release()
    video_s.release()

    # wh 
    
    IMGS_wh = []
    for filename in FILES_wh:
        #print(filename)

        img = cv2.imread(filename)
        IMGS_wh.append(img)
        
    h, w = IMGS_wh[0].shape[:2]

    h_s = int(h * (w_s / w))

    filename_mp4_wh = ALIGNDIR_DAT_PPMNG + "/" + sid + "_transition_wh" + ".mp4"        

    filename_mp4_wh_s = ALIGNDIR_DAT_PPMNG + "/" + sid + "_transition_wh_{:d}".format(w_s) + ".mp4"        
    
    video = cv2.VideoWriter(filename_mp4_wh, cv2.VideoWriter_fourcc('m', 'p', '4', 'v'), 10.0, (w, h))
    video_s = cv2.VideoWriter(filename_mp4_wh_s, cv2.VideoWriter_fourcc('m', 'p', '4', 'v'), 10.0, (w_s, h_s))
    
    for img in IMGS_wh:

        img_s = cv2.resize(img, dsize=(w_s, h_s), interpolation=cv2.INTER_NEAREST)
            
        video.write(img)
        video_s.write(img_s)
        
    video.release()
    video_s.release()
    
    print_html_footer(fp_html)    
        
    fp_html.close()



    #####
    
    fp = open(filename_txt, "wt")

    # PhenoCycler, FULL->HIGH
    print("ZOOM\t{:.4f}".format(0.1), file=fp)
    # INITIALIZATION
    print("SHIFT\t{:d}\t{:d}".format(dm0, dm0), file=fp)    
    print("ZOOM\t{:.4f}".format(ratio), file=fp)
    print("SHIFT\t{:d}\t{:d}".format(dx, dy), file=fp)
    print("ROTATE\t{:.2f}\t{:d}\t{:d}".format(degree, xc_vis, yc_vis), file=fp)
    # OPTIMIZATION
    print("ZOOM\t{:.4f}".format(max_pr), file=fp)
    print("SHIFT\t{:d}\t{:d}".format(max_dxc, max_dyc), file=fp) # CENTERING
    print("SHIFT\t{:d}\t{:d}".format(max_px, max_py), file=fp)
    print("ROTATE\t{:.2f}\t{:d}\t{:d}".format(max_pd, xc_vis+max_px, yc_vis+max_py), file=fp)
    print("SHIFT\t{:d}\t{:d}".format(-dm0, -dm0), file=fp)        

    fp.close()

    
#####

if __name__ == "__main__":

    main()

    
