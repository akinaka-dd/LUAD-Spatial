#!/usr/bin/env python3

import sys
import re
import pathlib
import shutil
import numpy as np
import cv2

import matplotlib.pyplot as plt

from util.image import *
from util.html import *

def main():

    # commad line arguments
    if len(sys.argv) < 4:
        print("ERROR: Command line arguments")
        sys.exit()

    # directory of SpaceRanger output
    dirname_spaceranger = sys.argv[1]

    # path to the full-resolution image
    filename_tif_ful_00 = sys.argv[2]

    # sample name
    sid = sys.argv[3]

    # path to the HE-stained image ("hires")
    filename_png_00 = dirname_spaceranger + "/" + "spatial/tissue_hires_image.png"

    # rotate and flip
    degree = ""
    if len(sys.argv) > 4:
        degree = sys.argv[3]

    flip = ""
    if len(sys.argv) > 4:
        flip = sys.argv[4]

    DEGREE = ""
    if degree == "ROTATE_90_CLOCKWISE":
        DEGREE = cv2.ROTATE_90_CLOCKWISE
    elif degree == "ROTATE_90_COUNTERCLOCKWISE":
        DEGREE = cv2.ROTATE_90_COUNTERCLOCKWISE

    FLIP = ""
    if flip == "FLIP_LR": # horizontal 
        FLIP = +1
    elif flip == "FLIP_UD": # vertical
        FLIP = 0
    elif flip == "FLIP_UD_LR":# both
        FLIP = -1

    # directories
    SPACEDIR = dirname_spaceranger
    WORKDIR = "."
    RESULTDIR = WORKDIR + "/" + sid
    pp = pathlib.Path(RESULTDIR)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)
    PREFIX = sid + "-" + "Visium"
    MASKDIR = RESULTDIR + "/" + PREFIX + "-MASK"
    MASKDIR_DAT = MASKDIR + "/" + "dat"    
    pp = pathlib.Path(MASKDIR_DAT)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)

    # load image
    img_00 = cv2.imread(filename_png_00)
    pp = pathlib.Path(filename_png_00)
    pp_name = pp.name
    pp_stem = pp.stem

    img_ful_00 = cv2.imread(filename_tif_ful_00)

    # html file
    filename_html = MASKDIR + "/" + PREFIX + ".html"        
    fp_html = open(filename_html, "wt")

    ww = "200px"
    print_html_header(fp_html)

    print("<h1>", file=fp_html)
    print("Visium SpaceRanger spatial" + " &mdash; " + sid, file=fp_html)
    print("</h1>", file=fp_html)

    #####

    # Copy the scalefactors(JSON) and position_list(CSV) files

    filename_sf = SPACEDIR + "/" + "spatial/scalefactors_json.json"
    filename_sf_prep = MASKDIR_DAT + "/" + "scalefactors_json.json"
    shutil.copy2(filename_sf, filename_sf_prep)
    print("COPY", filename_sf_prep)

    # # Space Ranger 1.3.0
    filename_pl = SPACEDIR + "/" + "spatial/tissue_positions_list.csv"
    filename_pl_prep = MASKDIR_DAT + "/" + "tissue_positions_list.csv"
    # Space Ranger 2.0.0   
    #filename_pl = SPACEDIR + "/" + "spatial/tissue_positions.csv"
    #filename_pl_prep = MASKDIR_DAT + "/" + "tissue_positions.csv"
    shutil.copy2(filename_pl, filename_pl_prep)
    print("COPY", filename_pl_prep)

    filename_ft = SPACEDIR + "/" + "filtered_feature_bc_matrix/features.tsv"
    filename_ft_prep = MASKDIR_DAT + "/" + "features.tsv"
    shutil.copy2(filename_ft, filename_ft_prep)
    print("COPY", filename_ft_prep)

    # Copy the original
    filename_png_dup = MASKDIR_DAT + "/" + PREFIX + "_00" + ".png"            
    shutil.copy2(filename_png_00, filename_png_dup)
    print("COPY", filename_png_dup)

    filename_tif_ful_dup = MASKDIR_DAT + "/" + PREFIX + "_ful_00" + ".tif"        
    shutil.copy2(filename_tif_ful_00, filename_tif_ful_dup)
    print("COPY", filename_tif_ful_dup)    

    print("<h2>", file=fp_html)
    print("Hires original", file=fp_html)
    print("</h2>", file=fp_html)

    filename_png_dup = MASKDIR_DAT + "/" + PREFIX + "_00" + ".png"
    filename_png_dup_rel = "dat" + "/" + PREFIX + "_00" + ".png"        
    print("<img src=" + "\"" + filename_png_dup_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    print("<h2>", file=fp_html)
    print("Visium masks", file=fp_html)
    print("</h2>", file=fp_html)

    img = img_00
    img_ful = img_ful_00
    if not DEGREE == "":
        img = cv2.rotate(img_00, DEGREE)
        img_ful = cv2.rotate(img_ful_00, DEGREE)
        
    if not FLIP == "":
        img = cv2.flip(img, FLIP)
        img_ful = cv2.flip(img_ful, FLIP)
        
    filename_png = MASKDIR_DAT + "/" + PREFIX + ".png"        
    cv2.imwrite(filename_png, img)

    filename_png_ful = MASKDIR_DAT + "/" + PREFIX + "_ful.png"    
    cv2.imwrite(filename_png_ful, img_ful)
    

    print("<h3>", file=fp_html)
    print("Rotating/Flipping [OPTIONAL]", file=fp_html)        
    print("</h3>", file=fp_html)

    filename_png = MASKDIR_DAT + "/" + PREFIX + ".png"
    filename_png_rel = "dat" + "/" + PREFIX + ".png"
    print("<img src=" + "\"" + filename_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    # Trimming of the edges
    img_tr = img.copy()    
    trim_dd = 20
    if sid == "FFPE_LUAD_2_C":
        img_tr[:, -100:] = (220, 220, 220)
    img_tr[:trim_dd, :] = (255, 255, 255)
    img_tr[-trim_dd:, :] = (255, 255, 255)
    img_tr[:, :trim_dd] = (255, 255, 255)
    img_tr[:, -trim_dd:] = (255, 255, 255)
    filename_tr = MASKDIR_DAT + "/" + PREFIX + "_tr.png"        
    cv2.imwrite(filename_tr, img_tr)
    print("TRIMMING", filename_tr)        

    # grayscale
    img_gray = cv2.cvtColor(img_tr, cv2.COLOR_BGR2GRAY) # BGR -> GRAYSCALE
    filename_gray = MASKDIR_DAT + "/" + PREFIX + "_tr_gray.png"        
    cv2.imwrite(filename_gray, img_gray)
    print("GRAYSCALE", filename_gray)        

    print("<h3>", file=fp_html)
    print("Trimming", file=fp_html)
    print("</h3>", file=fp_html)

    filename_tr = MASKDIR_DAT + "/" + PREFIX + "_tr.png"
    filename_tr_rel = "dat" + "/" + PREFIX + "_tr.png"
    print("<img src=" + "\"" + filename_tr_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    filename_gray = MASKDIR_DAT + "/" + PREFIX + "_tr_gray.png"
    filename_gray_rel = "dat" + "/" + PREFIX + "_tr_gray.png"
    print("<img src=" + "\"" + filename_gray_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    

    # Circle detection

    if sid == "FFPE_LUAD_2_C":
        circles = cv2.HoughCircles(img_gray, cv2.HOUGH_GRADIENT, dp=1, minDist=10, param1=25, param2=15, minRadius=6, maxRadius=10)
    elif sid == "FFPE_LUAD_4_C":
        circles = cv2.HoughCircles(img_gray, cv2.HOUGH_GRADIENT, dp=1, minDist=10, param1=25, param2=15, minRadius=8, maxRadius=10)
    else: # "FFPE_LUAD_3_B"
        circles = cv2.HoughCircles(img_gray, cv2.HOUGH_GRADIENT, dp=1, minDist=10, param1=25, param2=15, minRadius=10, maxRadius=10)

    img_gray_with_circle = img_gray.copy()
    img_gray_with_circle_bgr = cv2.cvtColor(img_gray_with_circle, cv2.COLOR_GRAY2BGR)
    img_gray_without_circle = img_gray.copy()
    img_gray_without_circle_bgr = cv2.cvtColor(img_gray_with_circle, cv2.COLOR_GRAY2BGR)    

    DICY = {}
    DICX = {}

    bb_x0 = 0
    bb_x1 = 0
    bb_y0 = 0
    bb_y1 = 0

    y_num,x_num = img_gray.shape

    bin_size = 30        
    y_num_bin0 = int(0.25*y_num/bin_size) # upperboud top
    x_num_bin0 = int(0.25*x_num/bin_size) # upperbound left
    y_num_bin1 = int(0.75*y_num/bin_size) # lowerboud bot
    x_num_bin1 = int(0.75*x_num/bin_size) # lowerbound right 
    bin_num_x = int(x_num/bin_size) # bin width
    bin_num_y = int(y_num/bin_size) # bin width

    # grid lines
    color_grid = 128
    for n in range(1, bin_num_x):
        xx = bin_size * n
        cv2.line(img_gray_with_circle_bgr, (xx, trim_dd), (xx, y_num-trim_dd), (color_grid, color_grid, color_grid), thickness=1)
    for n in range(1, bin_num_y):        
        yy = bin_size * n
        cv2.line(img_gray_with_circle_bgr, (trim_dd, yy), (x_num-trim_dd, yy), (color_grid, color_grid, color_grid), thickness=1)
    
    num_th_x0 = 10
    num_th_x1 = 10    
    num_th_y0 = 10
    num_th_y1 = 10
    
    if sid == "FFPE_LUAD_2_C":
        num_th_x0 = 12
        num_th_x1 = 12

    num_th = 10
    if circles is not None:
        for cx, cy, r, in circles.squeeze(axis=0).astype(int):
            cx_bin = int(cx/bin_size)
            cy_bin = int(cy/bin_size)
            if not cx_bin in DICX:
                DICX[cx_bin] = 0
            if not cy_bin in DICY:
                DICY[cy_bin] = 0
            DICX[cx_bin] += 1
            DICY[cy_bin] += 1
        # bb_x0
        for x,v in sorted(DICX.items()):
            if x < x_num_bin0:
                if v > num_th_x0:
                    bb_x0 = x
        # bb_x1
        for x,v in sorted(DICX.items(), reverse=True):
            if x > x_num_bin1:
                if v > num_th_x1:
                    bb_x1 = x
        # bb_y0
        for y,v in sorted(DICY.items()):
            if y < y_num_bin0:
                if v > num_th_y0:
                    bb_y0 = y
        # bb_y1
        for y,v in sorted(DICY.items(), reverse=True):
            if y > y_num_bin1:
                if v > num_th_y1:
                    bb_y1 = y
                    
    bb_x0 *= bin_size
    bb_x1 *= bin_size
    bb_y0 *= bin_size
    bb_y1 *= bin_size

    bb_x0 += bin_size
    bb_y0 += bin_size

    # lines in orange
    cv2.line(img_gray_with_circle_bgr, (bb_x0, trim_dd), (bb_x0, y_num-trim_dd), (0, 165, 255), thickness=2)
    cv2.line(img_gray_with_circle_bgr, (bb_x1, trim_dd), (bb_x1, y_num-trim_dd), (0, 165, 255), thickness=2)    

    cv2.line(img_gray_with_circle_bgr, (trim_dd, bb_y0), (x_num-trim_dd, bb_y0), (0, 165, 255), thickness=2)
    cv2.line(img_gray_with_circle_bgr, (trim_dd, bb_y1), (x_num-trim_dd, bb_y1), (0, 165, 255), thickness=2)

    r_white = 10
        
    binary_th = 160
    if sid == "FFPE_LUAD_2_C":
        binary_th = 80
        
    filename_fiducial_txt = MASKDIR_DAT + "/" + "fiducial_positions_list.txt"
    fp_fid = open(filename_fiducial_txt, "wt")

    if circles is not None:
        for cx, cy, r, in circles.squeeze(axis=0).astype(int):

            dr = 0

            cx_bin = int(cx/bin_size)
            cy_bin = int(cy/bin_size)

            cx_bin_freq = DICX[cx_bin]
            cy_bin_freq = DICY[cy_bin]

            flag = 1 # red
            if cx-r > bb_x0 and cx+dr < bb_x1 and cy-dr > bb_y0 and cy+dr < bb_y1: # center area
                flag = 0 # blue
            else:
                if cx < bb_x0: # left
                    if cy-dr < bb_y0 or cy+dr > bb_y1: # top and bot
                        if cx_bin_freq > num_th and cy_bin_freq > num_th:
                            flag = 2 # green
                        else:
                            flag = 4 # dark green
                    else:
                        if cx_bin_freq < num_th:
                            flag = 3 # magenta
                elif cx > bb_x1: # right
                    if cy-dr < bb_y0 or cy+dr > bb_y1: # top and bot
                        if cx_bin_freq > num_th and cy_bin_freq > num_th:
                            flag = 2 # green
                        else:
                            flag = 4 # dark green                            
                    elif cx_bin_freq < num_th:
                        flag = 3 # magenta
                        
            if flag == 2: # green
                print("{:d}\t{:d}".format(cx, cy), file=fp_fid)
                cv2.circle(img_gray_with_circle_bgr, (cx, cy), r, (0, 255, 0), 2)
                cv2.circle(img_gray_without_circle_bgr, (cx, cy), r_white, (255, 255, 255), thickness=-1)
                cv2.circle(img_gray_without_circle, (cx, cy), r_white, 255, thickness=-1)
            elif flag == 1: # red 
                print("{:d}\t{:d}".format(cx, cy), file=fp_fid)                
                cv2.circle(img_gray_with_circle_bgr, (cx, cy), r, (0, 0, 255), 2)
                cv2.circle(img_gray_without_circle_bgr, (cx, cy), r_white, (255, 255, 255), thickness=-1)
                cv2.circle(img_gray_without_circle, (cx, cy), r_white, 255, thickness=-1)
            elif flag == 0: # blue
                cv2.circle(img_gray_with_circle_bgr, (cx, cy), r, (255, 0, 0), 2)
            elif flag == 3: # magenta
                cv2.circle(img_gray_with_circle_bgr, (cx, cy), r, (255, 0, 255), 2)
                cv2.circle(img_gray_without_circle_bgr, (cx, cy), r_white, (255, 255, 255), thickness=-1)
                cv2.circle(img_gray_without_circle, (cx, cy), r_white, 255, thickness=-1)
            elif flag == 4: # dark green
                cv2.circle(img_gray_with_circle_bgr, (cx, cy), r, (0, 127, 0), 2)
                cv2.circle(img_gray_without_circle_bgr, (cx, cy), r_white, (255, 255, 255), thickness=-1)
                cv2.circle(img_gray_without_circle, (cx, cy), r_white, 255, thickness=-1)
                
    fp_fid.close()

    filename_gray_with_circle_bgr = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_with_circle_bgr.png"        
    cv2.imwrite(filename_gray_with_circle_bgr, img_gray_with_circle_bgr)

    filename_gray_without_circle_bgr = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_bgr.png"        
    cv2.imwrite(filename_gray_without_circle_bgr, img_gray_without_circle_bgr)

    filename_gray_without_circle = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle.png"        
    cv2.imwrite(filename_gray_without_circle, img_gray_without_circle)

    # -> white background
    index = np.all(img_gray_without_circle_bgr[:,:,:] >= 255-binary_th, axis=-1) # <= binary_th
    img_gray_without_circle_bgr_wb = img_gray_without_circle_bgr
    img_gray_without_circle_bgr_wb[index] = [255,255,255]
    img_gray_without_circle_wb = cv2.cvtColor(img_gray_without_circle_bgr_wb, cv2.COLOR_BGR2GRAY)
    filename_gray_without_circle_bgr_wb = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_bgr_wb.png"        
    cv2.imwrite(filename_gray_without_circle_bgr_wb, img_gray_without_circle_bgr_wb)

    print("CIRCLE (FIDUCIAL MARKERS)")

    print("<h3>", file=fp_html)
    print("Fiducial marker detection and removal", file=fp_html)
    print("</h3>", file=fp_html)

    filename_gray_with_circle_bgr = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_with_circle_bgr.png"
    filename_gray_with_circle_bgr_rel = "dat" + "/" + PREFIX + "_tr_gray_with_circle_bgr.png"
    print("<img src=" + "\"" + filename_gray_with_circle_bgr_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    filename_gray_without_circle_bgr = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_bgr.png"
    filename_gray_without_circle_bgr_rel = "dat" + "/" + PREFIX + "_tr_gray_without_circle_bgr.png"
    print("<img src=" + "\"" + filename_gray_without_circle_bgr_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    filename_gray_without_circle_bgr_wb = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_bgr_wb.png"
    filename_gray_without_circle_bgr_wb_rel = "dat" + "/" + PREFIX + "_tr_gray_without_circle_bgr_wb.png"
    print("<img src=" + "\"" + filename_gray_without_circle_bgr_wb_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    # Inversion
    img_gray_without_circle_nega = 255 - img_gray_without_circle
    filename_gray_without_circle_nega = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega.png"        
    cv2.imwrite(filename_gray_without_circle_nega, img_gray_without_circle_nega)

    img_gray_without_circle_wb_nega = 255 - img_gray_without_circle_wb
    filename_gray_without_circle_wb_nega = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_wb_nega.png"        
    cv2.imwrite(filename_gray_without_circle_wb_nega, img_gray_without_circle_wb_nega)

    print("<h3>", file=fp_html)
    print("Inversion", file=fp_html)
    print("</h3>", file=fp_html)

    filename_gray_without_circle_nega = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega.png"
    filename_gray_without_circle_nega_rel = "dat" + "/" + PREFIX + "_tr_gray_without_circle_nega.png"
    print("<img src=" + "\"" + filename_gray_without_circle_nega_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    filename_gray_without_circle_wb_nega = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_wb_nega.png"
    filename_gray_without_circle_wb_nega_rel = "dat" + "/" + PREFIX + "_tr_gray_without_circle_wb_nega.png"
    print("<img src=" + "\"" + filename_gray_without_circle_wb_nega_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    
    # Binarization

    # histgram

    hist_full = cv2.calcHist([img_gray_without_circle_nega], [0], None, [256], [0,256])
    filename_hist = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_hist" + ".png"
    filename_hist_rel = "dat" + "/" + PREFIX + "_tr_gray_without_circle_nega_hist" + ".png"
    
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['mathtext.rm'] = 'Arial'

    fig, ax = plt.subplots(figsize=(5,5))

    xs_full = []
    ys_full = []
    for x,d in enumerate(hist_full):
        y = d[0]
        if y > 0:
            xs_full.append(x)
            ys_full.append(y)
    
    ax.set_xlim([-3,258])
    ax.set_xticks([0, 32, 64, 96, 128, 160, 192, 224, 255])
    ax.set_yscale("log")            
    ax.set_ylim(0.8, 3e8)
    ax.set_yticks([1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8])
    ax.set_yticklabels(["1", "10", "$10^2$", "$10^3$", "$10^4$", "$10^5$", "$10^6$", "$10^7$", "$10^8$"])
    ax.set_xlabel("Pixel value")
    ax.set_ylabel("Frequency")
    ax.set_axisbelow(True)            
    ax.grid(which="major")
    ax.set_title("{}".format(sid), loc="left")
    
    ax.scatter(xs_full, ys_full, color="red", s=6, clip_on=False, zorder=1)
    fig.savefig(filename_hist, bbox_inches="tight", pad_inches=0.05, dpi=600, transparent=True)

    plt.close()
    
    # grayscale
    bin_th, img_gray_without_circle_nega_bin = cv2.threshold(img_gray_without_circle_nega, binary_th, 255, cv2.THRESH_BINARY) # GRAYSCALE -> BINARY
    filename_gray_without_circle_nega_bin = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_bin.png"        
    cv2.imwrite(filename_gray_without_circle_nega_bin, img_gray_without_circle_nega_bin)
    # color
    img_gray_without_circle_nega_bin_bgr = cv2.cvtColor(img_gray_without_circle_nega_bin, cv2.COLOR_GRAY2BGR)
    img_gray_without_circle_nega_bin_bgr[img_gray_without_circle_nega_bin == 255] = orangered 
    img_gray_without_circle_nega_bin_bgr[img_gray_without_circle_nega_bin == 0] = white
    img_gray_without_circle_nega_bin_bgr = bg_transp(img_gray_without_circle_nega_bin_bgr)
    
    filename_gray_without_circle_nega_bin_bgr = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_bin_bgr.png"                
    cv2.imwrite(filename_gray_without_circle_nega_bin_bgr, img_gray_without_circle_nega_bin_bgr)

    print("<h3>", file=fp_html)
    print("Binarization(>{:d})".format(binary_th), file=fp_html)
    print("</h3>", file=fp_html)

    print("<img src=" + "\"" + filename_hist_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)    
    
    filename_gray_without_circle_nega_bin = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_bin.png"
    filename_gray_without_circle_nega_bin_rel = "dat" + "/" + PREFIX + "_tr_gray_without_circle_nega_bin.png"
    print("<img src=" + "\"" + filename_gray_without_circle_nega_bin_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    
    filename_gray_without_circle_nega_bin_bgr = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_bin_bgr.png"
    filename_gray_without_circle_nega_bin_bgr_rel = "dat" + "/" + PREFIX + "_tr_gray_without_circle_nega_bin_bgr.png"
    print("<img src=" + "\"" + filename_gray_without_circle_nega_bin_bgr_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    # Dilation / G ... increases white area
    dy = 2
    dx = 2
    kernel = np.ones((dy,dx), np.uint8)
    
    img_gray_without_circle_nega_bin_dilate = cv2.dilate(img_gray_without_circle_nega_bin, kernel, iterations=2)
    filename_gray_without_circle_nega_bin_dilate = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_bin_dilate" + ".png"        
    cv2.imwrite(filename_gray_without_circle_nega_bin_dilate, img_gray_without_circle_nega_bin_dilate)

    # Erosion / G ... decreases white area
    dy = 2
    dx = 2
    kernel = np.ones((dy,dx), np.uint8)
    
    img_gray_without_circle_nega_bin_dilate_erode = cv2.erode(img_gray_without_circle_nega_bin_dilate, kernel, iterations=2)
    filename_gray_without_circle_nega_bin_dilate_erode = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_bin_dilate_erode" + ".png"        
    cv2.imwrite(filename_gray_without_circle_nega_bin_dilate_erode, img_gray_without_circle_nega_bin_dilate_erode)

    # Dilation / G ... increases white area
    
    dy = 2
    dx = 2
    kernel = np.ones((dy,dx), np.uint8)
    
    img_gray_without_circle_nega_bin_dilate_erode_dilate = cv2.dilate(img_gray_without_circle_nega_bin_dilate_erode, kernel, iterations=2)
    filename_gray_without_circle_nega_bin_dilate_erode_dilate = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_bin_dilate_erode_dilate" + ".png"        
    cv2.imwrite(filename_gray_without_circle_nega_bin_dilate_erode_dilate, img_gray_without_circle_nega_bin_dilate_erode_dilate)

    print("<h3>", file=fp_html)
    print("Dilation(2x2, x2), erosion(2x2, x2) and dilation(2x2, x2) [OPTIONAL]", file=fp_html)
    print("</h3>", file=fp_html)
    
    filename_gray_without_circle_nega_bin_dilate = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_bin_dilate.png"
    filename_gray_without_circle_nega_bin_dilate_rel = "dat" + "/" + PREFIX + "_tr_gray_without_circle_nega_bin_dilate.png"
    print("<img src=" + "\"" + filename_gray_without_circle_nega_bin_dilate_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    filename_gray_without_circle_nega_bin_dilate_erode = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_bin_dilate_erode.png"
    filename_gray_without_circle_nega_bin_dilate_erode_rel = "dat" + "/" + PREFIX + "_tr_gray_without_circle_nega_bin_dilate_erode.png"
    print("<img src=" + "\"" + filename_gray_without_circle_nega_bin_dilate_erode_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    filename_gray_without_circle_nega_bin_dilate_erode_dilate = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_bin_dilate_erode_dilate.png"
    filename_gray_without_circle_nega_bin_dilate_erode_dilate_rel = "dat" + "/" + PREFIX + "_tr_gray_without_circle_nega_bin_dilate_erode_dilate.png"
    print("<img src=" + "\"" + filename_gray_without_circle_nega_bin_dilate_erode_dilate_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    
    # Contours

    contours, hierarchy = cv2.findContours(img_gray_without_circle_nega_bin, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)    
    contour_max = max(contours, key=lambda x: cv2.contourArea(x))

    # holes
    hole_th = 800
    tiss_th = 100

    contours_tiss = []
    contours_hole = []
    for cnt in contours:
        area = cv2.contourArea(cnt, True)
        if area > 0 and abs(area) > hole_th:
            contours_hole.append(cnt)
        elif area < 0 and abs(area) > tiss_th:
            contours_tiss.append(cnt)
            
    # Region

    img_gray_without_circle_nega_bin_bgr_cont = cv2.cvtColor(img_gray_without_circle_nega_bin, cv2.COLOR_GRAY2BGR)
    cv2.drawContours(img_gray_without_circle_nega_bin_bgr_cont, contours_tiss, -1, color=(0,0,255), thickness=2)    
    cv2.drawContours(img_gray_without_circle_nega_bin_bgr_cont, contours_hole, -1, color=(255,0,0), thickness=2) # holes
    filename_gray_without_circle_nega_bin_bgr_cont = MASKDIR_DAT + "/" + PREFIX + "_" + "tr_gray_without_circle_nega_bin_bgr_cont" + ".png"    
    cv2.imwrite(filename_gray_without_circle_nega_bin_bgr_cont, img_gray_without_circle_nega_bin_bgr_cont)

    img_gray_without_circle_nega_bgr_cont = cv2.cvtColor(img_gray_without_circle_nega, cv2.COLOR_GRAY2BGR)
    cv2.drawContours(img_gray_without_circle_nega_bgr_cont, contours_tiss, -1, color=(0,0,255), thickness=2)    
    cv2.drawContours(img_gray_without_circle_nega_bgr_cont, contours_hole, -1, color=(255,0,0), thickness=2) # holes
    filename_gray_without_circle_nega_bgr_cont = MASKDIR_DAT + "/" + PREFIX + "_" + "tr_gray_without_circle_nega_bgr_cont" + ".png"    
    cv2.imwrite(filename_gray_without_circle_nega_bgr_cont, img_gray_without_circle_nega_bgr_cont)

    img_gray_without_circle_bgr_cont = cv2.cvtColor(img_gray_without_circle, cv2.COLOR_GRAY2BGR)
    cv2.drawContours(img_gray_without_circle_bgr_cont, contours_tiss, -1, color=(0,0,255), thickness=2)    
    cv2.drawContours(img_gray_without_circle_bgr_cont, contours_hole, -1, color=(255,0,0), thickness=2) # holes
    filename_gray_without_circle_bgr_cont = MASKDIR_DAT + "/" + PREFIX + "_" + "tr_gray_without_circle_bgr_cont" + ".png"    
    cv2.imwrite(filename_gray_without_circle_bgr_cont, img_gray_without_circle_bgr_cont)

    img_gray_without_circle_bgr_wb_cont = img_gray_without_circle_bgr_wb.copy()
    cv2.drawContours(img_gray_without_circle_bgr_wb_cont, contours_tiss, -1, color=(0,0,255), thickness=2)    
    cv2.drawContours(img_gray_without_circle_bgr_wb_cont, contours_hole, -1, color=(255,0,0), thickness=2) # holes
    img_gray_without_circle_bgr_wb_cont = bg_transp(img_gray_without_circle_bgr_wb_cont)
    filename_gray_without_circle_bgr_wb_cont = MASKDIR_DAT + "/" + PREFIX + "_" + "tr_gray_without_circle_bgr_wb_cont" + ".png"    
    cv2.imwrite(filename_gray_without_circle_bgr_wb_cont, img_gray_without_circle_bgr_wb_cont)

    print("<h3>", file=fp_html)
    print("Region detection", file=fp_html)
    print("</h3>", file=fp_html)

    filename_gray_without_circle_nega_bin_bgr_cont = MASKDIR_DAT + "/" + PREFIX + "_" + "tr_gray_without_circle_nega_bin_bgr_cont" + ".png"
    filename_gray_without_circle_nega_bin_bgr_cont_rel = "dat" + "/" + PREFIX + "_" + "tr_gray_without_circle_nega_bin_bgr_cont" + ".png"        
    print("<img src=" + "\"" + filename_gray_without_circle_nega_bin_bgr_cont_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    filename_gray_without_circle_nega_bgr_cont = MASKDIR_DAT + "/" + PREFIX + "_" + "tr_gray_without_circle_nega_bgr_cont" + ".png"
    filename_gray_without_circle_nega_bgr_cont_rel = "dat" + "/" + PREFIX + "_" + "tr_gray_without_circle_nega_bgr_cont" + ".png"        
    print("<img src=" + "\"" + filename_gray_without_circle_nega_bgr_cont_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    filename_gray_without_circle_bgr_cont = MASKDIR_DAT + "/" + PREFIX + "_" + "tr_gray_without_circle_bgr_cont" + ".png"
    filename_gray_without_circle_bgr_cont_rel = "dat" + "/" + PREFIX + "_" + "tr_gray_without_circle_bgr_cont" + ".png"        
    print("<img src=" + "\"" + filename_gray_without_circle_bgr_cont_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    filename_gray_without_circle_bgr_wb_cont = MASKDIR_DAT + "/" + PREFIX + "_" + "tr_gray_without_circle_bgr_wb_cont" + ".png"
    filename_gray_without_circle_bgr_wb_cont_rel = "dat" + "/" + PREFIX + "_" + "tr_gray_without_circle_bgr_wb_cont" + ".png"        
    print("<img src=" + "\"" + filename_gray_without_circle_bgr_wb_cont_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    

    # Mask generation
    
    img_gray_without_circle_nega_mask = np.zeros_like(img_gray)
    for k, v in enumerate(hierarchy[0]):
        next, prev, child, parent = v
        cnt = contours[k]
        area = cv2.contourArea(cnt, True)
        if parent == -1:
            if abs(area) > tiss_th:
                cv2.drawContours(img_gray_without_circle_nega_mask, [cnt], -1, color=255, thickness=-1)

    filename_gray_without_circle_nega_mask = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_mask.png"    
    cv2.imwrite(filename_gray_without_circle_nega_mask, img_gray_without_circle_nega_mask)

    img_gray_without_circle_nega_mask_bgr = cv2.cvtColor(img_gray_without_circle_nega_mask, cv2.COLOR_GRAY2BGR)
    img_gray_without_circle_nega_mask_bgr[img_gray_without_circle_nega_mask == 255] = orangered
    img_gray_without_circle_nega_mask_bgr[img_gray_without_circle_nega_mask == 0] = white
    img_gray_without_circle_nega_mask_bgr = bg_transp(img_gray_without_circle_nega_mask_bgr)
    filename_gray_without_circle_nega_mask_bgr = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_mask_bgr.png"                
    cv2.imwrite(filename_gray_without_circle_nega_mask_bgr, img_gray_without_circle_nega_mask_bgr)


    img_gray_without_circle_nega_mask_with_hole = np.zeros_like(img_gray)
    for k, v in enumerate(hierarchy[0]):
        next, prev, child, parent = v
        cnt = contours[k]
        area = cv2.contourArea(cnt, True)
        if area < 0 and abs(area) > tiss_th:
            cv2.drawContours(img_gray_without_circle_nega_mask_with_hole, [cnt], -1, color=255, thickness=-1)
        elif area > 0 and abs(area) > hole_th:
            cv2.drawContours(img_gray_without_circle_nega_mask_with_hole, [cnt], -1, color=0, thickness=-1)                    

    filename_gray_without_circle_nega_mask_with_hole = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_mask_with_hole.png"                
    cv2.imwrite(filename_gray_without_circle_nega_mask_with_hole, img_gray_without_circle_nega_mask_with_hole)

    img_gray_without_circle_nega_mask_with_hole_bgr = cv2.cvtColor(img_gray_without_circle_nega_mask_with_hole, cv2.COLOR_GRAY2BGR)
    img_gray_without_circle_nega_mask_with_hole_bgr[img_gray_without_circle_nega_mask_with_hole == 255] = orangered
    img_gray_without_circle_nega_mask_with_hole_bgr[img_gray_without_circle_nega_mask_with_hole == 0] = white
    img_gray_without_circle_nega_mask_with_hole_bgr = bg_transp(img_gray_without_circle_nega_mask_with_hole_bgr)
    filename_gray_without_circle_nega_mask_with_hole_bgr = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_mask_with_hole_bgr.png"                
    cv2.imwrite(filename_gray_without_circle_nega_mask_with_hole_bgr, img_gray_without_circle_nega_mask_with_hole_bgr)
    
    print("<h3>", file=fp_html)
    print("Mask generation", file=fp_html)
    print("</h3>", file=fp_html)

    filename_gray_without_circle_nega_mask = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_mask.png"
    filename_gray_without_circle_nega_mask_rel = "dat" + "/" + PREFIX + "_tr_gray_without_circle_nega_mask.png"
    print("<img src=" + "\"" + filename_gray_without_circle_nega_mask_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    filename_gray_without_circle_nega_mask_with_hole = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_mask_with_hole.png"            
    filename_gray_without_circle_nega_mask_with_hole_rel = "dat" + "/" + PREFIX + "_tr_gray_without_circle_nega_mask_with_hole.png"
    print("<img src=" + "\"" + filename_gray_without_circle_nega_mask_with_hole_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    filename_gray_without_circle_nega_mask_bgr = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_mask_bgr.png"            
    filename_gray_without_circle_nega_mask_bgr_rel = "dat" + "/" + PREFIX + "_tr_gray_without_circle_nega_mask_bgr.png"
    print("<img src=" + "\"" + filename_gray_without_circle_nega_mask_bgr_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    filename_gray_without_circle_nega_mask_with_hole_bgr = MASKDIR_DAT + "/" + PREFIX + "_tr_gray_without_circle_nega_mask_with_hole_bgr.png"            
    filename_gray_without_circle_nega_mask_with_hole_bgr_rel = "dat" + "/" + PREFIX + "_tr_gray_without_circle_nega_mask_with_hole_bgr.png"
    print("<img src=" + "\"" + filename_gray_without_circle_nega_mask_with_hole_bgr_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    
    print_html_footer(fp_html)
    
    fp_html.close()

    
#####

if __name__ == "__main__":
    main()
