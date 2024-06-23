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
from matplotlib.colors import LinearSegmentedColormap

from xml.etree import ElementTree as ET
from tifffile import TiffFile
import tifffile

from util.image import *
from util.html import *

def main():

    binary_th = 10

    if len(sys.argv) < 3:
        print("ERROR: Command line arguments")
        sys.exit()

    # path to QPTIFF file
    filename_qptiff = sys.argv[1]
    # sample name
    sid = sys.argv[2]

    pp = pathlib.Path(filename_qptiff)
    pp_name = pp.name
    pp_stem = pp.stem

    print(filename_qptiff)
    print(pp_name, pp_stem)

    # directories    
    WORKDIR = "."
    RESULTDIR = WORKDIR + "/" + sid
    pp = pathlib.Path(RESULTDIR)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)
    PREFIX = sid + "-" + "PhenoCycler"
    MASKDIR = RESULTDIR + "/" + PREFIX + "-" + "MASK"
    pp = pathlib.Path(MASKDIR)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)
    MASKDIR_DAT = MASKDIR + "/" + "dat"
    pp = pathlib.Path(MASKDIR_DAT)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)
    PAGESDIR = MASKDIR_DAT + "/" +  "PAGES"
    pp = pathlib.Path(PAGESDIR)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)
    MULTIDIR = MASKDIR_DAT + "/" + "MULTI"    
    pp = pathlib.Path(MULTIDIR)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)

    # load images 
    filename_pd = PAGESDIR + "/" + "page_description.txt"
    filename_reso = PAGESDIR + "/" + "resolution_nm.txt"

    fp_pd = open(filename_pd, "wt")
    fp_reso = open(filename_reso, "wt")

    IMAGES = []
    with TiffFile(filename_qptiff) as tif:

        count = 0
        for page in tif.series[0].pages:
            for tag in page.tags.values():
                name, value = tag.name, tag.value
                value_ss = str(value)
                print("{}\t{}".format(name, value_ss), file=fp_pd)

            root = ET.fromstring(page.description)
            name_org = root.find("Biomarker").text            
            name = re.sub("/+", "-", name_org)
            name = re.sub("-+", "-", name)

            if count == 0:
                resolution_nm = root.find("ScanProfile").find("ExperimentV4").find("Resolution_nm")
                print("resolution_nm", resolution_nm.text)
                print("resolution_nm" + "\t" + resolution_nm.text, file=fp_reso)

            img = tifffile.imread(filename_qptiff, key=count)
            IMAGES.append((name, name_org, img))

            count += 1


    fp_pd.close()
    fp_reso.close()
            
    images_num = len(IMAGES)

    # color map
    cm_name = "cm_custom"
    colors = ["#0020ff", "#0040ff", "#0080ff", "#00c0ff", "#00e0ff", #"#00ffff", 
              "#00ffe0", "#00ffc0", "#00ff80", "#00ff40", "#00ff20", #"#00ff00",   
              "#20ff00", "#40ff00", "#80ff00", "#c0ff00", "#e0ff00", #"#ffff00",
              "#ffe000", "#ffc000", "#ff8000", "#ff4000", "#ff2000"]
    cm = LinearSegmentedColormap.from_list(cm_name, colors, N=images_num)

    # page IDs for MULTI
    #PASSED = range(0, len(IMAGES))
    PASSED = [22, 4, 0, 1, 2, 10, 12, 21]

    # html file
    filename_html = MASKDIR + "/" + PREFIX + ".html"
    print(filename_html)
    fp_html = open(filename_html, "wt")

    ww = "200px"
    print_html_header(fp_html)
    
    print("<h1>", file=fp_html)
    print("PhenoCycler QPTIFF" + " &mdash; " + sid, file=fp_html)
    print("</h1>", file=fp_html)

    # PASSED, RGB
    filename_passed = PAGESDIR + "/" + PREFIX + "-PASSED.txt"
    filename_passed_rel = "dat" + "/" + "PAGES" + "/" + PREFIX + "-PASSED.txt"        

    filename_rgb = PAGESDIR + "/" + PREFIX + "-RGB.txt"        
    filename_rgb_rel = "dat" + "/" + "PAGES" + "/" + PREFIX + "-RGB.txt"

    # MULTI pages
    
    filename_MULTI_gray_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_01" + ".png"                
    filename_MULTI_gray_nega_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_nega_01" + ".png"                
    filename_MULTI_bgr_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_bgr_bb_01" + ".png"            
    filename_MULTI_bgr_wb_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_bgr_wb_01" + ".png"            

    print("<h2>", file=fp_html)
    print("Multi pages superimposed (" + "<a href=\"" + filename_passed_rel + "\">" + "antibodies" + "</a>" + ")", file=fp_html)
    print("</h2>", file=fp_html)

    print("<img src=" + "\"" + filename_MULTI_bgr_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)        
    print("<img src=" + "\"" + filename_MULTI_bgr_wb_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    print("<img src=" + "\"" + filename_MULTI_gray_nega_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)                
    print("<img src=" + "\"" + filename_MULTI_gray_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    fp_html.flush()    

    
    # MULTI
    # MULTI BGR BB (black back)
    
    filename_passed = PAGESDIR + "/" + PREFIX + "-PASSED.txt"
    fp_passed = open(filename_passed, "wt")

    print("Page\tName\tName(original)\t", fp_passed)

    count = 0
    for ci in PASSED:

        name, name_org, img = IMAGES[ci]  # black background

        print("{:02d}".format(ci+1) + "\t" + name + "\t" + name_org, file=fp_passed)

        cc = [int(v * 255) for v in reversed(cm(ci)[:3])] # RGB -> BGR
        [cc_b, cc_g, cc_r] = cc
        
        print("CC", ci, cc)

        img_th = np.where(img <= binary_th, 0, img)
        img_th_bgr = cv2.cvtColor(img_th, cv2.COLOR_GRAY2BGR)
        
        if count == 0:

            img_MULTI_bgr_bb = img_th_bgr.copy()            
            img_max = img_th.copy()            
            img_MULTI_bgr_bb[img_th > binary_th] = cc
            
        img_MULTI_bgr_bb[img_th > img_max] = cc
        img_max = np.where(img_th > img_max, img_th, img_max)                

        count += 1

    fp_passed.close()

    filename_MULTI_bgr_bb = MULTIDIR + "/" + PREFIX + "-" + "MULTI_bgr_bb" + ".png"    
    cv2.imwrite(filename_MULTI_bgr_bb, img_MULTI_bgr_bb)

    img_MULTI_bgr_bb_01 = cv2.resize(img_MULTI_bgr_bb, None, None, 0.1, 0.1, cv2.INTER_NEAREST)
    filename_MULTI_bgr_bb_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_bgr_bb_01" + ".png"    
    cv2.imwrite(filename_MULTI_bgr_bb_01, img_MULTI_bgr_bb_01)

    # white background (BB -> WB)
    img_MULTI_bgr_bb_gray = cv2.cvtColor(img_MULTI_bgr_bb, cv2.COLOR_BGR2GRAY) # BB
    img_MULTI_bgr_wb = img_MULTI_bgr_bb.copy()
    img_MULTI_bgr_wb[img_MULTI_bgr_bb_gray <= binary_th] = [255,255,255]
    
    filename_MULTI_bgr_wb = MULTIDIR + "/" + PREFIX + "-" + "MULTI_bgr_wb" + ".png"    
    cv2.imwrite(filename_MULTI_bgr_wb, img_MULTI_bgr_wb)

    img_MULTI_bgr_wb_01 = cv2.resize(img_MULTI_bgr_wb, None, None, 0.1, 0.1, cv2.INTER_NEAREST)
    img_MULTI_bgr_wb_01 = bg_transp(img_MULTI_bgr_wb_01)    
    filename_MULTI_bgr_wb_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_bgr_wb_01" + ".png"    
    cv2.imwrite(filename_MULTI_bgr_wb_01, img_MULTI_bgr_wb_01)


    ########

    # Mask by MULTI

    img_MULTI_gray_nega = img_MULTI_bgr_bb_gray # BB // the original

    # MULTI GRAY NEGA
    filename_MULTI_gray_nega = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega" + ".png"
    cv2.imwrite(filename_MULTI_gray_nega, img_MULTI_gray_nega)

    img_MULTI_gray_nega_01 = cv2.resize(img_MULTI_gray_nega, None, None, 0.1, 0.1, cv2.INTER_NEAREST)
    filename_MULTI_gray_nega_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_01" + ".png"
    filename_MULTI_gray_nega_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_nega_01" + ".png"                
    cv2.imwrite(filename_MULTI_gray_nega_01, img_MULTI_gray_nega_01)

    print("MULTI GRAY NEGA")

    # MULTI GRAY (POS)

    #img_MULTI_gray = cv2.bitwise_not(img_MULTI_gray_nega)  # white background
    img_MULTI_gray = 255 - img_MULTI_gray_nega  # white background
    filename_MULTI_gray = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray" + ".png"    
    cv2.imwrite(filename_MULTI_gray, img_MULTI_gray)

    img_MULTI_gray_01 = cv2.resize(img_MULTI_gray, None, None, 0.1, 0.1, cv2.INTER_NEAREST)
    img_MULTI_gray_01 = bg_transp(img_MULTI_gray_01)
    filename_MULTI_gray_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_01" + ".png"
    filename_MULTI_gray_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_01" + ".png"            
    cv2.imwrite(filename_MULTI_gray_01, img_MULTI_gray_01)

    print("MULTI GRAY (POS)")
    #print(img_MULTI_gray)

    # histgram
    hist_full = cv2.calcHist([img_MULTI_gray_nega], [0], None, [256], [0,256])
    filename_hist = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_hist" + ".png"
    filename_hist_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_nega_hist" + ".png"
    
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

    # MULTI GRAY NEGA -> BIN

    bin_th, img_MULTI_gray_nega_bin = cv2.threshold(img_MULTI_gray_nega, binary_th, 255, cv2.THRESH_BINARY) # GRAYSCALE -> BINARY    

    filename_MULTI_gray_nega_bin = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_bin" + ".png"    
    cv2.imwrite(filename_MULTI_gray_nega_bin, img_MULTI_gray_nega_bin)

    img_MULTI_gray_nega_bin_01 = cv2.resize(img_MULTI_gray_nega_bin, None, None, 0.1, 0.1, cv2.INTER_NEAREST)
    filename_MULTI_gray_nega_bin_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_01" + ".png"
    filename_MULTI_gray_nega_bin_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_01" + ".png"            
    cv2.imwrite(filename_MULTI_gray_nega_bin_01, img_MULTI_gray_nega_bin_01)

    # MULTI GRAY NEGA -> BIN -> BGR (BLUE)

    img_MULTI_gray_bin_bgr = cv2.cvtColor(img_MULTI_gray_nega_bin, cv2.COLOR_GRAY2BGR)    
    img_MULTI_gray_bin_bgr[img_MULTI_gray_nega_bin == 255] = royalblue
    img_MULTI_gray_bin_bgr[img_MULTI_gray_nega_bin == 0] = white
    filename_MULTI_gray_bin_bgr = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_bin_bgr" + ".png"    
    cv2.imwrite(filename_MULTI_gray_bin_bgr, img_MULTI_gray_bin_bgr)

    img_MULTI_gray_bin_bgr_01 = cv2.resize(img_MULTI_gray_bin_bgr, None, None, 0.1, 0.1, cv2.INTER_NEAREST)
    img_MULTI_gray_bin_bgr_01 = bg_transp(img_MULTI_gray_bin_bgr_01)    
    filename_MULTI_gray_bin_bgr_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_bin_bgr_01" + ".png"
    filename_MULTI_gray_bin_bgr_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_bin_bgr_01" + ".png"            
    cv2.imwrite(filename_MULTI_gray_bin_bgr_01, img_MULTI_gray_bin_bgr_01)

    print("<h2>", file=fp_html)
    print("PhenoCycler masks", file=fp_html)
    print("</h2>", file=fp_html)    

    print("<h3>Binarization(>10)</h3>", file=fp_html)         
    #print("<img src=" + "\"" + filename_MULTI_gray_nega_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html) 
    #print("<img src=" + "\"" + filename_MULTI_gray_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    #print("<img src=" + "\"" + filename_hist_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html) 
    
    print("<img src=" + "\"" + filename_MULTI_gray_nega_bin_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)                
    print("<img src=" + "\"" + filename_MULTI_gray_bin_bgr_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    fp_html.flush()        

    print("DILATION/EROSION")

    # Dilation (1) ... increases white area
    dy = 2
    dx = 2
    kernel = np.ones((dy,dx), np.uint8)
    iter_num = 5

    img_MULTI_gray_nega_bin_dilate = cv2.dilate(img_MULTI_gray_nega_bin, kernel, iterations=iter_num)
    filename_MULTI_gray_nega_bin_dilate = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_dilate" + ".png"    
    cv2.imwrite(filename_MULTI_gray_nega_bin_dilate, img_MULTI_gray_nega_bin_dilate)

    img_MULTI_gray_nega_bin_dilate_01 = cv2.resize(img_MULTI_gray_nega_bin_dilate, None, None, 0.1, 0.1, cv2.INTER_NEAREST)                
    filename_MULTI_gray_nega_bin_dilate_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_dilate_01" + ".png"
    filename_MULTI_gray_nega_bin_dilate_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_dilate_01" + ".png"                
    cv2.imwrite(filename_MULTI_gray_nega_bin_dilate_01, img_MULTI_gray_nega_bin_dilate_01)

    # Erosion ... decreases white area
    dy = 2
    dx = 2
    kernel = np.ones((dy,dx), np.uint8)
    iter_num = 10

    img_MULTI_gray_nega_bin_dilate_erode = cv2.erode(img_MULTI_gray_nega_bin_dilate, kernel, iterations=iter_num)
    filename_MULTI_gray_nega_bin_dilate_erode = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_dilate_erode" + ".png"
    cv2.imwrite(filename_MULTI_gray_nega_bin_dilate_erode, img_MULTI_gray_nega_bin_dilate_erode)

    img_MULTI_gray_nega_bin_dilate_erode_01 = cv2.resize(img_MULTI_gray_nega_bin_dilate_erode, None, None, 0.1, 0.1, cv2.INTER_NEAREST)                    
    filename_MULTI_gray_nega_bin_dilate_erode_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_dilate_erode_01" + ".png"
    filename_MULTI_gray_nega_bin_dilate_erode_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_dilate_erode_01" + ".png"                
    cv2.imwrite(filename_MULTI_gray_nega_bin_dilate_erode_01, img_MULTI_gray_nega_bin_dilate_erode_01)

    # Dilation (2) ... increases white area    
    dy = 2
    dx = 2
    kernel = np.ones((dy,dx), np.uint8)
    iter_num = 10

    img_MULTI_gray_nega_bin_dilate_erode_dilate = cv2.dilate(img_MULTI_gray_nega_bin_dilate_erode, kernel, iterations=iter_num)
    filename_MULTI_gray_nega_bin_dilate_erode_dilate = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_dilate_erode_dilate" + ".png"
    cv2.imwrite(filename_MULTI_gray_nega_bin_dilate_erode_dilate, img_MULTI_gray_nega_bin_dilate_erode_dilate)

    img_MULTI_gray_nega_bin_dilate_erode_dilate_01 = cv2.resize(img_MULTI_gray_nega_bin_dilate_erode_dilate, None, None, 0.1, 0.1, cv2.INTER_NEAREST)                    
    filename_MULTI_gray_nega_bin_dilate_erode_dilate_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_dilate_erode_dilate_01" + ".png"
    filename_MULTI_gray_nega_bin_dilate_erode_dilate_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_dilate_erode_dilate_01" + ".png"    
    cv2.imwrite(filename_MULTI_gray_nega_bin_dilate_erode_dilate_01, img_MULTI_gray_nega_bin_dilate_erode_dilate_01)

    print("<h3>Dilation(2x2, x5) and erosion(2x2, x10) and dilation(2x2, x10)</h3>", file=fp_html)
    print("<img src=" + "\"" + filename_MULTI_gray_nega_bin_dilate_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)        
    print("<img src=" + "\"" + filename_MULTI_gray_nega_bin_dilate_erode_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    print("<img src=" + "\"" + filename_MULTI_gray_nega_bin_dilate_erode_dilate_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)        

    # Contours

    print("FINDCONTOURS")

    contours, hierarchy = cv2.findContours(img_MULTI_gray_nega_bin_dilate_erode_dilate, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    contour_max = max(contours, key=lambda x: cv2.contourArea(x))

    # holes

    hole_th = 10000
    tiss_th = 10000    
    
    contours_hole_multi = []
    contours_tiss_multi = []
    for cnt in contours:
        area = cv2.contourArea(cnt, True)

        if area > 0 and abs(area) > hole_th:
            contours_hole_multi.append(cnt)
        elif area < 0 and abs(area) > tiss_th:            
            contours_tiss_multi.append(cnt)    

    # POS
    img_MULTI_gray_bgr_cont = cv2.cvtColor(img_MULTI_gray, cv2.COLOR_GRAY2BGR)
    cv2.drawContours(img_MULTI_gray_bgr_cont, contours_tiss_multi, -1, color=(0,0,255), thickness=15)    
    cv2.drawContours(img_MULTI_gray_bgr_cont, contours_hole_multi, -1, color=(255,0,0), thickness=15)
    filename_MULTI_gray_bgr_cont = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_cont" + ".png"
    cv2.imwrite(filename_MULTI_gray_bgr_cont, img_MULTI_gray_bgr_cont)

    img_MULTI_gray_bgr_cont_01 = cv2.resize(img_MULTI_gray_bgr_cont, None, None, 0.1, 0.1, cv2.INTER_NEAREST)
    img_MULTI_gray_bgr_cont_01 = bg_transp(img_MULTI_gray_bgr_cont_01)
    filename_MULTI_gray_bgr_cont_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_cont_01" + ".png"
    filename_MULTI_gray_bgr_cont_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_cont_01" + ".png"                
    cv2.imwrite(filename_MULTI_gray_bgr_cont_01, img_MULTI_gray_bgr_cont_01)

    # NEGA
    img_MULTI_gray_nega_bgr_cont = cv2.cvtColor(img_MULTI_gray_nega, cv2.COLOR_GRAY2BGR)
    cv2.drawContours(img_MULTI_gray_nega_bgr_cont, contours_tiss_multi, -1, color=(0,0,255), thickness=15)    
    cv2.drawContours(img_MULTI_gray_nega_bgr_cont, contours_hole_multi, -1, color=(255,0,0), thickness=15)    
    filename_MULTI_gray_nega_bgr_cont = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_cont" + ".png"    
    cv2.imwrite(filename_MULTI_gray_nega_bgr_cont, img_MULTI_gray_nega_bgr_cont)

    img_MULTI_gray_nega_bgr_cont_01 = cv2.resize(img_MULTI_gray_nega_bgr_cont, None, None, 0.1, 0.1, cv2.INTER_NEAREST)                            
    filename_MULTI_gray_nega_bgr_cont_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_cont_01" + ".png"
    filename_MULTI_gray_nega_bgr_cont_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_nega_cont_01" + ".png"            
    cv2.imwrite(filename_MULTI_gray_nega_bgr_cont_01, img_MULTI_gray_nega_bgr_cont_01)

    # NEGA BIN
    img_MULTI_gray_nega_bin_bgr_cont = cv2.cvtColor(img_MULTI_gray_nega_bin, cv2.COLOR_GRAY2BGR)
    cv2.drawContours(img_MULTI_gray_nega_bin_bgr_cont, contours_tiss_multi, -1, color=(0,0,255), thickness=15)    
    cv2.drawContours(img_MULTI_gray_nega_bin_bgr_cont, contours_hole_multi, -1, color=(255,0,0), thickness=15)    
    filename_MULTI_gray_nega_bin_bgr_cont = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_cont" + ".png"    
    cv2.imwrite(filename_MULTI_gray_nega_bin_bgr_cont, img_MULTI_gray_nega_bin_bgr_cont)

    img_MULTI_gray_nega_bin_bgr_cont_01 = cv2.resize(img_MULTI_gray_nega_bin_bgr_cont, None, None, 0.1, 0.1, cv2.INTER_NEAREST)                                
    filename_MULTI_gray_nega_bin_bgr_cont_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_cont_01" + ".png"
    filename_MULTI_gray_nega_bin_bgr_cont_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_cont_01" + ".png"            
    cv2.imwrite(filename_MULTI_gray_nega_bin_bgr_cont_01, img_MULTI_gray_nega_bin_bgr_cont_01)

    # NEGA BIN (boundary only)
    img_MULTI_gray_nega_bin_cont_boundary = np.zeros_like(img_MULTI_gray_nega_bin)

    for k, v in enumerate(hierarchy[0]):
        next, prev, child, parent = v
        cnt = contours[k]
        area = cv2.contourArea(cnt, True)
        if parent == -1:
            if abs(area) > tiss_th:
                cv2.drawContours(img_MULTI_gray_nega_bin_cont_boundary, [cnt], -1, color=255, thickness=15)
                
    filename_MULTI_gray_nega_bin_cont_boundary = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_cont_boundary" + ".png"    
    cv2.imwrite(filename_MULTI_gray_nega_bin_cont_boundary, img_MULTI_gray_nega_bin_cont_boundary)

    img_MULTI_gray_nega_bin_cont_boundary_01 = cv2.resize(img_MULTI_gray_nega_bin_cont_boundary, None, None, 0.1, 0.1, cv2.INTER_NEAREST)                                
    filename_MULTI_gray_nega_bin_cont_boundary_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_cont_boundary_01" + ".png"
    filename_MULTI_gray_nega_bin_cont_boundary_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_cont_boundary_01" + ".png"            
    cv2.imwrite(filename_MULTI_gray_nega_bin_cont_boundary_01, img_MULTI_gray_nega_bin_cont_boundary_01)

    # NEGA BIN (boundary only) with hole
    img_MULTI_gray_nega_bin_cont_boundary_with_hole = img_MULTI_gray_nega_bin_cont_boundary.copy()
    cv2.drawContours(img_MULTI_gray_nega_bin_cont_boundary_with_hole, contours_tiss_multi, -1, 255, thickness=15)    
    cv2.drawContours(img_MULTI_gray_nega_bin_cont_boundary_with_hole, contours_hole_multi, -1, 255, thickness=15)
    filename_MULTI_gray_nega_bin_cont_boundary_with_hole = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_cont_boundary_with_hole" + ".png"    
    cv2.imwrite(filename_MULTI_gray_nega_bin_cont_boundary_with_hole, img_MULTI_gray_nega_bin_cont_boundary_with_hole)

    img_MULTI_gray_nega_bin_cont_boundary_with_hole_01 = cv2.resize(img_MULTI_gray_nega_bin_cont_boundary_with_hole, None, None, 0.1, 0.1, cv2.INTER_NEAREST)
    filename_MULTI_gray_nega_bin_cont_boundary_with_hole_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_cont_boundary_with_hole_01" + ".png"
    filename_MULTI_gray_nega_bin_cont_boundary_with_hole_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_cont_boundary_with_hole_01" + ".png" 
    cv2.imwrite(filename_MULTI_gray_nega_bin_cont_boundary_with_hole_01, img_MULTI_gray_nega_bin_cont_boundary_with_hole_01)

    # NEGA BIN D/E/D
    img_MULTI_gray_nega_bin_dilate_erode_dilate_bgr_cont = cv2.cvtColor(img_MULTI_gray_nega_bin_dilate_erode_dilate, cv2.COLOR_GRAY2BGR)
    cv2.drawContours(img_MULTI_gray_nega_bin_dilate_erode_dilate_bgr_cont, contours_tiss_multi, -1, color=(0,0,255), thickness=15)    
    cv2.drawContours(img_MULTI_gray_nega_bin_dilate_erode_dilate_bgr_cont, contours_hole_multi, -1, color=(255,0,0), thickness=15)    
    filename_MULTI_gray_nega_bin_dilate_erode_dilate_bgr_cont = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_dilate_erode_dilate_cont" + ".png"    
    cv2.imwrite(filename_MULTI_gray_nega_bin_dilate_erode_dilate_bgr_cont, img_MULTI_gray_nega_bin_dilate_erode_dilate_bgr_cont)

    img_MULTI_gray_nega_bin_dilate_erode_dilate_bgr_cont_01 = cv2.resize(img_MULTI_gray_nega_bin_dilate_erode_dilate_bgr_cont, None, None, 0.1, 0.1, cv2.INTER_NEAREST)
    filename_MULTI_gray_nega_bin_dilate_erode_dilate_bgr_cont_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_dilate_erode_dilate_cont_01" + ".png"
    filename_MULTI_gray_nega_bin_dilate_erode_dilate_bgr_cont_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_nega_bin_dilate_erode_dilate_cont_01" + ".png" 
    cv2.imwrite(filename_MULTI_gray_nega_bin_dilate_erode_dilate_bgr_cont_01, img_MULTI_gray_nega_bin_dilate_erode_dilate_bgr_cont_01)

    print("<h3>Region detection</h3>", file=fp_html)
    print("<img src=" + "\"" + filename_MULTI_gray_nega_bin_dilate_erode_dilate_bgr_cont_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    
    print("<img src=" + "\"" + filename_MULTI_gray_nega_bgr_cont_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)        
    print("<img src=" + "\"" + filename_MULTI_gray_nega_bin_bgr_cont_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    print("<img src=" + "\"" + filename_MULTI_gray_bgr_cont_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)    
    #print("<img src=" + "\"" + filename_MULTI_gray_nega_bin_cont_boundary_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    #print("<img src=" + "\"" + filename_MULTI_gray_nega_bin_cont_boundary_with_hole_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)        

    # MASK (NEGA)

    # solid

    print("MASK NEGA SOLID")
    img_MULTI_gray_nega_mask = np.zeros_like(img)

    for k, v in enumerate(hierarchy[0]):
        next, prev, child, parent = v
        cnt = contours[k]
        area = cv2.contourArea(cnt, True)
        if parent == -1:
            if abs(area) > tiss_th:
                cv2.drawContours(img_MULTI_gray_nega_mask, [cnt], -1, color=255, thickness=-1)

    filename_MULTI_gray_nega_mask = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_mask.png"        
    cv2.imwrite(filename_MULTI_gray_nega_mask, img_MULTI_gray_nega_mask)

    img_MULTI_gray_nega_mask_01 = cv2.resize(img_MULTI_gray_nega_mask, None, None, 0.1, 0.1, cv2.INTER_NEAREST)                                    
    filename_MULTI_gray_nega_mask_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_mask_01.png"
    filename_MULTI_gray_nega_mask_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_nega_mask_01.png"
    cv2.imwrite(filename_MULTI_gray_nega_mask_01, img_MULTI_gray_nega_mask_01)

    print("MASK NEGA w/HOLE")
    
    # with holes (NEGA)
    img_MULTI_gray_nega_mask_with_hole = np.zeros_like(img_MULTI_gray_nega_mask)
    for k, v in enumerate(hierarchy[0]):

        next, prev, child, parent = v
        cnt = contours[k]
        area = cv2.contourArea(cnt, True)

        if area < 0 and abs(area) > tiss_th:
            cv2.drawContours(img_MULTI_gray_nega_mask_with_hole, [cnt], -1, color=255, thickness=-1)
        elif area > 0 and abs(area) > hole_th:
            cv2.drawContours(img_MULTI_gray_nega_mask_with_hole, [cnt], -1, color=0, thickness=-1)
            

    filename_MULTI_gray_nega_mask_with_hole = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_mask_with_hole.png"        
    cv2.imwrite(filename_MULTI_gray_nega_mask_with_hole, img_MULTI_gray_nega_mask_with_hole)

    img_MULTI_gray_nega_mask_with_hole_01 = cv2.resize(img_MULTI_gray_nega_mask_with_hole, None, None, 0.1, 0.1, cv2.INTER_NEAREST)         
    filename_MULTI_gray_nega_mask_with_hole_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_mask_with_hole_01.png"
    filename_MULTI_gray_nega_mask_with_hole_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_nega_mask_with_hole_01.png"        
    cv2.imwrite(filename_MULTI_gray_nega_mask_with_hole_01, img_MULTI_gray_nega_mask_with_hole_01)

    # with holes (POS)

    print("MASK POS w/HOLE")
    
    img_MULTI_gray_nega_mask_with_hole_not = cv2.bitwise_not(img_MULTI_gray_nega_mask_with_hole)

    # BGR
    img_MULTI_gray_nega_mask_with_hole_bgr = cv2.cvtColor(img_MULTI_gray_nega_mask_with_hole, cv2.COLOR_GRAY2BGR)
    img_MULTI_gray_nega_mask_bgr = cv2.cvtColor(img_MULTI_gray_nega_mask, cv2.COLOR_GRAY2BGR)
    
    ci = 0
    cc_r = int(cm(ci)[0]*255)
    cc_g = int(cm(ci)[1]*255)
    cc_b = int(cm(ci)[2]*255)
    cc = (cc_b, cc_g, cc_r)

    cc = (0, 255, 0)
    
    # with hole
    img_MULTI_gray_nega_mask_with_hole_bgr[img_MULTI_gray_nega_mask_with_hole == 255] = royalblue # green
    img_MULTI_gray_nega_mask_with_hole_bgr[img_MULTI_gray_nega_mask_with_hole == 0] = white
    filename_MULTI_gray_nega_mask_with_hole_bgr = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_mask_with_hole_bgr.png"        
    cv2.imwrite(filename_MULTI_gray_nega_mask_with_hole_bgr, img_MULTI_gray_nega_mask_with_hole_bgr)

    img_MULTI_gray_nega_mask_with_hole_bgr_01 = cv2.resize(img_MULTI_gray_nega_mask_with_hole_bgr, None, None, 0.1, 0.1, cv2.INTER_NEAREST)
    img_MULTI_gray_nega_mask_with_hole_bgr_01 = bg_transp(img_MULTI_gray_nega_mask_with_hole_bgr_01)    
    filename_MULTI_gray_nega_mask_with_hole_bgr_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_mask_with_hole_bgr_01.png"
    filename_MULTI_gray_nega_mask_with_hole_bgr_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_nega_mask_with_hole_bgr_01.png"            
    cv2.imwrite(filename_MULTI_gray_nega_mask_with_hole_bgr_01, img_MULTI_gray_nega_mask_with_hole_bgr_01)

    # without hole
    img_MULTI_gray_nega_mask_bgr[img_MULTI_gray_nega_mask == 255] = royalblue # green
    img_MULTI_gray_nega_mask_bgr[img_MULTI_gray_nega_mask == 0] = white
    filename_MULTI_gray_nega_mask_bgr = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_mask_bgr.png"        
    cv2.imwrite(filename_MULTI_gray_nega_mask_bgr, img_MULTI_gray_nega_mask_bgr)

    img_MULTI_gray_nega_mask_bgr_01 = cv2.resize(img_MULTI_gray_nega_mask_bgr, None, None, 0.1, 0.1, cv2.INTER_NEAREST)
    img_MULTI_gray_nega_mask_bgr_01 = bg_transp(img_MULTI_gray_nega_mask_bgr_01)    
    filename_MULTI_gray_nega_mask_bgr_01 = MULTIDIR + "/" + PREFIX + "-" + "MULTI_gray_nega_mask_bgr_01.png"
    filename_MULTI_gray_nega_mask_bgr_01_rel = "dat" + "/" + "MULTI" + "/" + PREFIX + "-" + "MULTI_gray_nega_mask_bgr_01.png"        
    cv2.imwrite(filename_MULTI_gray_nega_mask_bgr_01, img_MULTI_gray_nega_mask_bgr_01)

    print("<h3>Mask generation</h3>", file=fp_html)                
    print("<img src=" + "\"" + filename_MULTI_gray_nega_mask_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)        
    print("<img src=" + "\"" + filename_MULTI_gray_nega_mask_with_hole_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    print("<img src=" + "\"" + filename_MULTI_gray_nega_mask_bgr_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)        
    print("<img src=" + "\"" + filename_MULTI_gray_nega_mask_with_hole_bgr_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    # QPTIFF Pages

    DOPAGES = True
    if DOPAGES:

        print("Takes time..., \"DOPAGES=False\" reuses the previous results.")

        filename_rgb = PAGESDIR + "/" + PREFIX + "-RGB.txt"
        fp_rgb = open(filename_rgb, "wt")
        print("Page\tName\tName(original)\tR\tG\tB", file=fp_rgb)

        for ci in range(0, len(IMAGES)):

            name, name_org, img = IMAGES[ci]  # black background
            img_PAGE_gray_nega = img.copy()

            PAGEDIR = PAGESDIR + "/" + "PAGE{:02d}".format(ci+1) + "-" + name
            qq = pathlib.Path(PAGEDIR)
            if not qq.exists():
                qq.mkdir(exist_ok=True, parents=True)
            
            # GRAY NEGA
            img_PAGE_gray_nega_01 = cv2.resize(img_PAGE_gray_nega, None, None, 0.1, 0.1, cv2.INTER_NEAREST)

            # BGR NEGA, black background (BB)
            img_PAGE_gray_nega_bgr = cv2.cvtColor(img_PAGE_gray_nega, cv2.COLOR_GRAY2BGR)            
            
            # GRAY POS, white background (WB)
            img_PAGE_gray = cv2.bitwise_not(img_PAGE_gray_nega)  
            img_PAGE_gray_01 = cv2.resize(img_PAGE_gray, None, None, 0.1, 0.1, cv2.INTER_NEAREST)
            
            # BGR POS, white background (WB)
            img_PAGE_gray_bgr = cv2.cvtColor(img_PAGE_gray, cv2.COLOR_GRAY2BGR)            
            
            cc = [int(v * 255) for v in reversed(cm(ci)[:3])] # RGB -> BGR
            [cc_b, cc_g, cc_r] = cc

            print("{:02d}".format(ci+1) + "\t" + name + "\t" + name_org + "\t" + "{:.3f}".format(cc_r) + "\t" + "{:.3f}".format(cc_g) + "\t" + "{:.3f}".format(cc_b), file=fp_rgb)

            fp_rgb.flush()
            
            img_PAGE_bgr_bb = np.zeros_like(img_PAGE_gray) # background=(0,0,0)
            img_PAGE_bgr_bb = cv2.cvtColor(img_PAGE_bgr_bb, cv2.COLOR_GRAY2BGR)

            # black background (BB)            
            img_PAGE_bgr_bb = img_PAGE_gray_nega_bgr.copy()
            index = np.all(img_PAGE_gray_nega_bgr > binary_th, axis=-1)                        
            img_PAGE_bgr_bb[index] = cc

            img_PAGE_bgr_bb_01 = cv2.resize(img_PAGE_bgr_bb, None, None, 0.1, 0.1, cv2.INTER_NEAREST)
            
            # white background (BB -> WB)
            img_PAGE_bgr_wb = img_PAGE_bgr_bb.copy()
            index = np.all(img_PAGE_bgr_bb <= binary_th, axis=-1)            
            img_PAGE_bgr_wb[index] = [255,255,255]
            
            img_PAGE_bgr_wb_01 = cv2.resize(img_PAGE_bgr_wb, None, None, 0.1, 0.1, cv2.INTER_NEAREST)

            # Histgrams
            hist_full = cv2.calcHist([img_PAGE_gray_nega], [0], None, [256], [0,256])
            hist_fg = cv2.calcHist([img_PAGE_gray_nega], [0], img_MULTI_gray_nega_mask_with_hole, [256], [0,256])
            hist_bg = cv2.calcHist([img_PAGE_gray_nega], [0], img_MULTI_gray_nega_mask_with_hole_not, [256], [0,256])

            filename_hist = PAGEDIR + "/" + PREFIX + "-" + "PAGE{:02d}".format(ci+1) + "_hist" + "-" + name + ".png"
            
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

            xs_fg = []
            ys_fg = [] 
            for x,d in enumerate(hist_fg):
                y = d[0]
                if y > 0:
                    xs_fg.append(x)
                    ys_fg.append(y)

            xs_bg = []
            ys_bg = [] 
            for x,d in enumerate(hist_bg):
                y = d[0]
                if y > 0:
                    xs_bg.append(x)
                    ys_bg.append(y)
                    
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

            ax.set_title("P{:02d} ".format(ci+1)+name, loc="left")

            ax.scatter(xs_fg, ys_fg, color="red", marker="D", s=2, label="FG", clip_on=False, zorder=3)
            ax.scatter(xs_bg, ys_bg, color="blue", marker="D", s=2, label="BG", clip_on=False, zorder=2)
            ax.scatter(xs_full, ys_full, color="green", s=6, label="FG+BG", clip_on=False, zorder=1)            

            ax.legend(loc="upper right", markerscale=2, handletextpad=-0.5, framealpha=0.0, frameon=False)

            fig.savefig(filename_hist, bbox_inches="tight", pad_inches=0.05, dpi=600, transparent=True)
            print(filename_hist)

            plt.close()


            # FILES

            # GRAY NEGA
            filename_PAGE_gray_nega_01 = PAGEDIR + "/" + PREFIX + "-" + "PAGE{:02d}".format(ci+1) + "_gray_nega" + "-" + name + "_01" + ".png"
            cv2.imwrite(filename_PAGE_gray_nega_01, img_PAGE_gray_nega_01)

            filename_PAGE_gray_nega = PAGEDIR + "/" + PREFIX + "-" + "PAGE{:02d}".format(ci+1) + "_gray_nega" + "-" + name + ".png"    
            cv2.imwrite(filename_PAGE_gray_nega, img_PAGE_gray_nega)


            # GRAY POS
            img_PAGE_gray_01 = bg_transp(img_PAGE_gray_01)
            filename_PAGE_gray_01 = PAGEDIR + "/" + PREFIX + "-" + "PAGE{:02d}".format(ci+1) + "_gray" + "-" + name + "_01" + ".png"
            cv2.imwrite(filename_PAGE_gray_01, img_PAGE_gray_01)

            filename_PAGE_gray = PAGEDIR + "/" + PREFIX + "-" + "PAGE{:02d}".format(ci+1) + "_gray" + "-" + name + ".png"    
            cv2.imwrite(filename_PAGE_gray, img_PAGE_gray)
            
            # BGR POS, black background (BB)
            filename_PAGE_bgr_bb_01 = PAGEDIR + "/" + PREFIX + "-" + "PAGE{:02d}".format(ci+1) + "_bgr_bb" + "-" + name + "_01" + ".png"            
            cv2.imwrite(filename_PAGE_bgr_bb_01, img_PAGE_bgr_bb_01)
            
            #filename_PAGE_bgr_bb = PAGEDIR + "/" + PREFIX + "-" + "PAGE{:02d}".format(ci+1) + "_bgr_bb" + "-" + name + ".png"    
            #cv2.imwrite(filename_PAGE_bgr_bb, img_PAGE_bgr_bb)
            #print(filename_PAGE_bgr_bb)


            # white background (BB -> WB)
            img_PAGE_bgr_wb_01 = bg_transp(img_PAGE_bgr_wb_01)            
            filename_PAGE_bgr_wb_01 = PAGEDIR + "/" + PREFIX + "-" + "PAGE{:02d}".format(ci+1) + "_bgr_wb" + "-" + name + "_01" + ".png"            
            cv2.imwrite(filename_PAGE_bgr_wb_01, img_PAGE_bgr_wb_01)            

            #filename_PAGE_bgr_wb = PAGEDIR + "/" + PREFIX + "-" + "PAGE{:02d}".format(ci+1) + "_bgr_wb" + "-" + name + ".png"    
            #cv2.imwrite(filename_PAGE_bgr_wb, img_PAGE_bgr_wb)
            #print(filename_PAGE_bgr_wb)

            
        fp_rgb.close()


    print("<h2>", file=fp_html)
    print("QPTIFF Pages", file=fp_html)
    print("</h2>", file=fp_html)

    for ci in range(0, len(IMAGES)):

        name, name_org, img_PAGE_gray_nega = IMAGES[ci]

        PAGEDIR = "PAGE{:02d}".format(ci+1) + "-" + name

        filename_PAGE_gray_nega_01_rel = "dat" + "/" + "PAGES" + "/" + PAGEDIR + "/" + PREFIX + "-" + "PAGE{:02d}".format(ci+1) + "_gray_nega" + "-" + name + "_01" + ".png"
        filename_PAGE_gray_01_rel = "dat" + "/" + "PAGES" + "/" + PAGEDIR + "/" + PREFIX + "-" + "PAGE{:02d}".format(ci+1) + "_gray" + "-" + name + "_01" + ".png"
        filename_PAGE_bgr_bb_01_rel = "dat" + "/" + "PAGES" + "/" + PAGEDIR + "/" + PREFIX + "-" + "PAGE{:02d}".format(ci+1) + "_bgr_bb" + "-" + name + "_01" + ".png"
        filename_PAGE_bgr_wb_01_rel = "dat" + "/" + "PAGES" + "/" + PAGEDIR + "/" + PREFIX + "-" + "PAGE{:02d}".format(ci+1) + "_bgr_wb" + "-" + name + "_01" + ".png"
        filename_hist_rel = "dat" + "/" + "PAGES" + "/" + PAGEDIR + "/" + PREFIX + "-" + "PAGE{:02d}".format(ci+1) + "_hist" + "-" + name + ".png"

        print("<h3>", file=fp_html)
        print("Page{:02d}".format(ci+1) + ":" + name_org, file=fp_html)
        print("</h3>", file=fp_html)

        print("<img src=" + "\"" + filename_PAGE_gray_nega_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
        print("<img src=" + "\"" + filename_PAGE_gray_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\""  + ">", file=fp_html)
        print("<img src=" + "\"" + filename_PAGE_bgr_bb_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\""  + ">", file=fp_html) 
        print("<img src=" + "\"" + filename_PAGE_bgr_wb_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\""  + ">", file=fp_html)
        print("<img src=" + "\"" + filename_hist_rel + "\"" + " width=\"" + ww + "\"" + ">", file=fp_html)

        


        
    
    print_html_footer(fp_html)    

    fp_html.close()

    #####

if __name__ == "__main__":

    main()

    
