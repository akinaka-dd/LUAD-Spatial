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
import matplotlib as mpl

import subprocess

import pandas as pd

from util.image import *
from util.html import *

from util.load_data import *
from util.projection import *

def main():

    if len(sys.argv) < 3:
        print("sys.argv")
        sys.exit()

    a0 = sys.argv[0]
    pp = pathlib.Path(a0)
    SRCDIR = str(pp.parent)
        
    sid = sys.argv[1]
    filename_seg = sys.argv[2]
    
    rot90 = ""
    if len(sys.argv) >= 4:
        rot90 = sys.argv[3]

    DEGREE = ""
    if rot90 == "ROTATE_90_CLOCKWISE":
        DEGREE = cv2.ROTATE_90_CLOCKWISE
    elif rot90 == "ROTATE_90_COUNTERCLOCKWISE":
        DEGREE = cv2.ROTATE_90_COUNTERCLOCKWISE

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
    MASKDIR_p_DAT_MULTI = MASKDIR_p_DAT + "/" + "MULTI"
    print("MASKDIR_p_DAT_MULTI", MASKDIR_p_DAT_MULTI)
    
    # PAGESDIR
    PAGESDIR = MASKDIR_p_DAT + "/" +  "PAGES"

    # ALIGNDIR
    ALIGNDIR = RESULTDIR + "/" + PREFIX + "-" + "ALIGN"
    ALIGNDIR_DAT = ALIGNDIR + "/" + "dat"

    # NORMDIR    
    NORMDIR = RESULTDIR + "/" + PREFIX_v + "-" + "NORM"
    NORMDIR_DAT = NORMDIR + "/" + "dat"
    
    # CLSTDIR
    CLSTDIR = RESULTDIR + "/" + PREFIX_v + "-" + "CLST"
    CLSTDIR_DAT = CLSTDIR + "/" + "dat"

    # ROIDIR 
    ROIDIR = RESULTDIR + "/" + PREFIX + "-" + "ROI"
    pp = pathlib.Path(ROIDIR)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)

    # ROIDIR_DAT
    ROIDIR_DAT = ROIDIR + "/" + "dat"
    pp = pathlib.Path(ROIDIR_DAT)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)

    #####
    # LOAD masks
    img_visium_mask_wh, img_phenocycler_mask_wh_01, img_phenocycler_mask_wh,  img_visium_mask, img_phenocycler_mask_01, img_phenocycler_mask \
        = load_masks(MASKDIR_v_DAT, PREFIX_v, MASKDIR_p_DAT_MULTI, PREFIX_p)
    
    #####
    # LOAD PhenoCycler resolution
    filename_reso = PAGESDIR + "/" + "resolution_nm.txt"
    resolution_nm, resolution_um = load_phenocycler_resolution(filename_reso)

    #####    
    # LOAD PhenoCycler pages , IMAGE[]
    filename_rgb = PAGESDIR + "/" + PREFIX_p + "-RGB.txt"
    print(filename_rgb)
    IMAGES, PAGEID, PAGEID_pos = load_phenocycler_pages(PAGESDIR, PREFIX_p, filename_rgb)

    #####
    # LOAD visium text files
    SCALEFACTORS, spot_diameter_fullres, tissue_hires_scalef, fiducial_diameter_fullres, SPOT, GRID_rc, GRID_cr, FID, FMP = \
        load_visium_text_files(DEGREE, MASKDIR_v_DAT)
    scale = tissue_hires_scalef
    
    #####
    # LOAD visium cluster text file
    filename_seurat_clusters_txt = CLSTDIR_DAT + "/" + PREFIX_v + "-SCT_seurat_clusters.txt"
    CLUSTER = load_visium_clusters(filename_seurat_clusters_txt)

    filename_seurat_clusters1_txt = CLSTDIR_DAT + "/" + PREFIX_v + "-SCT_seurat_clusters1.txt"    
    CLUSTER1 = load_visium_clusters(filename_seurat_clusters1_txt)

    #####
    # LOAD align text file
    filename_txt = ALIGNDIR_DAT + "/" + sid + "_matrix" + ".txt"
    MATRIX, rotate_total = load_matrix(filename_txt)

    #####    
    # LOAD Visium SCT
    # SCT@scale.data
    filename_sct_sd = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-scale_data.txt"
    df_sct_sd = pd.read_csv(filename_sct_sd, index_col=0)

    filename_sct_sd1 = CLSTDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-scale_data1.txt"
    df_sct_sd1 = pd.read_csv(filename_sct_sd1, index_col=0)

    #####
    # HTML 
    filename_html = ROIDIR + "/" + PREFIX + "-ROI" + ".html"
    print(filename_html)
    fp_html = open(filename_html, "wt")
    ww = "200px"    
    print_html_header(fp_html)        

    print("<h1>", file=fp_html)
    print("ROIs of Visium and PhenoCycler data" + " &mdash; " + sid, file=fp_html)
    print("</h1>", file=fp_html)

    fp_html.flush()

    #############################    
    
    # PROJECTION

    (pid, name, name_org) = IMAGES[0]
    PAGEDIR = PAGESDIR + "/" + "PAGE{:s}".format(pid) + "-" + name    
    
    filename_PAGE_gray_nega_nn = PAGEDIR + "/" + PREFIX_p + "-" + "PAGE{:s}".format(pid) + "_gray_nega" + "-" + name + "_01" + ".png"
    filename_MULTI_gray = MASKDIR_p_DAT_MULTI + "/" + PREFIX_p + "-" + "MULTI_gray.png"
    filename_MULTI_bgr = MASKDIR_p_DAT_MULTI + "/" + PREFIX_p + "-" + "MULTI_bgr_wb.png"
    filename_MULTI_bgr_bb = MASKDIR_p_DAT_MULTI + "/" + PREFIX_p + "-" + "MULTI_bgr_bb.png"
    
    filename_visium_high = MASKDIR_v_DAT + "/" + sid + "-Visium.png"
    filename_visium_ful = MASKDIR_v_DAT + "/" + PREFIX_v + "_ful.png"

    img_MULTI_gray = cv2.imread(filename_MULTI_gray, cv2.IMREAD_GRAYSCALE)
    img_MULTI_bgr = cv2.imread(filename_MULTI_bgr, cv2.IMREAD_UNCHANGED)
    img_MULTI_bgr_bb = cv2.imread(filename_MULTI_bgr_bb, cv2.IMREAD_UNCHANGED)    
    img_visium = cv2.imread(filename_visium_high, cv2.IMREAD_UNCHANGED)

    img_vv = img_visium
    img_pp = img_MULTI_gray

    img_ppv = img_pp.copy()
    img_vvp = img_vv.copy()
    
    img_pp = cv2.cvtColor(img_pp, cv2.COLOR_GRAY2BGRA)
    img_ppv = cv2.cvtColor(img_ppv, cv2.COLOR_GRAY2BGRA)    
    img_vvp = cv2.cvtColor(img_vvp, cv2.COLOR_BGR2BGRA)
    img_MULTI_bgra = cv2.cvtColor(img_MULTI_bgr, cv2.COLOR_BGR2BGRA)
    img_MULTI_bgra_bb = cv2.cvtColor(img_MULTI_bgr_bb, cv2.COLOR_BGR2BGRA)
    
    img_vv_ful = cv2.imread(filename_visium_ful, cv2.IMREAD_UNCHANGED)        
    
    hv, wv, _ = img_vv.shape    
    hv_f, wv_f, _ = img_vv_ful.shape        
    hp, wp, _ = img_pp.shape
    hv_h, wv_h = img_visium_mask.shape            

    # FRAME (Visium)
    
    x0, y0 = f_v2p(MATRIX, scale, 0, 0)
    x1, y1 = f_v2p(MATRIX, scale, wv_f, 0)
    x2, y2 = f_v2p(MATRIX, scale, wv_f, hv_f)
    x3, y3 = f_v2p(MATRIX, scale, 0, hv_f)

    xvp_min = min(x0, x1, x2, x3)
    xvp_max = max(x0, x1, x2, x3)
    yvp_min = min(y0, y1, y2, y3)
    yvp_max = max(y0, y1, y2, y3)

    xx_p = xvp_max - xvp_min
    yy_p = yvp_max - yvp_min

    xvp_min_dd = 0
    if xvp_min < 0:
        xvp_min_dd = 0 - xvp_min
    xvp_max_dd = 0
    if xvp_max > wp:
        xvp_max_dd = xvp_max - wp
        
    yvp_min_dd = 0
    if yvp_min < 0:
        yvp_min_dd = 0 - yvp_min
    yvp_max_dd = 0
    if yvp_max > hp:
        yvp_max_dd = yvp_max - hp
        
    xxx = wp + xvp_min_dd + xvp_max_dd
    yyy = hp + yvp_min_dd + yvp_max_dd

    ww_num = 200
    ww = "{:d}px".format(ww_num)
    ww2 = "{:d}px".format(int(ww_num * xxx/wp))

    ww2a = "230px"


    img_ppv_ext = np.ones((yyy, xxx))*255
    img_ppv_ext = cv2.cvtColor(img_ppv_ext.astype(np.uint8), cv2.COLOR_GRAY2BGRA)
    img_ppv_ext[yvp_min_dd:yvp_min_dd+hp, xvp_min_dd:xvp_min_dd+wp] = img_ppv

    # FRAME (PhenoCycler)
    
    (x0p, y0p) = (0 + xvp_min_dd, 0 + yvp_min_dd)
    (x1p, y1p) = (wp + xvp_min_dd, 0 + yvp_min_dd)
    (x2p, y2p) = (wp + xvp_min_dd, hp + yvp_min_dd)    
    (x3p, y3p) = (0 + xvp_min_dd, hp + yvp_min_dd)
    
    cv2.line(img_ppv_ext, (x0p, y0p), (x1p, y1p), [0, 0, 0, 255], thickness=20)
    cv2.line(img_ppv_ext, (x1p, y1p), (x2p, y2p), [0, 0, 0, 255], thickness=20)    
    cv2.line(img_ppv_ext, (x2p, y2p), (x3p, y3p), [0, 0, 0, 255], thickness=20)
    cv2.line(img_ppv_ext, (x3p, y3p), (x0p, y0p), [0, 0, 0, 255], thickness=20)

    # FRAME (Visium)
    
    x0d = x0 + xvp_min_dd
    x1d = x1 + xvp_min_dd
    x2d = x2 + xvp_min_dd
    x3d = x3 + xvp_min_dd

    y0d = y0 + yvp_min_dd
    y1d = y1 + yvp_min_dd
    y2d = y2 + yvp_min_dd
    y3d = y3 + yvp_min_dd
    
    cv2.line(img_ppv_ext, (x0d, y0d), (x1d, y1d), [0, 0, 0, 255], thickness=20)
    cv2.line(img_ppv_ext, (x1d, y1d), (x2d, y2d), [0, 0, 0, 255], thickness=20)    
    cv2.line(img_ppv_ext, (x2d, y2d), (x3d, y3d), [0, 0, 0, 255], thickness=20)
    cv2.line(img_ppv_ext, (x3d, y3d), (x0d, y0d), [0, 0, 0, 255], thickness=20)
    
    img_vv = bg_transp(img_vv)
    img_vv_ful = bg_transp(img_vv_ful)        
    
    # diameters/radius of spots
    
    xo, yo = f_v2p(MATRIX, scale, 0, 0)
    xds, yds = f_v2p(MATRIX, scale, spot_diameter_fullres, 0)
    xdf, ydf = f_v2p(MATRIX, scale, fiducial_diameter_fullres, 0)    

    dp_s = math.sqrt((xds-xo)*(xds-xo) + (yds-yo)*(yds-yo))
    rp_s = int(dp_s/2)

    dp_f = math.sqrt((xdf-xo)*(xdf-xo) + (ydf-yo)*(ydf-yo))
    rp_f = int(dp_f/2)
    
    for bc, (xv_f, yv_f, in_tissue, ci, ri) in SPOT.items():

        xx = xv_f
        yy = yv_f

        xp_f, yp_f = f_v2p(MATRIX, scale, xx, yy)

        xv_h = int(xv_f * scale)
        yv_h = int(yv_f * scale)
        v_mask = img_visium_mask_wh[yv_h, xv_h]                
        
        if v_mask != 0:

            if in_tissue == 1:
                cv2.circle(img_ppv, (xp_f, yp_f), rp_s, [0, 0, 255, 255], thickness=30)
                cv2.circle(img_ppv_ext, (xp_f+xvp_min_dd, yp_f+yvp_min_dd), rp_s, [0, 0, 255, 255], thickness=30)                        
            else:
                cv2.circle(img_ppv, (xp_f, yp_f), rp_s, [255, 0, 0, 255], thickness=10)
                cv2.circle(img_ppv_ext, (xp_f+xvp_min_dd, yp_f+yvp_min_dd), rp_s, [255, 0, 0, 255], thickness=10)
            
    # Fiducial maker positions
    
    for (xc_h, yc_h) in FMP:

        xc_f = xc_h / scale
        yc_f = yc_h / scale

        xp_f, yp_f = f_v2p(MATRIX, scale, xc_f, yc_f)

        cv2.circle(img_ppv, (xp_f, yp_f), rp_f, [0, 0, 0, 255], thickness=10)
        cv2.circle(img_ppv_ext, (xp_f+xvp_min_dd, yp_f+yvp_min_dd), rp_f, [0, 0, 0, 255], thickness=10)
    
    img_ppv = bg_transp(img_ppv)        
            
    filename_ppv = ROIDIR_DAT + "/" + PREFIX_p + "-" + "Visium-" + "map.png"
    cv2.imwrite(filename_ppv, img_ppv)

    img_ppv_01 = cv2.resize(img_ppv, None, None, 0.1, 0.1, cv2.INTER_NEAREST)
    img_ppv_01 = bg_transp(img_ppv_01)        
    filename_ppv_01 = ROIDIR_DAT + "/" + PREFIX_p + "-" + "Visium-" + "map_01.png"
    filename_ppv_01_rel = "dat" + "/" + PREFIX_p + "-" + "Visium-" + "map_01.png"                
    cv2.imwrite(filename_ppv_01, img_ppv_01)

    img_ppv_ext = bg_transp(img_ppv_ext)            
    
    filename_ppv_ext = ROIDIR_DAT + "/" + PREFIX_p + "-" + "Visium-" + "map_ext.png"
    cv2.imwrite(filename_ppv_ext, img_ppv_ext)

    img_ppv_ext_01 = cv2.resize(img_ppv_ext, None, None, 0.1, 0.1, cv2.INTER_NEAREST)
    img_ppv_ext_01 = bg_transp(img_ppv_ext_01)                
    filename_ppv_ext_01 = ROIDIR_DAT + "/" + PREFIX_p + "-" + "Visium-" + "map_ext_01.png"
    filename_ppv_ext_01_rel = "dat" + "/" + PREFIX_p + "-" + "Visium-" + "map_ext_01.png"                
    cv2.imwrite(filename_ppv_ext_01, img_ppv_ext_01)
    
    print("<h2>", file=fp_html)
    print("Visium on PhenoCycler", file=fp_html)
    print("</h2>", file=fp_html)

    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_ppv_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_ppv_ext_01_rel + "\"" + " width=\"" + ww2 + "\"" + " border=\"0\"" + ">", file=fp_html)
    
    # ROTATE

    print("ROTATE")

    def rotate_img(img, dg):
        h, w, _ = img_ppv_ext_01.shape
        xc = int(w/2)
        yc = int(h/2)
        h_max = h
        w_max = w
        img_rotated, _ = ROTATE(img, img, xc, yc, h_max, w_max, dg)

        return img_rotated

    img_ppv_ext_rotated_01 = rotate_img(img_ppv_ext_01, rotate_total)

    filename_ppv_ext_rotated_01 = ROIDIR_DAT + "/" + PREFIX_p + "-" + "Visium-" + "map_ext_rotated_01.png"
    filename_ppv_ext_rotated_01_rel = "dat" + "/" + PREFIX_p + "-" + "Visium-" + "map_ext_rotated_01.png"                
    cv2.imwrite(filename_ppv_ext_rotated_01, img_ppv_ext_rotated_01)
    
    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_ppv_ext_rotated_01_rel + "\"" + " width=\"" + ww2 + "\"" + " border=\"0\"" + ">", file=fp_html)

    fp_html.flush()

    ######    ######    ######    ######

    # P on V

    # FRAME PhenoCycler

    scale_h = 1.0

    x0, y0 = f_p2v(MATRIX, scale_h, 0, 0)
    x1, y1 = f_p2v(MATRIX, scale_h, wp, 0)
    x2, y2 = f_p2v(MATRIX, scale_h, wp, hp)
    x3, y3 = f_p2v(MATRIX, scale_h, 0, hp)
    
    cv2.line(img_vvp, (x0, y0), (x1, y1), [0, 0, 0, 255], thickness=2)
    cv2.line(img_vvp, (x1, y1), (x2, y2), [0, 0, 0, 255], thickness=2)    
    cv2.line(img_vvp, (x2, y2), (x3, y3), [0, 0, 0, 255], thickness=2)
    cv2.line(img_vvp, (x3, y3), (x0, y0), [0, 0, 0, 255], thickness=2)

    xpv_min = min(x0, x1, x2, x3)
    xpv_max = max(x0, x1, x2, x3)
    ypv_min = min(y0, y1, y2, y3)
    ypv_max = max(y0, y1, y2, y3)

    xx_p = xpv_max - xpv_min
    yy_p = ypv_max - ypv_min

    print("MIN/MAX", xpv_min, xpv_max, ypv_min, ypv_max)

    xpv_min_dd = 0
    if xpv_min < 0:
        xpv_min_dd = 0 - xpv_min
    xpv_max_dd = 0
    if xpv_max > wv_h:
        xpv_max_dd = xpv_max - wv_h
        
    ypv_min_dd = 0
    if ypv_min < 0:
        ypv_min_dd = 0 - ypv_min
    ypv_max_dd = 0
    if ypv_max > hv_h:
        ypv_max_dd = ypv_max - hv_h
        
    xxx = wv_h + xpv_min_dd + xpv_max_dd
    yyy = hv_h + ypv_min_dd + ypv_max_dd

    ww_num = 200
    ww = "{:d}px".format(ww_num)
    ww3 = "{:d}px".format(int(ww_num * xxx/wv_h))

    # Visium H&E

    print("<h2>", file=fp_html)
    print("Visium H&E", file=fp_html)
    print("</h2>", file=fp_html)

    filename_vv_HE = ROIDIR_DAT + "/" + PREFIX_v + "-HE.png"    
    cv2.imwrite(filename_vv_HE, img_vv)
    filename_vv_HE_rel = "dat" + "/" + PREFIX_v + "-HE.png"                
    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vv_HE_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    
    fp_html.flush()


    # Visium cluster ID

    print("<h2>", file=fp_html)
    print("Visium Cluster ID", file=fp_html)
    print("</h2>", file=fp_html)

    print("<h3>", file=fp_html)
    print("SCTransform with clip.range=c(-20,50) option (slice)", file=fp_html)    
    print("</h3>", file=fp_html)
    
    print('SCTransform(slice, assay="Spatial", return.only.var.genes=FALSE, verbose=TRUE, clip.range=c(-20,50))', file=fp_html)
    print("<br>", file=fp_html)
    print("CV2", file=fp_html)    
    print("<br>", file=fp_html)    

    cluster_num = len(set((CLUSTER.values())))
    print("CLUSTER_NUM", cluster_num)

    cmap = plt.get_cmap("rainbow", cluster_num)
    R, G, B, A = cmap(0 / cluster_num)
    print(R, G, B, A)
    
    rv_s = int(spot_diameter_fullres / 2)
    rv_s_h = int(scale*spot_diameter_fullres / 2)

    img_vvp_ful = img_vv_ful.copy()
    img_vvp_ful = cv2.cvtColor(img_vvp_ful, cv2.COLOR_BGR2BGRA)
    img_vvp_ful = cv2.cvtColor(img_vvp_ful, cv2.COLOR_BGR2BGRA)
    
    img_vv_cluster = img_vv.copy()
    img_vv_cluster_ful = img_vv_ful.copy()


    def getRGB(filename_seurat_clusters_rgb_txt):

        fp_rgb = open(filename_seurat_clusters_rgb_txt)
        line = fp_rgb.readline()
        RGB = {}
        for line in fp_rgb:
            line = line.rstrip()
            vals = line.split(",")
            cid = int(vals[0])-1
            rgb = vals[1]

            m = re.match("#(\S\S)(\S\S)(\S\S)", rgb)
            if m:
                r = m.group(1)
                g = m.group(2)        
                b = m.group(3)
                r = int(r, 16)
                g = int(g, 16)
                b = int(b, 16)
            else:
                print("ERROR RGB")
                sys.exit()

            RGB[cid] = (r, g, b)

        return RGB
    

    filename_seurat_clusters_rgb_txt = CLSTDIR_DAT + "/" + PREFIX_v + "-SCT_seurat_clusters_RGB.txt"        
    RGB = getRGB(filename_seurat_clusters_rgb_txt)

    count = 0
    for bc, (xv_f, yv_f, in_tissue, ci, ri) in SPOT.items():

        xx = xv_f
        yy = yv_f

        xv_h = int(xv_f * scale)
        yv_h = int(yv_f * scale)

        if not bc in CLUSTER:
            continue

        v_mask = img_visium_mask_wh[yv_h, xv_h]        
        
        cluster_id = CLUSTER[bc]
        #R, G, B, A = cmap(cluster_id / (cluster_num-1) )
        (r, g, b) = RGB[cluster_id]
        a = 1.0
        (B, G, R, A) = (b/255, g/255, r/255, a)

        if in_tissue == 1:
            if v_mask != 0:
                cv2.circle(img_vv_cluster, (xv_h, yv_h), rv_s_h, [B*255, G*255, R*255, 255], thickness=-1)
                cv2.circle(img_vv_cluster_ful, (xv_f, yv_f), rv_s, [B*255, G*255, R*255, 255], thickness=-1)

        count += 1

    print("BARCODE2", count)

    colors = []
    for k in range(0, cluster_num):
        #c = cmap(k / cluster_num)
        (r, g, b) = RGB[k]
        c = (b/255, g/255, r/255, 1.0)
        colors.append(c)
        
    plt.rcParams['font.family'] = "Arial"             
    fig, ax = plt.subplots()

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

    count = 0
    for c in colors:
        ax.plot([], [], marker="s", color=c, label=str(count), linestyle="none")
        count += 1
        
    legend = ax.legend(frameon=False, handletextpad=0, ncol=1, columnspacing=1)    
    legend_fig = legend.figure
    legend_fig.canvas.draw()
    bbox = legend.get_window_extent().transformed(legend_fig.dpi_scale_trans.inverted())

    filename_vv_cluster_legend = ROIDIR_DAT + "/" + PREFIX_v + "-cluster_legend.png"            
    fig.savefig(filename_vv_cluster_legend, bbox_inches=bbox, pad_inches=0.05, dpi=600, transparent=True)
    plt.clf()    
    plt.close()

    filename_vv_cluster = ROIDIR_DAT + "/" + PREFIX_v + "-cluster.png"    
    cv2.imwrite(filename_vv_cluster, img_vv_cluster)
    filename_vv_cluster_rel = "dat" + "/" + PREFIX_v + "-cluster.png"                
    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vv_cluster_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    filename_vv_cluster_ful = ROIDIR_DAT + "/" + PREFIX_v + "-cluster_ful.png"    
    cv2.imwrite(filename_vv_cluster_ful, img_vv_cluster_ful)
    filename_vv_cluster_ful_rel = "dat" + "/" + PREFIX_v + "-cluster_ful.png"                
    #print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vv_cluster_ful_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    filename_vv_cluster_legend = ROIDIR_DAT + "/" + PREFIX_v + "-cluster_legend.png"    
    filename_vv_cluster_legend_rel = "dat" + "/" + PREFIX_v + "-cluster_legend.png"                
    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vv_cluster_legend_rel + "\"" + " width=\"" + "30px" + "\"" + " border=\"1\"" + ">", file=fp_html)

    fp_html.flush()    

    print("<br>", file=fp_html)
    print("<br>", file=fp_html)
    print("Seurat::SpatialDimPlot", file=fp_html)        
    print("<br>", file=fp_html)        

    filename_png = CLSTDIR_DAT + "/" + PREFIX_v + "-SpatialDimPlot.png"
    filename_png_roi = ROIDIR_DAT + "/" + PREFIX_v + "-SpatialDimPlot.png"
    shutil.copy2(filename_png, filename_png_roi)
    filename_png_rel = "dat" + "/" + PREFIX_v + "-SpatialDimPlot.png"
    print("<img src=" + "\"" + filename_png_rel + "\"" + " width=\"" + ww2a + "\"" + " border=\"0\"" + ">", file=fp_html)

    print("<br>", file=fp_html)        

    filename_seurat_clusters_txt = CLSTDIR_DAT + "/" + PREFIX_v + "-SCT_seurat_clusters.txt"
    filename_seurat_clusters_txt_roi = ROIDIR_DAT + "/" + PREFIX_v + "-SCT_seurat_clusters.txt"    
    shutil.copy2(filename_seurat_clusters_txt, filename_seurat_clusters_txt_roi)    
    filename_seurat_clusters_txt_rel = "dat" + "/" + PREFIX_v + "-SCT_seurat_clusters.txt"    
    print("<a href=\"" + filename_seurat_clusters_txt_rel + "\">" + "seurat_clusters (SCT)" + "</a>", file=fp_html)

    print("<br>", file=fp_html)    
    
   
    fp_html.flush()


    # Cluster ID 2 (Seurat SpatialDimPlot)

    '''
    if False:
        command = SRCDIR + "/exec_clustering2.r" + " " + sid + " " + ROIDIR_DAT 
        print(command)
        p = subprocess.run(command, shell=True, check=True)

    '''


    print("<h3>", file=fp_html)
    print("SCTransform with default options (return.only.var.genes=FALSE) (slice1)", file=fp_html)
    print("</h3>", file=fp_html)
    
    print('SCTransform(slice1, assay = "Spatial", verbose = FALSE)', file=fp_html)
    print("<br>", file=fp_html)
    print("CV2", file=fp_html)    
    print("<br>", file=fp_html)    

    cluster1_num = len(set((CLUSTER1.values())))
    print("CLUSTER1_NUM", cluster1_num)

    cmap1 = plt.get_cmap("rainbow", cluster1_num)
    R, G, B, A = cmap1(0 / cluster1_num)
    print(R, G, B, A)
    
    rv_s = int(spot_diameter_fullres / 2)
    rv_s_h = int(scale*spot_diameter_fullres / 2)

    img_vvp_ful1 = img_vv_ful.copy()
    img_vvp_ful1 = cv2.cvtColor(img_vvp_ful, cv2.COLOR_BGR2BGRA)
    img_vvp_ful1 = cv2.cvtColor(img_vvp_ful, cv2.COLOR_BGR2BGRA)
    
    img_vv_cluster1 = img_vv.copy()
    img_vv_cluster1_ful = img_vv_ful.copy()

    filename_seurat_clusters1_rgb_txt = CLSTDIR_DAT + "/" + PREFIX_v + "-SCT_seurat_clusters1_RGB.txt"
    RGB1 = getRGB(filename_seurat_clusters1_rgb_txt)    

    count = 0
    for bc, (xv_f, yv_f, in_tissue, ci, ri) in SPOT.items():

        xx = xv_f
        yy = yv_f

        xv_h = int(xv_f * scale)
        yv_h = int(yv_f * scale)

        if not bc in CLUSTER:
            continue

        v_mask = img_visium_mask_wh[yv_h, xv_h]        
        
        cluster_id = CLUSTER1[bc]
        #R, G, B, A = cmap1(cluster_id / (cluster1_num-1) )
        (r, g, b) = RGB1[cluster_id]
        a = 1.0
        (B, G, R, A) = (b/255, g/255, r/255, a)

        if in_tissue == 1:
            if v_mask != 0:
                cv2.circle(img_vv_cluster1, (xv_h, yv_h), rv_s_h, [B*255, G*255, R*255, 255], thickness=-1)
                cv2.circle(img_vv_cluster1_ful, (xv_f, yv_f), rv_s, [B*255, G*255, R*255, 255], thickness=-1)

        count += 1

    print("BARCODE1", count)

    colors1 = []
    for k in range(0, cluster1_num):
        #c = cmap1(k / cluster1_num)
        (b, g, r) = RGB1[k]
        c = (b/255, g/255, r/255, 1.0)
        colors1.append(c)
        
    plt.rcParams['font.family'] = "Arial"             
    fig, ax = plt.subplots()

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

    count = 0
    for c in colors1:
        ax.plot([], [], marker="s", color=c, label=str(count), linestyle="none")
        count += 1
        
    legend = ax.legend(frameon=False, handletextpad=0, ncol=1, columnspacing=1)    
    legend_fig = legend.figure
    legend_fig.canvas.draw()
    bbox = legend.get_window_extent().transformed(legend_fig.dpi_scale_trans.inverted())

    filename_vv_cluster_legend1 = ROIDIR_DAT + "/" + PREFIX_v + "-cluster_legend1.png"            
    fig.savefig(filename_vv_cluster_legend1, bbox_inches=bbox, pad_inches=0.05, dpi=600, transparent=True)
    plt.clf()    
    plt.close()

    filename_vv_cluster1 = ROIDIR_DAT + "/" + PREFIX_v + "-cluster1.png"    
    cv2.imwrite(filename_vv_cluster1, img_vv_cluster1)
    filename_vv_cluster1_rel = "dat" + "/" + PREFIX_v + "-cluster1.png"                
    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vv_cluster1_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    filename_vv_cluster1_ful = ROIDIR_DAT + "/" + PREFIX_v + "-cluster1_ful.png"    
    cv2.imwrite(filename_vv_cluster1_ful, img_vv_cluster1_ful)
    filename_vv_cluster1_ful_rel = "dat" + "/" + PREFIX_v + "-cluster1_ful.png"                
    #print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vv_cluster1_ful_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    filename_vv_cluster_legend1 = ROIDIR_DAT + "/" + PREFIX_v + "-cluster_legend1.png"    
    filename_vv_cluster_legend1_rel = "dat" + "/" + PREFIX_v + "-cluster_legend1.png"                
    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vv_cluster_legend1_rel + "\"" + " width=\"" + "30px" + "\"" + " border=\"1\"" + ">", file=fp_html)
    
    fp_html.flush()

    print("<br>", file=fp_html)    
    print("<br>", file=fp_html)
    print("Seurat::SpatialDimPlot", file=fp_html)        
    print("<br>", file=fp_html)        

    filename_png1 = CLSTDIR_DAT + "/" + PREFIX_v + "-SpatialDimPlot1.png"
    filename_png1_roi = ROIDIR_DAT + "/" + PREFIX_v + "-SpatialDimPlot1.png"
    shutil.copy2(filename_png1, filename_png1_roi)
    filename_png1_rel = "dat" + "/" + PREFIX_v + "-SpatialDimPlot1.png"
    print("<img src=" + "\"" + filename_png1_rel + "\"" + " width=\"" + ww2a + "\"" + " border=\"0\"" + ">", file=fp_html)
    
    print("<br>", file=fp_html)
    
    #filename_seurat_clusters1_txt = ROIDIR_DAT + "/" + PREFIX_v + "-SCT_seurat_clusters1.txt"
    filename_seurat_clusters1_txt = CLSTDIR_DAT + "/" + PREFIX_v + "-SCT_seurat_clusters1.txt"
    filename_seurat_clusters1_txt_roi = ROIDIR_DAT + "/" + PREFIX_v + "-SCT_seurat_clusters1.txt"    
    shutil.copy2(filename_seurat_clusters1_txt, filename_seurat_clusters1_txt_roi)    
    filename_seurat_clusters1_txt_rel = "dat" + "/" + PREFIX_v + "-SCT_seurat_clusters1.txt"
    print("<a href=\"" + filename_seurat_clusters1_txt_rel + "\">" + "seurat_clusters (SCT)" + "</a>", file=fp_html)

    print("<br>", file=fp_html)    
    
    fp_html.flush()

    
    # PhenoCycler on Visium (hires and fulres)
    
    print("<h2>", file=fp_html)
    print("PhenoCycler (MULTI) on Visium", file=fp_html)
    print("</h2>", file=fp_html)

    # P

    img_MULTI_bgra_bb_roi = img_MULTI_bgra_bb.copy()
    img_MULTI_bgra_bb_roi_01 = cv2.resize(img_MULTI_bgra_bb_roi, None, None, 0.1, 0.1, cv2.INTER_NEAREST)

    hx, wx, _ = img_MULTI_bgra_bb_roi.shape
    cv2.rectangle(img_MULTI_bgra_bb_roi, (int(wx*2/5), int(hx*2/5)), (int(wx*3/5), int(hx*3/5)), (255, 255, 255, 255), thickness=10)
    filename_MULTI_map0 = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "map0.png"
    cv2.imwrite(filename_MULTI_map0, img_MULTI_bgra_bb_roi)
    filename_MULTI_map0_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "map0.png"

    hx, wx, _ = img_MULTI_bgra_bb_roi_01.shape
    cv2.rectangle(img_MULTI_bgra_bb_roi_01, (int(wx*2/5), int(hx*2/5)), (int(wx*3/5), int(hx*3/5)), (255, 255, 255, 255), thickness=3)
    filename_MULTI_map0_01 = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "map0_01.png"
    cv2.imwrite(filename_MULTI_map0_01, img_MULTI_bgra_bb_roi_01)
    filename_MULTI_map0_01_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "map0_01.png"

    #print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_MULTI_map0_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_MULTI_map0_01_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)        
    
    fp_html.flush()
    
    # P -> V
    
    PONV = True
    if PONV:

        img_vvp_dd = img_vvp.copy()
        img_vvp_dd_ful = img_vvp_ful.copy()

        h_vvp, w_vvp, _ = img_vvp_dd.shape
        h_vvp_ful, w_vvp_ful, _ = img_vvp_dd_ful.shape

        print("SHAPE1", h_vvp, w_vvp)

        for xp in range(int(wp*2/5), int(wp*3/5)):
            for yp in range(int(hp*2/5), int(hp*3/5)):

                xv, yv = f_p2v(MATRIX, scale_h, xp, yp)
                xv_f, yv_f = f_p2v(MATRIX, scale, xp, yp)            

                if xv_f < w_vvp_ful and xv_f >= 0 and yv_f < h_vvp_ful and yv_f >= 0:
                    #v = img_pp[yp][xp]
                    v = img_MULTI_bgra[yp][xp]

                    B, G, R, A = v
                    
                    v_mask = img_visium_mask_wh[yv, xv]
                    p_mask = img_phenocycler_mask_wh[yp, xp]

                    if p_mask != 0:

                        if B <= 230 or G <= 230 or R <= 230:

                            img_vvp_dd_ful[yv_f][xv_f] = v
                            img_vvp_dd[yv][xv] = v


                """
                if xv < w_vvp and xv >= 0 and yv < h_vvp and yv >= 0:
                    #v = img_pp[yp][xp]
                    v = img_MULTI_bgra[yp][xp]               
                    img_vvp_dd[yv][xv] = v
                """

        (x0, y0) = (int(wp*2/5), int(hp*2/5))
        (x1, y1) = (int(wp*3/5), int(hp*3/5))
                
        xva, yva = f_p2v(MATRIX, scale_h, x0, y0)
        xvb, yvb = f_p2v(MATRIX, scale_h, x1, y0)
        xvc, yvc = f_p2v(MATRIX, scale_h, x1, y1)
        xvd, yvd = f_p2v(MATRIX, scale_h, x0, y1)

        xva_f, yva_f = f_p2v(MATRIX, scale, x0, y0)
        xvb_f, yvb_f = f_p2v(MATRIX, scale, x1, y0)
        xvc_f, yvc_f = f_p2v(MATRIX, scale, x1, y1)
        xvd_f, yvd_f = f_p2v(MATRIX, scale, x0, y1)
        
        
        cv2.line(img_vvp_dd, (xva, yva), (xvb, yvb), [255, 255, 255, 255], thickness=3)
        cv2.line(img_vvp_dd, (xvb, yvb), (xvc, yvc), [255, 255, 255, 255], thickness=3)
        cv2.line(img_vvp_dd, (xvc, yvc), (xvd, yvd), [255, 255, 255, 255], thickness=3)
        cv2.line(img_vvp_dd, (xvd, yvd), (xva, yva), [255, 255, 255, 255], thickness=3)
        
        cv2.line(img_vvp_dd_ful, (xva_f, yva_f), (xvb_f, yvb_f), [255, 255, 255, 255], thickness=10)
        cv2.line(img_vvp_dd_ful, (xvb_f, yvb_f), (xvc_f, yvc_f), [255, 255, 255, 255], thickness=10)
        cv2.line(img_vvp_dd_ful, (xvc_f, yvc_f), (xvd_f, yvd_f), [255, 255, 255, 255], thickness=10)
        cv2.line(img_vvp_dd_ful, (xvd_f, yvd_f), (xva_f, yva_f), [255, 255, 255, 255], thickness=10)

        
        img_vvp_dd = bg_transp(img_vvp_dd)

        filename_vvp_MULTI = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "map.png"
        cv2.imwrite(filename_vvp_MULTI, img_vvp_dd)

        filename_vvp_MULTI_ful = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "map_ful.png"
        cv2.imwrite(filename_vvp_MULTI_ful, img_vvp_dd_ful)


    filename_vvp_MULTI_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "map.png"                        
    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vvp_MULTI_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)    

    filename_vvp_MULTI_ful_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "map_ful.png"                           
    #print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vvp_MULTI_ful_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

    
    fp_html.flush()


    # PAGES

    if False:

        print("<h2>", file=fp_html)
        print("PhenoCycler (PAGE) on Visium", file=fp_html)
        print("</h2>", file=fp_html)

        # DAPI blue
        # E-cadherin 04 green
        # Collagen-IV 15 brown
        # CD4 02 red
        # CD8 16 white
        # CD20 09 yellow
        # CD68 08 cyan
        # ACTA2/a-SMA 06 magenta

        # CD68 08
        # CD11c 14
        # CD56 34
        # CD20 09
        # CD3e 29
        # CD8a 16 (CD8)
        # CD4 02
        # FOXP3 31

        # CD31 04

        (pid, name, name_org) = IMAGES[4] # CD44

        PAGEDIR = PAGESDIR + "/" + "PAGE{:s}".format(pid) + "-" + name    
        filename_PAGE_gray_nega = PAGEDIR + "/" + PREFIX_p + "-" + "PAGE{:s}".format(pid) + "_gray_nega" + "-" + name + ".png"
        filename_PAGE_gray_nega_01 = PAGEDIR + "/" + PREFIX_p + "-" + "PAGE{:s}".format(pid) + "_gray_nega" + "-" + name + "_01" + ".png"

        img_PAGE_gray_nega = cv2.imread(filename_PAGE_gray_nega, cv2.IMREAD_GRAYSCALE)    
        img_PAGE_gray_nega_01 = cv2.imread(filename_PAGE_gray_nega_01, cv2.IMREAD_GRAYSCALE)

        if PONV:

            img_vvp_dd = img_vvp.copy()
            img_vvp_dd_ful = img_vvp_ful.copy()

            h_vvp, w_vvp, _ = img_vvp_dd.shape
            h_vvp_ful, w_vvp_ful, _ = img_vvp_dd_ful.shape

            print("SHAPE2", h_vvp, w_vvp)

            for xp in range(int(wp*2/5), int(wp*3/5)):
                for yp in range(int(hp*2/5), int(hp*3/5)):

                    xv, yv = f_p2v(MATRIX, scale_h, xp, yp)
                    xv_f, yv_f = f_p2v(MATRIX, scale, xp, yp)            

                    if xv_f < w_vvp_ful and xv_f >= 0 and yv_f < h_vvp_ful and yv_f >= 0:
                        #v = img_pp[yp][xp]
                        v = img_MULTI_bgra[yp][xp]

                        u = img_PAGE_gray_nega[yp][xp]

                        v_mask = img_visium_mask_wh[yv, xv]
                        p_mask = img_phenocycler_mask_wh[yp, xp]

                        if p_mask != 0:

                            if u > 30:

                                img_vvp_dd_ful[yv_f][xv_f] = [255-u, 255-u, 255, 255]
                                img_vvp_dd[yv][xv] = [255-u, 255-u, 255, 255]

                            elif u > 0:

                                img_vvp_dd_ful[yv_f][xv_f] = [255-u, 255-u, 255, 255]
                                img_vvp_dd[yv][xv] = [255-u, 255-u, 255, 255]



            img_vvp_dd = bg_transp(img_vvp_dd)

            filename_vvp_MULTI = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "page.png"
            cv2.imwrite(filename_vvp_MULTI, img_vvp_dd)

            #####

            filename_vvp_MULTI_ful = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "page_ful.png"
            cv2.imwrite(filename_vvp_MULTI_ful, img_vvp_dd_ful)

        filename_vvp_MULTI_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "page.png"                        
        print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vvp_MULTI_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

        filename_vvp_MULTI_ful_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "page_ful.png"                    
        #print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vvp_MULTI_ful_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)


        fp_html.flush()



    # P -> V, Pages, V coordinates

    if False:
        print("<h2>", file=fp_html)
        print("PhenoCycler (MULTI) on Visium, Visium coordinates", file=fp_html)
        print("</h2>", file=fp_html)

        if PONV:

            img_vvp_dd = img_vvp.copy()
            img_vvp_dd_ful = img_vvp_ful.copy()

            h_vvp, w_vvp, _ = img_vvp_dd.shape
            h_vvp_ful, w_vvp_ful, _ = img_vvp_dd_ful.shape

            print("SHAPE3", h_vvp, w_vvp)

            for xv_f in range(int(wv_f*2/5), int(wv_f*3/5)):
                for yv_f in range(int(hv_f*2/5), int(hv_f*3/5)):

                    #xv, yv = f_p2v(MATRIX, scale_h, xp, yp)
                    xv = int(xv_f * scale)
                    yv = int(yv_f * scale)

                    xv = min(xv, wv-1)
                    yv = min(yv, hv-1)                

                    #xv_f, yv_f = f_p2v(MATRIX, scale, xp, yp)

                    xp, yp = f_v2p(MATRIX, scale, xv_f, yv_f)

                    if xp < wp and xp >= 0 and yp < hp and yp >= 0:                

                        #v = img_pp[yp][xp]
                        v = img_MULTI_bgra[yp][xp]
                        B, G, R, A = v

                        u = img_PAGE_gray_nega[yp][xp]

                        v_mask = img_visium_mask_wh[yv, xv]
                        p_mask = img_phenocycler_mask_wh[yp, xp]

                        if p_mask != 0:

                            '''
                            if u > 30:
                                img_vvp_dd_ful[yv_f][xv_f] = [255-u, 255-u, 255, 255]
                                img_vvp_dd[yv][xv] = [255-u, 255-u, 255, 255]
                            elif u > 0:
                                img_vvp_dd_ful[yv_f][xv_f] = [255-u, 255-u, 255, 255]
                                img_vvp_dd[yv][xv] = [255-u, 255-u, 255, 255]
                            '''

                            if B <= 230 or G <= 230 or R <= 230:
                                img_vvp_dd_ful[yv_f][xv_f] = v
                                img_vvp_dd[yv][xv] = v

            img_vvp_dd = bg_transp(img_vvp_dd)

            filename_vvp_MULTI = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "page_vc.png"
            cv2.imwrite(filename_vvp_MULTI, img_vvp_dd)

            filename_vvp_MULTI_ful = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "page_vc_ful.png"
            cv2.imwrite(filename_vvp_MULTI_ful, img_vvp_dd_ful)


        filename_vvp_MULTI_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "page_vc.png"                
        print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vvp_MULTI_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

        filename_vvp_MULTI_ful_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "page_vc_ful.png"
        #print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vvp_MULTI_ful_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)    

        fp_html.flush()


    # Clip Region

    print("<h2>", file=fp_html)
    print("ROIs", file=fp_html)
    print("</h2>", file=fp_html)


    # CELL TYPE

    CELLTYPE = True

    if CELLTYPE:

        print("CELLTYPE")
        
        df_seg = pd.read_csv(filename_seg, sep="\t")
        print(df_seg)

        col_names = df_seg.columns

        col_names_rev = []

        for cc in col_names:

            if cc.startswith("Class"):
                cc_r = cc
            elif cc.startswith("ROI"):
                cc_r = cc                
            elif cc.startswith("Centroid X"):
                cc_r = "Centroid X"
            elif cc.startswith("Centroid Y"):
                cc_r = "Centroid Y"
            elif cc.startswith("Detection probability"):
                cc_r = cc
            elif cc.startswith("Nucleus: Area"):
                cc_r = "Nucleus Area"
            elif cc.startswith("Nucleus: Length"):
                cc_r = "Nucleus Length"
            else:

                m = re.match("(.+?):", cc)
                if m:
                    cc_r = m.group(1)
                else:
                    cc_r = cc

            print(cc, " --> ", cc_r)

            col_names_rev.append(cc_r)

        df_seg.columns = col_names_rev

        print(df_seg)


        #####

        def SetClass(row):

            binary_th = 20
            
            xc = row["Centroid X"]
            yc = row["Centroid Y"]

            row["X"] = int(xc/resolution_um)
            row["Y"] = int(yc/resolution_um)

            if row["E-Cadherin"] >= binary_th:
                row["Class"] = "Epithelial cell"
                

            elif row["Collagen-IV"] >= binary_th:
                row["Class"] = "Stromal cell"
            elif row["a-SMA-ACTA2"] >= binary_th:
                row["Class"] = "Stromal cell"
            elif row["CD31"] >= binary_th:
                row["Class"] = "Stromal cell"
                

            elif row["CD68"] >= binary_th:
                row["Class"] = "Macrophage"
            elif row["CD11c"] >= binary_th:
                row["Class"] = "Dendritic cell"
            elif row["CD20"] >= binary_th:
                row["Class"] = "B cell"

            elif row["CD3e"] >= binary_th:

                if row["CD8"] >= binary_th: # CD8a?
                    row["Class"] = "CD8+ T cell"

                    if row["FOXP3"] >= binary_th:
                        row["Class"] = "CD8+ T regulatory cell"

                    if row["CD45RO"] >= binary_th:
                        row["Class"] = "Memory CD8+ T cell"
                    else:
                        row["Class"] = "Naive CD8+ T cell"
                        
                elif row["CD4"] >= binary_th:
                    row["Class"] = "CD4+ T cell"

                    if row["FOXP3"] >= binary_th:
                        row["Class"] = "CD4+ T regulatory cell"
                    
                    elif row["CD45RO"] >= binary_th:
                        row["Class"] = "Memory CD4+ T cell"
                    else:
                        row["Class"] = "Naive CD4+ T cell"

                else:
                    row["Class"] = "Other T cell"

            #elif row["CD56"] >= binary_th:
            #    row["Class"] = "NK cell"

            else:
                row["Class"] = "Others"
            

            return row

        
        
        def SetClass1(row):

            binary_th = 10
            
            xc = row["Centroid X"]
            yc = row["Centroid Y"]

            row["X"] = int(xc/resolution_um)
            row["Y"] = int(yc/resolution_um)

            if row["CD68"] >= binary_th:
                row["Class"] = "Macrophage"
            elif row["CD11c"] >= binary_th:
                row["Class"] = "Dendritic cell"
            elif row["CD20"] >= binary_th:
                row["Class"] = "B cell"

            elif row["CD3e"] >= binary_th:

                if row["CD8"] >= binary_th: # CD8a?
                    row["Class"] = "CD8+ T cell"

                    #if row["FOXP3"] >= binary_th:
                    #    row["Class"] = "CD8+ T regulatory cell"

                    if row["CD45RO"] >= binary_th:
                        row["Class"] = "Memory CD8+ T cell"
                    else:
                        row["Class"] = "Naive CD8+ T cell"
                        
                elif row["CD4"] >= binary_th:
                    row["Class"] = "CD4+ T cell"

                    if row["FOXP3"] >= binary_th:
                        row["Class"] = "CD4+ T regulatory cell"
                    
                    elif row["CD45RO"] >= binary_th:
                        row["Class"] = "Memory CD4+ T cell"
                    else:
                        row["Class"] = "Naive CD4+ T cell"

                else:
                    row["Class"] = "Other T cell"


            elif row["E-Cadherin"] >= binary_th:
                row["Class"] = "Epithelial cell"

            elif row["Collagen-IV"] >= binary_th:
                row["Class"] = "Stromal cell"
            elif row["a-SMA-ACTA2"] >= binary_th:
                row["Class"] = "Stromal cell"
            elif row["CD31"] >= binary_th:
                row["Class"] = "Stromal cell"

            #elif row["CD56"] >= binary_th:
            #    row["Class"] = "NK cell"

            else:
                row["Class"] = "Others"
            

            return row

    #####    
    # cell type

    class_names_dic = {"Macrophage" : "Macrophage (CD68)",
                       "Dendritic cell" : "Dendritic cell (CD11c)",
                       #"NK cell" : "NK cell (CD56)", 
                       "B cell" : "B cell (CD20)",
                       #"T cell" : "T cell (CD3e)",
                       "Memory CD8+ T cell" : "Memory CD8+ T cell (CD3e/CD8/CD45RO)",
                       "Naive CD8+ T cell" : "Naive CD8+ T cell (CD3e/CD8)",
                       "Memory CD4+ T cell" : "Memory CD4+ T cell (CD3e/CD4/CD45RO)",
                       "Naive CD4+ T cell" : "Naive CD4+ T cell (CD3e/CD4)",
                       "CD8+ T regulatory cell" : "CD8+ T regulatory cell (CD3e/CD8/FOXP3)",
                       "CD4+ T regulatory cell" : "CD4+ T regulatory cell (CD3e/CD4/FOXP3)",
                       "Other T cell" : "Other T cell (CD3e)",
                       "Epithelial cell" : "Epithelial cell (E-Cadherin)",
                       "Stromal cell" : "Stromal cell (Collagen-IV,a-SMA,CD31)", 
                       "Others" : "Others",
    }

    class_names = [
        "Epithelial cell",
        "Stromal cell", 
        "Macrophage",
        "Dendritic cell",
        "B cell",
        #"CD8+ T regulatory cell",
        "Memory CD8+ T cell",
        "Naive CD8+ T cell",
        "CD4+ T regulatory cell",
        "Memory CD4+ T cell",
        "Naive CD4+ T cell",
        "Other T cell",
        #"NK cell",
        "Others",
    ]        

    colors_class_name = {
        "Epithelial cell" : "red",
        "Stromal cell" : "pink", 
        "Macrophage" : "blue", 
        "Dendritic cell" : "yellow",
         "B cell" : "orange", 
        #"CD8+ T regulatory cell" : "blue",
        "Memory CD8+ T cell" : "deepskyblue",
        "Naive CD8+ T cell" : "cyan",
        "CD4+ T regulatory cell" : "limegreen", 
        "Memory CD4+ T cell" : "lime",
        "Naive CD4+ T cell" : "greenyellow",
        "Other T cell" : "lightseagreen",
        "NK cell" : "darkviolet",
        "Others" : "lightgrey", 
    }

    

    class_names_num = len(class_names)
    cmap_ct = plt.get_cmap("rainbow", class_names_num)
        
    colors_ct = []
    for k in range(0, class_names_num):
        c = cmap_ct(k / class_names_num)
        colors_ct.append(c)

    colors_ct = list(reversed(colors_ct))


    named_colors = mpl.colors.get_named_colors_mapping()

    def getRGBA(class_name):

        color_name = colors_class_name[class_name]
        color_hex = named_colors[color_name]
        color_rgb = mpl.colors.to_rgb(color_hex)
        color_rgba = color_rgb + (1.0,)
        R, G, B, A = color_rgba # 0 - 1.0

        return color_rgba
    
    colors_ct_name = [
        "red", # Epithelial cell        
        "pink", # Stromal cell        
        "blue", # Macrophage        
        "yellow", # Dendritic cell
        "orange", # B cell
        #"blue", # CD8+ T regulatory cell
        "deepskyblue", # Memory CD8+ T cell
        "cyan", # Naive CD8+ T cell
        "limegreen", # CD4+ T regulatory cell
        "lime", # Memory CD4+ T cell
        "greenyellow", # Naive CD4+ T cell
        "lightseagreen", # Other T cell
        #"darkviolet", # NK cell
        "lightgrey", # Others
        ]

    colors_ct = []

    '''
    for x in colors_ct_name:

        color_hex = named_colors[x]
        color_rgb = mpl.colors.to_rgb(color_hex)
        color_rgba = color_rgb + (1.0,)
        
        print(x, color_rgba)

        colors_ct.append(color_rgba)
    '''
    
    for x in class_names:
        color_rgba = getRGBA(x)
        colors_ct.append(color_rgba)

    #####

    (pid, name, name_org) = IMAGES[0] 

    PAGEDIR = PAGESDIR + "/" + "PAGE{:s}".format(pid) + "-" + name    
    filename_PAGE_gray_nega = PAGEDIR + "/" + PREFIX_p + "-" + "PAGE{:s}".format(pid) + "_gray_nega" + "-" + name + ".png"
    
    img_PAGE_gray_nega = cv2.imread(filename_PAGE_gray_nega, cv2.IMREAD_GRAYSCALE)    

    def ClipRegion_PP(img_p, dg, xc, yc, x_size, y_size):

        M = cv2.getRotationMatrix2D([xc, yc], dg, 1.0)
        M[0][2] += -1*xc + x_size/2
        M[1][2] += -1*yc + y_size/2        

        img_rr = cv2.warpAffine(img_p, M, [x_size, y_size])

        return img_rr
        
        
    def ClipRegion(img_v, img_p, x0v, y0v, x1v, y1v): # obsoleted
        
        x0p, y0p = f_v2p(MATRIX, scale, x0v, y0v)
        x1p, y1p = f_v2p(MATRIX, scale, x1v, y0v)
        x2p, y2p = f_v2p(MATRIX, scale, x1v, y1v)
        x3p, y3p = f_v2p(MATRIX, scale, x0v, y1v)        

        img_v_clipped = img_v[y0v:(y1v+1), x0v:(x1v+1)]

        
        xc = int((x0p+x2p)/2)
        yc = int((y0p+y2p)/2)

        x_size = int(math.sqrt((x1p-x0p)*(x1p-x0p) + (y1p-y0p)*(y1p-y0p)))
        y_size = int(math.sqrt((x3p-x0p)*(x3p-x0p) + (y3p-y0p)*(y3p-y0p)))
        
        img_p_clipped = ClipRegion_PP(img_p, rotate_total, xc, yc, x_size, y_size)

        print(img_v.shape)        
        print(img_v_clipped.shape)

        # Visium
        filename_vvp_ROI_ful = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI_ful.png"
        cv2.imwrite(filename_vvp_ROI_ful, img_v_clipped)

        filename_vvp_ROI_ful_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI_ful.png"
        print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vvp_ROI_ful_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)    

        
        # PhenoCycler
        filename_ppv_ROI_ful = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI_pp_ful.png"
        cv2.imwrite(filename_ppv_ROI_ful, img_p_clipped)

        filename_ppv_ROI_ful_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI_pp_ful.png"
        print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_ppv_ROI_ful_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)    


    def ClipRegion_V(img_v, x0v, y0v, x1v, y1v):
        
        x0p, y0p = f_v2p(MATRIX, scale, x0v, y0v)
        x1p, y1p = f_v2p(MATRIX, scale, x1v, y0v)
        x2p, y2p = f_v2p(MATRIX, scale, x1v, y1v)
        x3p, y3p = f_v2p(MATRIX, scale, x0v, y1v)        

        img_v_clipped = img_v[y0v:(y1v+1), x0v:(x1v+1)]

        print(img_v.shape)
        
        return img_v_clipped

    
    def ClipRegion_P(img_p, x0v, y0v, x1v, y1v):
        
        x0p, y0p = f_v2p(MATRIX, scale, x0v, y0v)
        x1p, y1p = f_v2p(MATRIX, scale, x1v, y0v)
        x2p, y2p = f_v2p(MATRIX, scale, x1v, y1v)
        x3p, y3p = f_v2p(MATRIX, scale, x0v, y1v)        

        xc = int((x0p+x2p)/2)
        yc = int((y0p+y2p)/2)

        xmin = min(x0p, x1p, x2p, x3p)
        xmax = max(x0p, x1p, x2p, x3p)
        
        ymin = min(y0p, y1p, y2p, y3p)
        ymax = max(y0p, y1p, y2p, y3p)
        

        x_size = int(math.sqrt((x1p-x0p)*(x1p-x0p) + (y1p-y0p)*(y1p-y0p)))
        y_size = int(math.sqrt((x3p-x0p)*(x3p-x0p) + (y3p-y0p)*(y3p-y0p)))
        
        img_p_clipped = ClipRegion_PP(img_p, rotate_total, xc, yc, x_size, y_size)

        print(img_p_clipped.shape)

        return img_p_clipped

        # PhenoCycler
        filename_ppv_ROI_ful = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI_pp_ful.png"
        cv2.imwrite(filename_ppv_ROI_ful, img_p_clipped)

        filename_ppv_ROI_ful_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI_pp_ful.png"
        print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_ppv_ROI_ful_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)    
    

        
    def ClipRegion_Pcell(img_p, df_seg, x0v, y0v, x1v, y1v):

        print("PCELL")

        print(img_p.shape)

        x0p, y0p = f_v2p(MATRIX, scale, x0v, y0v)
        x1p, y1p = f_v2p(MATRIX, scale, x1v, y0v)
        x2p, y2p = f_v2p(MATRIX, scale, x1v, y1v)
        x3p, y3p = f_v2p(MATRIX, scale, x0v, y1v)

        print(x0p, y0p, x1p, y1p, x2p, y2p, x3p, y3p)        

        xc = int((x0p+x2p)/2)
        yc = int((y0p+y2p)/2)

        x_size = int(math.sqrt((x1p-x0p)*(x1p-x0p) + (y1p-y0p)*(y1p-y0p)))
        y_size = int(math.sqrt((x3p-x0p)*(x3p-x0p) + (y3p-y0p)*(y3p-y0p)))

        x_size_p = x_size
        y_size_p = y_size

        x_size_v = abs(x1v-x0v)
        y_size_v = abs(y1v-y0v)

        pv_ratio = x_size_p / x_size_v
        
        xmin = min(x0p, x1p, x2p, x3p)
        xmax = max(x0p, x1p, x2p, x3p)
        
        ymin = min(y0p, y1p, y2p, y3p)
        ymax = max(y0p, y1p, y2p, y3p)
        
        print(xmin, xmax)
        print(ymin, ymax)

        xmin_um = xmin*resolution_um
        xmax_um = xmax*resolution_um
        
        ymin_um = ymin*resolution_um
        ymax_um = ymax*resolution_um

        print("UM", xmin_um, xmax_um)
        print("UM", ymin_um, ymax_um)
        
        df_seg = df_seg[(df_seg["Centroid X"] >= xmin_um) & (df_seg["Centroid X"] <= xmax_um) & (df_seg["Centroid Y"] >= ymin_um) & (df_seg["Centroid Y"] <= ymax_um)]
        print(df_seg)

        # CELL TYPE
        df_seg = df_seg.apply(SetClass1, axis=1)        
        print(df_seg)
        df_seg_class = df_seg[df_seg["Class"] != "Others"]
        print(df_seg_class)
        

        def cross_product(x0, y0, x1, y1, xt, yt):

            v1x = x1 - x0
            v1y = y1 - y0

            vtx = xt - x0
            vty = yt - y0

            cp = v1x*vty - v1y*vtx

            return cp

        def region_inside(x0p, y0p, x1p, y1p, x2p, y2p, x3p, y3p, xt, yt):

            cp01 = cross_product(x0p, y0p, x1p, y1p, xt, yt)
            cp12 = cross_product(x1p, y1p, x2p, y2p, xt, yt)
            cp23 = cross_product(x2p, y2p, x3p, y3p, xt, yt)
            cp30 = cross_product(x3p, y3p, x0p, y0p, xt, yt)

            rval = False

            if cp01 > 0 and cp12 > 0 and cp23 > 0 and cp30 > 0:
                rval = True

            return rval


        print("CELLPOS LOOP")

        img_p_cell = img_p.copy()

        # cell type

        CELL_COUNT = {}
        CELL_ROWS = []
        
        col_names = df_seg.columns.values
        col = col_names[0]
        print(col, end="")
        vv = []
        for col in col_names:
            print("\t" + col, end="")
            vv.append(col)
        print()

        CELL_ROWS.append(vv)

        for i, row in df_seg.iterrows():        
            xt = int(row["X"])
            yt = int(row["Y"])

            class_name = row["Class"]

            r = region_inside(x0p, y0p, x1p, y1p, x2p, y2p, x3p, y3p, xt, yt)
            if r:

                dd = len(df_seg[(df_seg["X"] == xt) & (df_seg["Y"] == yt)])                                    

                if dd == 1:
                    
                    R, G, B, A = getRGBA(class_name)

                    color = [B*255, G*255, R*255, 255]

                    if class_name == "Others":
                        cv2.circle(img_p_cell, (xt, yt), 8, color, thickness=-1)
                    else:
                        cv2.circle(img_p_cell, (xt, yt), 8, color, thickness=-1)
                    
                    ("CIRCLE", xt, yt, color)

                    if not class_name in CELL_COUNT:
                        CELL_COUNT[class_name] = 0

                    CELL_COUNT[class_name] += 1

                    vv = []                    
                    v = row[0]
                    vv.append(v)
                    for i in range(1, len(row)):
                        v = row[i]
                        vv.append(v)

                    CELL_ROWS.append(vv)

                elif dd > 1:
                    print("CIRCLE DUPULICATED")
                    sys.exit()
                else:
                    pass
                
        img_p_cell = ClipRegion_PP(img_p_cell, rotate_total, xc, yc, x_size, y_size)

        print(img_p_cell.shape)

        #####

        return img_p_cell, CELL_COUNT, CELL_ROWS, pv_ratio

    ##### #####

    def ClipROI(x0, y0, x1, y1, df_sct_sd, CLUSTER, filename_seurat_clusters_rgb_txt, name=""):

        print("<h3>", file=fp_html)
        if name == "":
            print("ROI-{:02d}".format(roi_id), file=fp_html)
        else:
            print("ROI-{:02d} ({})".format(roi_id, name), file=fp_html)
        print("</h3>", file=fp_html)

        # Visium

        print("<h4>", file=fp_html)
        print("Visium H&E", file=fp_html)    
        print("</h4>", file=fp_html)

        cmap_sct = plt.get_cmap("rainbow", 100)
        
        def PlotGOI(df, gene_name, with_he=True): # gene of interest

            print(gene_name)

            if not gene_name in df.index:
                print("NOT EXISTS", gene_name)
                return

            if with_he: # with HE image

                img_v_clipped_spots_gn = img_v_clipped.copy()

            else: # with white back

                img_v_clipped_spots_wb = np.ones_like(img_v_clipped)
                img_v_clipped_spots_wb = img_v_clipped_spots_wb * 255
                img_v_clipped_spots_wb = cv2.cvtColor(img_v_clipped_spots_wb, cv2.COLOR_BGR2BGRA)
                img_v_clipped_spots_wb = bg_transp(img_v_clipped_spots_wb)

                img_v_clipped_spots_gn = img_v_clipped_spots_wb

            SP = []
            sig_list = []
            for bc, (xv_f, yv_f, in_tissue, ci, ri) in SPOT.items():

                if not bc in CLUSTER:
                    continue
                if xv_f >= x0 and xv_f <= x1 and yv_f >= y0 and yv_f <= y1: # inside
                    pass
                elif xv_f >= x0-rv_s and xv_f <= x1+rv_s and yv_f >= y0-rv_s and yv_f <= y1+rv_s: # on the edge
                    pass
                else: # outside
                    continue
            
                xv_h = int(xv_f * scale)
                yv_h = int(yv_f * scale)
                v_mask = img_visium_mask_wh[yv_h, xv_h]

                sig_sct = df[bc][gene_name]
        
                #R, G, B, A = cmap(sig_sct / max_sct)        
            
                if in_tissue == 1:
                    sig_list.append(sig_sct)
                    SP.append( [bc, xv_f, yv_f, sig_sct] )

            max_sct = 3            
            min_sct = -2

            for bc, xv_f, yv_f, sig_sct in SP:

                if sig_sct > max_sct:
                    sig_sct = max_sct
                elif sig_sct < min_sct:
                    sig_sct = min_sct

                R, G, B, A = cmap( (sig_sct - min_sct) /(max_sct - min_sct))        
                cv2.circle(img_v_clipped_spots_gn, (xv_f-x0, yv_f-y0), rv_s, [B*255, G*255, R*255, 255], thickness=-1)
                
            if with_he:
                
                filename_vvp_ROI_visium_spots_gn = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots-{}.png".format(roi_id, gene_name)
                cv2.imwrite(filename_vvp_ROI_visium_spots_gn, img_v_clipped_spots_gn)
                filename_vvp_ROI_visium_spots_gn_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots-{}.png".format(roi_id, gene_name)

            else:

                filename_vvp_ROI_visium_spots_gn = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots-{}_wb.png".format(roi_id, gene_name)
                cv2.imwrite(filename_vvp_ROI_visium_spots_gn, img_v_clipped_spots_gn)
                filename_vvp_ROI_visium_spots_gn_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots-{}_wb.png".format(roi_id, gene_name)


            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("<b>{:}</b>".format(gene_name), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vvp_ROI_visium_spots_gn_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

            # colorbar

            plt.rcParams['font.family'] = "Arial"

            fig, cbar = plt.subplots(figsize=(1, 8))

            vmin = min_sct
            vmax = max_sct
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

            #cmap_sct = plt.get_cmap("rainbow", 100)
            cb = mpl.colorbar.Colorbar(
                ax=cbar,
                mappable=mpl.cm.ScalarMappable(norm=norm, cmap=cmap_sct),
                orientation="vertical",
                extend="both",
                extendrect=True,
                )

            cb.set_label(gene_name, fontsize=24, loc="center")
            cb.set_ticks(ticks=[-2, -1, 0, 1, 2, 3])                        
            cb.ax.tick_params(labelsize=24) 

            filename_vvp_ROI_visium_spots_gn_cb = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots-{}_colorbar.png".format(roi_id, gene_name)

            if with_he:
                fig.savefig(filename_vvp_ROI_visium_spots_gn_cb, bbox_inches="tight", pad_inches=0.05, dpi=600, transparent=True)
                
            plt.clf()            
            plt.close()
            
            filename_vvp_ROI_visium_spots_gn_cb_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots-{}_colorbar.png".format(roi_id, gene_name)
            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            #print("<b>{:}</b>".format(gene_name), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vvp_ROI_visium_spots_gn_cb_rel + "\"" + " height=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

            fp_html.flush()                    
            
        ##### end 
        

        img_v_clipped = ClipRegion_V(img_vv_ful, x0, y0, x1, y1)

        filename_vvp_ROI_visium = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium.png".format(roi_id)
        cv2.imwrite(filename_vvp_ROI_visium, img_v_clipped)
        filename_vvp_ROI_visium_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium.png".format(roi_id)

        print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
        print("<figcaption>", file=fp_html)
        print("<b>{:}</b>".format("H&E stain"), file=fp_html)            
        print("</figcaption>", file=fp_html)                        
        print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vvp_ROI_visium_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
        print("<figcaption>", file=fp_html)
        
        print("</figcaption>", file=fp_html)                
        print("</figure>", file=fp_html)

        fp_html.flush()        

        # H&E with spot positions

        img_v_clipped_spots_he = img_v_clipped.copy()

        cmap = plt.get_cmap("rainbow", 256)        
        
        # with Spots // cluster ID

        print("SPOTS")
        
        img_v_clipped_spots = img_v_clipped.copy()

        if False:
            img_v_clipped_spots_position = np.ones_like(img_v_clipped)
            img_v_clipped_spots_position = img_v_clipped_spots_position * 255
            img_v_clipped_spots_position = cv2.cvtColor(img_v_clipped_spots_position, cv2.COLOR_BGR2BGRA)
            img_v_clipped_spots_position = bg_transp(img_v_clipped_spots_position)

            img_v_clipped_spots_position2 = np.ones_like(img_v_clipped)
            img_v_clipped_spots_position2 = img_v_clipped_spots_position2 * 255
            img_v_clipped_spots_position2 = cv2.cvtColor(img_v_clipped_spots_position2, cv2.COLOR_BGR2BGRA)
            img_v_clipped_spots_position2 = bg_transp(img_v_clipped_spots_position2)

            img_v_clipped_spots_position3 = np.ones_like(img_v_clipped)
            img_v_clipped_spots_position3 = img_v_clipped_spots_position3 * 255
            img_v_clipped_spots_position3 = cv2.cvtColor(img_v_clipped_spots_position3, cv2.COLOR_BGR2BGRA)
            img_v_clipped_spots_position3 = bg_transp(img_v_clipped_spots_position3)

        RGB = getRGB(filename_seurat_clusters_rgb_txt)

        SPOTS_POS = []

        for bc, (xv_f, yv_f, in_tissue, ci, ri) in SPOT.items():

            if not bc in CLUSTER:
                continue

            if xv_f >= x0 and xv_f <= x1 and yv_f >= y0 and yv_f <= y1: # inside
                pass
            elif xv_f >= x0-rv_s and xv_f <= x1+rv_s and yv_f >= y0-rv_s and yv_f <= y1+rv_s: # on the edge
                pass
            else: # outside
                continue
            
            xv_h = int(xv_f * scale)
            yv_h = int(yv_f * scale)
            v_mask = img_visium_mask_wh[yv_h, xv_h]        
        
            cluster_id = CLUSTER[bc]
            #R, G, B, A = cmap(cluster_id / (cluster_num-1) )
            (r, g, b) = RGB[cluster_id]
            a = 1.0
            (B, G, R, A) = (b/255, g/255, r/255, a)
            
            if in_tissue == 1:
                cv2.circle(img_v_clipped_spots_he, (xv_f-x0, yv_f-y0), rv_s, [255, 255, 255, 255], thickness=1)                                
                cv2.circle(img_v_clipped_spots, (xv_f-x0, yv_f-y0), rv_s, [B*255, G*255, R*255, 255], thickness=-1)
                if False:
                    cv2.circle(img_v_clipped_spots_position, (xv_f-x0, yv_f-y0), rv_s, [0, 0, 255, 255], thickness=1)
                    cv2.circle(img_v_clipped_spots_position2, (xv_f-x0, yv_f-y0), rv_s, [190, 190, 190, 255], thickness=1)
                    cv2.circle(img_v_clipped_spots_position3, (xv_f-x0, yv_f-y0), rv_s, [0, 0, 0, 255], thickness=1)

                SPOTS_POS.append((xv_f-x0, yv_f-y0))

        # spots on HE
        filename_vvp_ROI_visium_spots_he = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots_he.png".format(roi_id)
        cv2.imwrite(filename_vvp_ROI_visium_spots_he, img_v_clipped_spots_he)
        filename_vvp_ROI_visium_spots_he_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots_he.png".format(roi_id)
                
        # Spot position
        if False:
            filename_vvp_ROI_visium_spots_position = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots_position.png".format(roi_id)
            cv2.imwrite(filename_vvp_ROI_visium_spots_position, img_v_clipped_spots_position)
            filename_vvp_ROI_visium_spots_position_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots_position.png".format(roi_id)
            #print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vvp_ROI_visium_spots_position_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

            filename_vvp_ROI_visium_spots_position2 = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots_position2.png".format(roi_id)
            cv2.imwrite(filename_vvp_ROI_visium_spots_position2, img_v_clipped_spots_position2)
            filename_vvp_ROI_visium_spots_position2_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots_position2.png".format(roi_id)
            #print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vvp_ROI_visium_spots_position2_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

            filename_vvp_ROI_visium_spots_position3 = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots_position3.png".format(roi_id)
            cv2.imwrite(filename_vvp_ROI_visium_spots_position3, img_v_clipped_spots_position3)
            filename_vvp_ROI_visium_spots_position3_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots_position3.png".format(roi_id)
            #print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vvp_ROI_visium_spots_position3_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)

        # on HE
        print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
        print("<figcaption>", file=fp_html)
        print("<b>{:}</b>".format("Spot position"), file=fp_html)            
        print("</figcaption>", file=fp_html)                        
        print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vvp_ROI_visium_spots_he_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
        print("<figcaption>", file=fp_html)
        print("</figcaption>", file=fp_html)                
        print("</figure>", file=fp_html)

        if False:
            # red
            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("<b>{:}</b>".format("Spot position"), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vvp_ROI_visium_spots_position_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

            # gray
            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("<b>{:}</b>".format("Spot position"), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vvp_ROI_visium_spots_position2_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)   

            # black
            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("<b>{:}</b>".format("Spot position"), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vvp_ROI_visium_spots_position3_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

       
        # cluster ID
        filename_vvp_ROI_visium_spots = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots.png".format(roi_id)
        cv2.imwrite(filename_vvp_ROI_visium_spots, img_v_clipped_spots)
        filename_vvp_ROI_visium_spots_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots.png".format(roi_id)

        print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
        print("<figcaption>", file=fp_html)
        print("<b>{:}</b>".format("Cluster ID"), file=fp_html)            
        print("</figcaption>", file=fp_html)                        
        print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vvp_ROI_visium_spots_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
        print("<figcaption>", file=fp_html)
        
        print("</figcaption>", file=fp_html)                
        print("</figure>", file=fp_html)

        fp_html.flush()                        

        # wirh Spots // SCT scale.data

        gene_list = ["CDH1", "COL4A1", "CD4", "CD8A", "MS4A1", "CD68", "ACTA2", "NCAM1", "FOXP3", "PECAM1", "CD3E", "IDO1", ]

        # with HE
        for gene_name in gene_list:
            PlotGOI(df_sct_sd, gene_name)

        # colorbar without gene name
        if False:

            max_sct = 3        
            min_sct = -2

            gene_name = ""

            plt.rcParams['font.family'] = "Arial"

            fig, cbar = plt.subplots(figsize=(1, 8))

            vmin = min_sct
            vmax = max_sct
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

            #cmap_sct = plt.get_cmap("rainbow", 100)
            cb = mpl.colorbar.Colorbar(
                ax=cbar,
                mappable=mpl.cm.ScalarMappable(norm=norm, cmap=cmap_sct),
                orientation="vertical",
                extend="both",
                extendrect=True,
                )

            cb.set_label("SCTransform residual", fontsize=24, loc="center")
            cb.set_ticks(ticks=[-2, -1, 0, 1, 2, 3])                
            cb.ax.tick_params(labelsize=24) 

            filename_vvp_ROI_visium_spots_cb = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots_colorbar.png".format(roi_id)

            fig.savefig(filename_vvp_ROI_visium_spots_cb, bbox_inches="tight", pad_inches=0.05, dpi=600, transparent=True)
            plt.clf()        
            plt.close()

            filename_vvp_ROI_visium_spots_cb_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots_colorbar.png".format(roi_id)

            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            #print("<b>{:}</b>".format(gene_name), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vvp_ROI_visium_spots_cb_rel + "\"" + " height=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)

            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

            print("<br>", file=fp_html)

            fp_html.flush()                
        
        ######################
        
        # WHITE BACK

        if False:

            # H&E stain
            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("<b>{:}</b>".format("H&E stain"), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vvp_ROI_visium_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

            # on HE
            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("<b>{:}</b>".format("Spot position"), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vvp_ROI_visium_spots_he_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

            # Spot position
            # red
            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("<b>{:}</b>".format("Spot position"), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vvp_ROI_visium_spots_position_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

            # gray
            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("<b>{:}</b>".format("Spot position"), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vvp_ROI_visium_spots_position2_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)   

            # black
            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("<b>{:}</b>".format("Spot position"), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vvp_ROI_visium_spots_position3_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

            # cluster ID
            img_v_clipped_spots_wb = np.ones_like(img_v_clipped)
            img_v_clipped_spots_wb = img_v_clipped_spots_wb * 255
            img_v_clipped_spots_wb = cv2.cvtColor(img_v_clipped_spots_wb, cv2.COLOR_BGR2BGRA)
            img_v_clipped_spots_wb = bg_transp(img_v_clipped_spots_wb)

            for bc, (xv_f, yv_f, in_tissue, ci, ri) in SPOT.items():

                if not bc in CLUSTER:
                    continue

                if xv_f >= x0 and xv_f <= x1 and yv_f >= y0 and yv_f <= y1: # inside
                    pass
                elif xv_f >= x0-rv_s and xv_f <= x1+rv_s and yv_f >= y0-rv_s and yv_f <= y1+rv_s: # on the edge
                    pass
                else: # outside
                    continue

                xv_h = int(xv_f * scale)
                yv_h = int(yv_f * scale)
                v_mask = img_visium_mask_wh[yv_h, xv_h]        

                cluster_id = CLUSTER[bc]
                #R, G, B, A = cmap(cluster_id / (cluster_num-1) )
                (r, g, b) = RGB[cluster_id]
                a = 1.0
                (B, G, R, A) = (b/255, g/255, r/255, a)


                if in_tissue == 1:
                    cv2.circle(img_v_clipped_spots_wb, (xv_f-x0, yv_f-y0), rv_s, [B*255, G*255, R*255, 255], thickness=-1)

            filename_vvp_ROI_visium_spots = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots_wb.png".format(roi_id)
            cv2.imwrite(filename_vvp_ROI_visium_spots, img_v_clipped_spots_wb) # _WB
            filename_vvp_ROI_visium_spots_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_visium_spots_wb.png".format(roi_id)

            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("<b>{:}</b>".format("Cluster ID"), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vvp_ROI_visium_spots_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

            # with WB
            for gene_name in gene_list:
                PlotGOI(df_sct_sd, gene_name, with_he=False)

            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            #print("<b>{:}</b>".format(gene_name), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vvp_ROI_visium_spots_cb_rel + "\"" + " height=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)


        ########################            
        ########################

        # PhenoCycler

        p_list = [0, 4, 15, 2, 16, 9, 8, 6] + [14, 34, 29, 31, 7] + [20, 21, 10, 12] + [22, 26, 5, 18, 19, 28] + [3, 17]
        #p_list = [0, 4, 15]

        print("<h4>", file=fp_html)
        print("PhenoCycler", file=fp_html)    
        print("</h4>", file=fp_html)


        # Cell position and type

        (pid, name, name_org) = IMAGES[0]
        print(pid, name)
        PAGEDIR = PAGESDIR + "/" + "PAGE{:s}".format(pid) + "-" + name    
        filename_PAGE_gray_nega = PAGEDIR + "/" + PREFIX_p + "-" + "PAGE{:s}".format(pid) + "_gray_nega" + "-" + name + ".png"
        img_PAGE_gray_nega = cv2.imread(filename_PAGE_gray_nega, cv2.IMREAD_GRAYSCALE)


        # LOWERBOUND
        img_p_cell = img_PAGE_gray_nega.copy()
        img_p_cell = cv2.cvtColor(img_p_cell, cv2.COLOR_GRAY2BGRA)        
        lowerbound = 10
        img_p_cell[:, :, 3] = np.where(np.all(img_p_cell[:,:,0:3] < lowerbound, axis=-1), 0, 255) # if the pixel is almost black --> transparent   
        img_p_cell, CELL_COUNT, CELL_ROWS, pv_ratio = ClipRegion_Pcell(img_p_cell, df_seg, x0, y0, x1, y1)

        if False:

            filename_ppv_ROI_cp_ful = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}-celltype.png".format(roi_id)
            cv2.imwrite(filename_ppv_ROI_cp_ful, img_p_cell)
            filename_ppv_ROI_cp_ful_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}-celltype.png".format(roi_id)

            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("<b>{:}</b>".format("Cell position/type"), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_ppv_ROI_cp_ful_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            #print("p: {:.2e}<br>p_adj: {:.2e}<br>Log2FC: {:.2f}".format(pval, pval_adj, log2fc), file=fp_html)    
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

            fp_html.flush()


        # cell type with spot position
        img_p_cell_spots = img_p_cell.copy()
        for (xv, yv) in SPOTS_POS:
            xc = int(pv_ratio*xv)
            yc = int(pv_ratio*yv)
            rv = int(pv_ratio*rv_s)
            cv2.circle(img_p_cell_spots, (xc, yc), rv, [0, 0, 0, 255], thickness=2)

        filename_ppv_ROI_cp_ful_spots = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}-celltype_spots.png".format(roi_id)
        cv2.imwrite(filename_ppv_ROI_cp_ful_spots, img_p_cell_spots)
        filename_ppv_ROI_cp_ful_spots_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}-celltype_spots.png".format(roi_id)

        print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
        print("<figcaption>", file=fp_html)
        print("<b>{:}</b>".format("Cell position/type"), file=fp_html)            
        print("</figcaption>", file=fp_html)                        
        print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_ppv_ROI_cp_ful_spots_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
        print("<figcaption>", file=fp_html)
        #print("p: {:.2e}<br>p_adj: {:.2e}<br>Log2FC: {:.2f}".format(pval, pval_adj, log2fc), file=fp_html)    
        print("</figcaption>", file=fp_html)                
        print("</figure>", file=fp_html)

        fp_html.flush()
        

        # QuPath Data

        filename_ppv_ROI_qupath = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}-qupath.txt".format(roi_id)
        filename_ppv_ROI_qupath_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}-qupath.txt".format(roi_id)
        filename_ppv_ROI_qupath_wo_path = PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}-qupath.txt".format(roi_id)
        
        fp_qp = open(filename_ppv_ROI_qupath, "wt")

        for vv in CELL_ROWS:
            print(vv[0], end="", file=fp_qp)
            for v in vv[1:]:
                print("\t", end="", file=fp_qp)                                
                print(v, end="", file=fp_qp)

            print(file=fp_qp)

        fp_qp.close()

        fp_html.flush()

        
        # Pcell legend 

        plt.rcParams['font.family'] = "Arial"             
        fig, ax = plt.subplots()

        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)

        count = 0
        #for color, class_name in zip(colors_ct, class_names):
        for class_name in class_names:
            color = getRGBA(class_name)
            #ax.plot([], [], marker="o", color=color, label=class_names_dic[class_name], linestyle="none")
            ax.plot([], [], marker="o", color=color, markeredgewidth=0.5, markeredgecolor="black", label=class_names_dic[class_name], linestyle="none")

            count += 1
        
        legend = ax.legend(frameon=False, handletextpad=0, ncol=1, columnspacing=1, labelspacing=0.05)    
            
        legend_fig = legend.figure
        legend_fig.canvas.draw()
        bbox = legend.get_window_extent().transformed(legend_fig.dpi_scale_trans.inverted())
        
        filename_vv_celltype_legend = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}-celltype_legend.png".format(roi_id)            
        fig.savefig(filename_vv_celltype_legend, bbox_inches=bbox, pad_inches=0.05, dpi=600, transparent=True)
        plt.clf()        
        plt.close()

        filename_vv_celltype_legend_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}-celltype_legend.png".format(roi_id)                           

        print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
        print("<figcaption>", file=fp_html)
        #print("<b>{:}</b>".format(gene_name), file=fp_html)            
        print("</figcaption>", file=fp_html)                        
        print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vv_celltype_legend_rel + "\"" + " height=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
        print("<figcaption>", file=fp_html)
        print("</figcaption>", file=fp_html)                
        print("</figure>", file=fp_html)
        
        fp_html.flush()

        
        # Pcell count

        print("BARCHART")

        class_vals = []
        
        for class_name in class_names:
            
            val = 0

            if class_name in CELL_COUNT:
                val = CELL_COUNT[class_name]

            class_vals.append(val)

            print(class_name, val)

            
        plt.rcParams['font.family'] = "Arial"                    
        fig, ax = plt.subplots(figsize=(8, 4))

        ax.barh(list(reversed(class_names)), list(reversed(class_vals)), color=list(reversed(colors_ct)), log=True)
        fig.show()
        
        #ax.set_xlabel("Number of cells")

        ax.tick_params(axis="x", labelsize=18)
        ax.tick_params(axis="y", labelsize=18)        

        filename_vv_celltype_barchart = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}-celltype_barchart.png".format(roi_id)            
        fig.savefig(filename_vv_celltype_barchart, bbox_inches="tight", pad_inches=0.05, dpi=600, transparent=True)
        plt.clf()        
        plt.close()

        filename_vv_celltype_barchart_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}-celltype_barchart.png".format(roi_id)                           

        print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
        print("<figcaption>", file=fp_html)
        #print("<b>{:}</b>".format(gene_name), file=fp_html)            
        print("</figcaption>", file=fp_html)                        
        print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vv_celltype_barchart_rel + "\"" + " height=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
        print("<figcaption>", file=fp_html)
        print("</figcaption>", file=fp_html)                
        print("</figure>", file=fp_html)
        
        fp_html.flush()

        print("HEATMAP")

        # DAPI : 0 blue
        # CD68 : 8 cyan
        # CD4 : 2 red
        # CD8 : 16 white
        # FOXP3 : 31 gold
        # E-Cadherin : 4 green
        # a-SMA-ACTA2 : 6 magenta
        # CD3e : 29 ???
        # CD20 : 9 yellow
        # CD31 : 3 pink?
        # Ki67 : 28 orange
        # IDO1 : 17 purple
        # Collagen-IV : 15 darkred
        
        q_list = [8, 2, 16, 31, 4, 6]
        q_list = [8, 2, 16, 4, 6]
        q_list = [8, 2, 16, 31, 29]
        q_list = [8, 2, 16, 31]
        q_list = [0, 8, 2, 16, 31, 4, 6, 9, 3, 28, 17, 15]
        #q_list = [0, 8, 2, 16, 31, 4, 6, 9, 3, 28,]
        colors_q = {0 : [255, 0, 0, 255], # DAPI
                    8 : [255, 255, 0, 255], # CD68
                    2 : [0, 0, 255, 255], # CD4
                    16 : [255, 255, 255, 255], # CD8
                    4 : [0, 255, 0, 255], # E-Cadherin
                    6 : [255, 0, 255, 255], # a-SMA-ACTA2
                    #31 : [0, 220, 220, 255], # FOXP3
                    31 : [0, 215, 255, 255], # FOXP3
                    #29 : [0, 255, 0, 255], # CD3e
                    9 : [0, 255, 255, 255], # CD20
                    3 : [203, 195, 245, 255], # CD31
                    28 : [51, 135, 239, 255], # Ki67
                    15 : [12, 20, 117, 255], # Collagen-IV
                    17 : [124, 20, 117, 255], # IDO1, purple
                    }

        abody_names_dic = {
            "DAPI" : "DAPI",
            "CD68" : "CD68",
            "CD4" : "CD4",
            "CD8" : "CD8",
            "E-Cadherin" : "E-cadherin",
            "a-SMA-ACTA2" : "a-SMA/ACTA2",
            "FOXP3" : "FOXP3",
            #"CD3e" : "CD3e",
            "CD20" : "CD20",
            "CD31" : "CD31",
            "Ki67" : "Ki67",
            "Collagen-IV" : "Collagen-IV",
            "IDO1" : "IDO1",
        }

        abody_names = [
            "DAPI",
            "CD68",
            "CD4",
            "CD8",
            "E-Cadherin",
            "a-SMA-ACTA2",
            "FOXP3",
            #"CD3e",
            "CD20",
            "CD31",
            "Ki67",
            "Collagen-IV",
            "IDO1",
        ]

        colors_abody_name = {
            "DAPI" : [255, 0, 0, 255], 
            "CD68" : [255, 255, 0, 255], 
            "CD4" : [0, 0, 255, 255], 
            "CD8" : [255, 255, 255, 255],
            "E-Cadherin" : [0, 255, 0, 255],
            "a-SMA-ACTA2" : [255, 0, 255, 255],
            "FOXP3" : [0, 215, 255, 255],  
            #"CD3e" : [0, 255, 0, 255],
            "CD20" : [0, 255, 255, 255],
            "CD31" : [203, 195, 245, 255],
            "Ki67" : [51, 135, 239, 255],
            "Collagen-IV" : [12, 20, 117, 255], 
            "IDO1" : [124, 20, 117, 255], 
        }
        
        img_q_clipped = np.zeros_like(img_PAGE_gray_nega)
        img_q_clipped = ClipRegion_P(img_q_clipped, x0, y0, x1, y1)
        img_q_clipped = cv2.cvtColor(img_q_clipped, cv2.COLOR_GRAY2BGRA)

        binary_th = 25
        count = 0
        for kk in q_list:

            (pid, name, name_org) = IMAGES[kk]

            print(pid, name)

            PAGEDIR = PAGESDIR + "/" + "PAGE{:s}".format(pid) + "-" + name    
            filename_PAGE_gray_nega = PAGEDIR + "/" + PREFIX_p + "-" + "PAGE{:s}".format(pid) + "_gray_nega" + "-" + name + ".png"
            img_PAGE_gray_nega = cv2.imread(filename_PAGE_gray_nega, cv2.IMREAD_GRAYSCALE)    

            img_p_clipped = ClipRegion_P(img_PAGE_gray_nega, x0, y0, x1, y1)
            
            color = colors_q[kk]
            
            img_th = np.where(img_p_clipped <= binary_th, 0, img_p_clipped)

            if count == 0:
                img_max = img_th.copy()
                img_q_clipped[img_th > 0] = color                

            img_q_clipped[img_th > binary_th] = color
            img_max = np.where(img_th > img_max, img_th, img_max)

            count += 1

        if False:
            # heatmap without spots
            filename_ppv_ROI_hm_ful = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_heatmap.png".format(roi_id)
            cv2.imwrite(filename_ppv_ROI_hm_ful, img_q_clipped)

            filename_ppv_ROI_hm_ful_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_heatmap.png".format(roi_id)

            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("<b>{:}</b>".format("heatmap"), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_ppv_ROI_hm_ful_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

            fp_html.flush()            
        
        # heatmap with spots
        img_q_clipped_spots = img_q_clipped.copy()
        for (xv, yv) in SPOTS_POS:
            xc = int(pv_ratio*xv)
            yc = int(pv_ratio*yv)
            rv = int(pv_ratio*rv_s)
            #cv2.circle(img_q_clipped_spots, (xc, yc), rv, [190, 190, 190, 255], thickness=2)
            cv2.circle(img_q_clipped_spots, (xc, yc), rv, [255, 255, 255, 255], thickness=2)            

        filename_ppv_ROI_hm_ful_spots = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_heatmap_spots.png".format(roi_id)
        cv2.imwrite(filename_ppv_ROI_hm_ful_spots, img_q_clipped_spots)
        filename_ppv_ROI_hm_ful_spots_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_heatmap_spots.png".format(roi_id)

        print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
        print("<figcaption>", file=fp_html)
        print("<b>{:}</b>".format("heatmap"), file=fp_html)            
        print("</figcaption>", file=fp_html)                        
        print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_ppv_ROI_hm_ful_spots_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
        print("<figcaption>", file=fp_html)
        print("</figcaption>", file=fp_html)                
        print("</figure>", file=fp_html)
        
        fp_html.flush()

        # heatmap legend 

        plt.rcParams['font.family'] = "Arial"             
        fig, ax = plt.subplots()

        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)

        count = 0
        #for abody_name in abody_names:        
        for kk in q_list:

            (pid, name, name_org) = IMAGES[kk]
            
            abody_name = name
            #color = getRGBA(abody_name)

            [B, G, R, A] = [x/255 for x in colors_abody_name[abody_name]]
            color = [R, G, B, A]
            '''
            [B, G, R, A] = colors_abody_name[abody_name]            
            B = int(B/255)
            G = int(G/255)
            R = int(R/255)
            A = int(A/255)
            color = [R, G, B, A ]
            '''

            if all([c > 0.9 for c in color]):
                ax.plot([], [], marker="o", color=color, markeredgewidth=0.5, markeredgecolor="black", label=abody_names_dic[abody_name], linestyle="none")
            else:
                #ax.plot([], [], marker="o", color=color, label=abody_names_dic[abody_name], linestyle="none")
                ax.plot([], [], marker="o", color=color, markeredgewidth=0.5, markeredgecolor="black", label=abody_names_dic[abody_name], linestyle="none")
                
            count += 1
        
        legend = ax.legend(frameon=False, handletextpad=0, ncol=1, columnspacing=1, labelspacing=0.05)    
            
        legend_fig = legend.figure
        legend_fig.canvas.draw()
        bbox = legend.get_window_extent().transformed(legend_fig.dpi_scale_trans.inverted())
        
        filename_vv_heatmap_legend = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}-heatmap_legend.png".format(roi_id)            
        fig.savefig(filename_vv_heatmap_legend, bbox_inches=bbox, pad_inches=0.05, dpi=600, transparent=True)
        plt.clf()        
        plt.close()

        filename_vv_heatmap_legend_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}-heatmap_legend.png".format(roi_id)                           

        print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
        print("<figcaption>", file=fp_html)
        #print("<b>{:}</b>".format(gene_name), file=fp_html)            
        print("</figcaption>", file=fp_html)                        
        print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_vv_heatmap_legend_rel + "\"" + " height=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
        print("<figcaption>", file=fp_html)
        print("</figcaption>", file=fp_html)                
        print("</figure>", file=fp_html)
        
        fp_html.flush()

        
        # Signal

        for kk in p_list:

            (pid, name, name_org) = IMAGES[kk]

            print(pid, name)

            PAGEDIR = PAGESDIR + "/" + "PAGE{:s}".format(pid) + "-" + name    
            filename_PAGE_gray_nega = PAGEDIR + "/" + PREFIX_p + "-" + "PAGE{:s}".format(pid) + "_gray_nega" + "-" + name + ".png"
            img_PAGE_gray_nega = cv2.imread(filename_PAGE_gray_nega, cv2.IMREAD_GRAYSCALE)    

            img_p_clipped = ClipRegion_P(img_PAGE_gray_nega, x0, y0, x1, y1)
            img_p_clipped_gray = img_p_clipped.copy()
            img_p_clipped = cv2.cvtColor(img_p_clipped, cv2.COLOR_GRAY2BGRA)

            #img_p_clipped[:, :, 3] = np.where(np.all(img_p_clipped == 0, axis=-1), 0, 255)

            # LOWERBOUND
            #lowerbound = 10
            lowerbound = 0

            img_p_clipped_color = img_p_clipped.copy()

            # original gray scale
            if False:
                img_p_clipped[:, :, 3] = np.where(np.all(img_p_clipped[:,:,0:3] < lowerbound, axis=-1), 0, 255) # if the pixel is almost black --> transparent   

                filename_ppv_ROI_ful = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_PAGE{:s}-{:s}.png".format(roi_id, pid, name)
                cv2.imwrite(filename_ppv_ROI_ful, img_p_clipped)
                filename_ppv_ROI_ful_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_PAGE{:s}-{:s}.png".format(roi_id, pid, name)

                print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
                print("<figcaption>", file=fp_html)
                print("<b>{:}</b>".format(name), file=fp_html)            
                print("</figcaption>", file=fp_html)                        
                print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_ppv_ROI_ful_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
                print("<figcaption>", file=fp_html)
                print("</figcaption>", file=fp_html)                
                print("</figure>", file=fp_html)   

                fp_html.flush()

            # color

            # if the pixel is almost black --> transparent               
            lowerbound = 10
            #index = np.where(np.all(img_p_clipped_color[:,:,0:3] < lowerbound, axis=-1))
            #img_p_clipped_color[index] = [255, 255, 255, 0]

            colmap = plt.get_cmap("jet_r")
            img_p_clipped_color = (colmap(img_p_clipped_gray/255) * 255).astype(np.uint8) #[:,:,:3]
            #img_p_clipped_color = cv2.cvtColor(img_p_clipped_color, cv2.COLOR_RGBA2BGRA)

            ## spot

            img_p_clipped_color_spots = img_p_clipped_color.copy()            

            for (xv, yv) in SPOTS_POS:
                xc = int(pv_ratio*xv)
                yc = int(pv_ratio*yv)
                rv = int(pv_ratio*rv_s)
                cv2.circle(img_p_clipped_color_spots, (xc, yc), rv, [255, 255, 255, 255], thickness=2)

            # color
            if False:            
                filename_ppv_ROI_rgba_ful = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_PAGE{:s}-{:s}_rgba.png".format(roi_id, pid, name)
                cv2.imwrite(filename_ppv_ROI_rgba_ful, img_p_clipped_color)
                filename_ppv_ROI_rgba_ful_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_PAGE{:s}-{:s}_rgba.png".format(roi_id, pid, name)

                print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
                print("<figcaption>", file=fp_html)
                print("<b>{:}</b>".format(name), file=fp_html)            
                print("</figcaption>", file=fp_html)                        
                print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_ppv_ROI_rgba_ful_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
                print("<figcaption>", file=fp_html)
                print("</figcaption>", file=fp_html)                
                print("</figure>", file=fp_html)

                fp_html.flush()                        

            # color with spots
            filename_ppv_ROI_rgba_spots_ful = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_PAGE{:s}-{:s}_rgba_spots.png".format(roi_id, pid, name)
            cv2.imwrite(filename_ppv_ROI_rgba_spots_ful, img_p_clipped_color_spots)
            filename_ppv_ROI_rgba_spots_ful_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_PAGE{:s}-{:s}_rgba_spots.png".format(roi_id, pid, name)
            
            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("<b>{:}</b>".format(name), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_ppv_ROI_rgba_spots_ful_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

            fp_html.flush()            
            

        # colorbar in 256-level grayscale

        max_gray = 255
        min_gray = 10

        gene_name = ""

        plt.rcParams['font.family'] = "Arial"

        # with bb=10, RGBA
        fig, cbar = plt.subplots(figsize=(1, 8))
        vmin = 0
        vmax = max_gray
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        cmap_gray = plt.get_cmap("jet", 256)
        cmap_gray.set_under("white")
        cb = mpl.colorbar.Colorbar(
            ax=cbar,
            mappable=mpl.cm.ScalarMappable(norm=norm, cmap=cmap_gray),
            orientation="vertical",
            #extend="min",
            #extendrect=True,
            )
        cb.set_label("Pixel value", fontsize=24, loc="center")
        bb = 0
        #cb.set_ticks(ticks=[0, 0.25*(255-bb), 0.5*(255-bb), 0.75*(255-bb), 255], labels=[10, 63, 127, 191, 255])
        cb.set_ticks(ticks=[0, 0.25*(255-bb), 0.5*(255-bb), 0.75*(255-bb), 255], labels=[0, 63, 127, 191, 255])        
        cb.ax.tick_params(labelsize=24) 
        filename_ppv_ROI_cb_rgba = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_colorbar_rgba.png".format(roi_id)
        fig.savefig(filename_ppv_ROI_cb_rgba, bbox_inches="tight", pad_inches=0.05, dpi=600, transparent=True)
        plt.clf()        
        plt.close()

        filename_ppv_ROI_cb_rgba_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_colorbar_rgba.png".format(roi_id)
        print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
        print("<figcaption>", file=fp_html)
        #print("<b>{:}</b>".format(gene_name), file=fp_html)            
        print("</figcaption>", file=fp_html)                        
        print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_ppv_ROI_cb_rgba_rel + "\"" + " height=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
        print("<figcaption>", file=fp_html)
        print("</figcaption>", file=fp_html)                
        print("</figure>", file=fp_html)
        
        # without bb=10 (regular one)
        if False:
            fig, cbar = plt.subplots(figsize=(1, 8))
            vmin = 0
            vmax = max_gray
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            cmap_gray = plt.get_cmap("binary_r", 256)
            cmap_gray.set_under("white")
            cb = mpl.colorbar.Colorbar(
                ax=cbar,
                mappable=mpl.cm.ScalarMappable(norm=norm, cmap=cmap_gray),
                orientation="vertical",
                #extend="min",
                #extendrect=True,
                )
            cb.set_label("Pixel value", fontsize=24, loc="center")
            bb = 0
            #cb.set_ticks(ticks=[0, 0.25*(255-bb), 0.5*(255-bb), 0.75*(255-bb), 255], labels=[10, 63, 127, 191, 255])
            cb.set_ticks(ticks=[0, 0.25*(255-bb), 0.5*(255-bb), 0.75*(255-bb), 255], labels=[0, 63, 127, 191, 255])        
            cb.ax.tick_params(labelsize=24) 
            filename_ppv_ROI_cb = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_colorbar.png".format(roi_id)
            fig.savefig(filename_ppv_ROI_cb, bbox_inches="tight", pad_inches=0.05, dpi=600, transparent=True)
            plt.clf()        
            plt.close()

            filename_ppv_ROI_cb_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_colorbar.png".format(roi_id)
            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            #print("<b>{:}</b>".format(gene_name), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_ppv_ROI_cb_rel + "\"" + " height=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

        # with bb=10
        if False:
            fig, cbar = plt.subplots(figsize=(1, 8))
            vmin = 0
            vmax = max_gray
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            cmap_gray = plt.get_cmap("binary_r", 256)
            cmap_gray.set_under("white")
            cb = mpl.colorbar.Colorbar(
                ax=cbar,
                mappable=mpl.cm.ScalarMappable(norm=norm, cmap=cmap_gray),
                orientation="vertical",
                extend="min",
                extendrect=True,
                )
            cb.set_label("Pixel value", fontsize=24, loc="center")
            bb = 10
            #cb.set_ticks(ticks=[0, 0.25*(255-bb), 0.5*(255-bb), 0.75*(255-bb), 255], labels=[10, 63, 127, 191, 255])
            cb.set_ticks(ticks=[0, 0.25*(255-bb), 0.5*(255-bb), 0.75*(255-bb), 255], labels=[10, 63, 127, 191, 255])        
            cb.ax.tick_params(labelsize=24) 
            filename_ppv_ROI_cb_bb = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_colorbar_bb.png".format(roi_id)
            fig.savefig(filename_ppv_ROI_cb_bb, bbox_inches="tight", pad_inches=0.05, dpi=600, transparent=True)
            plt.clf()        
            plt.close()

            filename_ppv_ROI_cb_bb_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_colorbar_bb.png".format(roi_id)
            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            #print("<b>{:}</b>".format(gene_name), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_ppv_ROI_cb_bb_rel + "\"" + " height=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

        fp_html.flush()            

        # HEATMAP superimposed

        fp_html.flush()            
        
        # SPOTS

        if False:
            img_q_clipped_spots_position = np.ones_like(img_q_clipped)
            img_q_clipped_spots_position = img_q_clipped_spots_position * 255
            img_q_clipped_spots_position = cv2.cvtColor(img_q_clipped_spots_position, cv2.COLOR_BGR2BGRA)
            img_q_clipped_spots_position = bg_transp(img_q_clipped_spots_position)

            img_q_clipped_spots_position2 = np.ones_like(img_q_clipped_spots_position)
            img_q_clipped_spots_position3 = np.ones_like(img_q_clipped_spots_position)


            for (xv, yv) in SPOTS_POS:
                xc = int(pv_ratio*xv)
                yc = int(pv_ratio*yv)
                rv = int(pv_ratio*rv_s)
                cv2.circle(img_q_clipped_spots_position, (xc, yc), rv, [0, 0, 255, 255], thickness=2)
                cv2.circle(img_q_clipped_spots_position2, (xc, yc), rv, [190, 190, 190, 255], thickness=2)
                cv2.circle(img_q_clipped_spots_position3, (xc, yc), rv, [0, 0, 0, 255], thickness=2)

            filename_ppv_ROI_phenocycler_spots_position = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_phenocycler_spots_position.png".format(roi_id)
            cv2.imwrite(filename_ppv_ROI_phenocycler_spots_position, img_q_clipped_spots_position)
            filename_ppv_ROI_phenocycler_spots_position_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_phenocycler_spots_position.png".format(roi_id)

            filename_ppv_ROI_phenocycler_spots_position2 = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_phenocycler_spots_position2.png".format(roi_id)
            cv2.imwrite(filename_ppv_ROI_phenocycler_spots_position2, img_q_clipped_spots_position2)
            filename_ppv_ROI_phenocycler_spots_position2_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_phenocycler_spots_position2.png".format(roi_id)

            filename_ppv_ROI_phenocycler_spots_position3 = ROIDIR_DAT + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_phenocycler_spots_position3.png".format(roi_id)
            cv2.imwrite(filename_ppv_ROI_phenocycler_spots_position3, img_q_clipped_spots_position3)
            filename_ppv_ROI_phenocycler_spots_position3_rel = "dat" + "/" + PREFIX_v + "-" + "PhenoCycler_MULTI-" + "ROI{:02d}_phenocycler_spots_position3.png".format(roi_id)

            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("<b>{:}</b>".format("Spot position"), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_ppv_ROI_phenocycler_spots_position_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("<b>{:}</b>".format("Spot position"), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_ppv_ROI_phenocycler_spots_position2_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

            print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("<b>{:}</b>".format("Spot position"), file=fp_html)            
            print("</figcaption>", file=fp_html)                        
            print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_ppv_ROI_phenocycler_spots_position3_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
            print("<figcaption>", file=fp_html)
            print("</figcaption>", file=fp_html)                
            print("</figure>", file=fp_html)

        
        # QuPath Data
        print("<br>", file=fp_html)
        print("<a href=\"" + filename_ppv_ROI_qupath_rel + "\">" + "QuPath data " + "(" + filename_ppv_ROI_qupath_wo_path + ")" + "</a>", file=fp_html)
        
        
        fp_html.flush()
        
            
        
    def getRegion(x0, y0, x1, y1, rv_s, mm):

        dd = int(rv_s * mm)

        xx = []
        yy = []
        
        for bc, (xv_f, yv_f, in_tissue, ci, ri) in SPOT.items():

            if xv_f >= x0 and xv_f <= x1 and yv_f >= y0 and yv_f <= y1: # inside

                xx.append(xv_f)
                yy.append(yv_f)

        if len(xx) > 0:

            x0r = min(xx) - rv_s - dd
            y0r = min(yy) - rv_s - dd

            x1r = max(xx) + rv_s + dd
            y1r = max(yy) + rv_s + dd


        return (x0r, y0r, x1r, y1r)


    print('SCTransform(slice, assay="Spatial", return.only.var.genes=FALSE, verbose=TRUE, clip.range=c(-20,50))', file=fp_html)
    print("<br>", file=fp_html)
    img_vv_cluster_roi = img_vv_cluster.copy()
    filename_vv_cluster_roi = ROIDIR_DAT + "/" + PREFIX_v + "-cluster_roi.png"    
    filename_vv_cluster_roi_rel = "dat" + "/" + PREFIX_v + "-cluster_roi.png"                
    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vv_cluster_roi_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    
    print("<br>", file=fp_html)
    
    print('SCTransform(slice1, assay = "Spatial", verbose = FALSE)', file=fp_html)
    print("<br>", file=fp_html)
    img_vv_cluster1_roi = img_vv_cluster1.copy()
    filename_vv_cluster1_roi = ROIDIR_DAT + "/" + PREFIX_v + "-cluster1_roi.png"    
    filename_vv_cluster1_roi_rel = "dat" + "/" + PREFIX_v + "-cluster1_roi.png"                
    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vv_cluster1_roi_rel + "\"" + " width=\"" + ww + "\"" + " border=\"1\"" + ">", file=fp_html)
    
    
    def DrawROI(img, x0, y0, x1, y1):

        x0h = int(x0 * scale)
        y0h = int(y0 * scale)
        x1h = int(x1 * scale)
        y1h = int(y1 * scale)

        cv2.rectangle(img, (x0h,y0h), (x1h,y1h), (255, 255, 255, 255), thickness=2)
        

    mm = 0.2  # margin

    if sid == "FFPE_LUAD_3_B":

        # Cluster3-0
        roi_id = 0
        (x0, y0, x1, y1) = (2886, 7498, 3206, 7854)
        (x0r, y0r, x1r, y1r) = getRegion(x0, y0, x1, y1, rv_s, mm)
        ClipROI(x0r, y0r, x1r, y1r, df_sct_sd1, CLUSTER1, filename_seurat_clusters1_rgb_txt, "3-0")
        DrawROI(img_vv_cluster_roi, x0r, y0r, x1r, y1r)        
        DrawROI(img_vv_cluster1_roi, x0r, y0r, x1r, y1r)

        # Cluster3-1
        roi_id = 1
        (x0, y0, x1, y1) = (2310, 6778, 2642, 7122)
        (x0r, y0r, x1r, y1r) = getRegion(x0, y0, x1, y1, rv_s, mm)
        ClipROI(x0r, y0r, x1r, y1r, df_sct_sd1, CLUSTER1, filename_seurat_clusters1_rgb_txt, "3-1")
        DrawROI(img_vv_cluster_roi, x0r, y0r, x1r, y1r)                
        DrawROI(img_vv_cluster1_roi, x0r, y0r, x1r, y1r)

        # Cluster3-2
        roi_id = 2
        (x0, y0, x1, y1) = (2198, 7125, 2541, 7468)
        (x0r, y0r, x1r, y1r) = getRegion(x0, y0, x1, y1, rv_s, mm)
        ClipROI(x0r, y0r, x1r, y1r, df_sct_sd1, CLUSTER1, filename_seurat_clusters1_rgb_txt, "3-2")
        DrawROI(img_vv_cluster_roi, x0r, y0r, x1r, y1r)                
        DrawROI(img_vv_cluster1_roi, x0r, y0r, x1r, y1r)
        
        # Cluster5-0
        roi_id = 3    
        (x0, y0, x1, y1) = (3106, 6642, 3446, 7006)
        (x0r, y0r, x1r, y1r) = getRegion(x0, y0, x1, y1, rv_s, mm)
        ClipROI(x0r, y0r, x1r, y1r, df_sct_sd1, CLUSTER1, filename_seurat_clusters1_rgb_txt, "5-0")
        DrawROI(img_vv_cluster_roi, x0r, y0r, x1r, y1r)                
        DrawROI(img_vv_cluster1_roi, x0r, y0r, x1r, y1r)
        
        # Cluster5-1
        roi_id = 4
        (x0, y0, x1, y1) = (1862, 5922, 2194, 6282)
        (x0r, y0r, x1r, y1r) = getRegion(x0, y0, x1, y1, rv_s, mm)
        ClipROI(x0r, y0r, x1r, y1r, df_sct_sd1, CLUSTER1, filename_seurat_clusters1_rgb_txt, "5-1")
        DrawROI(img_vv_cluster_roi, x0r, y0r, x1r, y1r)                
        DrawROI(img_vv_cluster1_roi, x0r, y0r, x1r, y1r)
        
        # Cluster5-2    
        roi_id = 5
        (x0, y0, x1, y1) = (2538, 5786, 2874, 6162)
        (x0r, y0r, x1r, y1r) = getRegion(x0, y0, x1, y1, rv_s, mm)
        ClipROI(x0r, y0r, x1r, y1r, df_sct_sd1, CLUSTER1, filename_seurat_clusters1_rgb_txt, "5-2")
        DrawROI(img_vv_cluster_roi, x0r, y0r, x1r, y1r)                
        DrawROI(img_vv_cluster1_roi, x0r, y0r, x1r, y1r)
        
        # Cluster6-0
        roi_id = 6
        (x0, y0, x1, y1) = (4358, 5654, 4814, 6098)
        (x0r, y0r, x1r, y1r) = getRegion(x0, y0, x1, y1, rv_s, mm)
        ClipROI(x0r, y0r, x1r, y1r, df_sct_sd1, CLUSTER1, filename_seurat_clusters1_rgb_txt, "6-0")
        DrawROI(img_vv_cluster_roi, x0r, y0r, x1r, y1r)                
        DrawROI(img_vv_cluster1_roi, x0r, y0r, x1r, y1r)

        # Cluster6-1
        roi_id = 7
        (x0, y0, x1, y1) = (3792, 5652, 4228, 6088)
        (x0r, y0r, x1r, y1r) = getRegion(x0, y0, x1, y1, rv_s, mm)
        ClipROI(x0r, y0r, x1r, y1r, df_sct_sd1, CLUSTER1, filename_seurat_clusters1_rgb_txt, "6-1")
        DrawROI(img_vv_cluster_roi, x0r, y0r, x1r, y1r)                
        DrawROI(img_vv_cluster1_roi, x0r, y0r, x1r, y1r)
        
        # Cluster6-2
        roi_id = 8
        (x0, y0, x1, y1) = (3794, 5177, 4233, 5616)
        (x0r, y0r, x1r, y1r) = getRegion(x0, y0, x1, y1, rv_s, mm)
        ClipROI(x0r, y0r, x1r, y1r, df_sct_sd1, CLUSTER1, filename_seurat_clusters1_rgb_txt, "6-2")
        DrawROI(img_vv_cluster_roi, x0r, y0r, x1r, y1r)                
        DrawROI(img_vv_cluster1_roi, x0r, y0r, x1r, y1r)

        # Cluster6-3
        roi_id = 9
        (x0, y0, x1, y1) = (4264, 4669, 4710, 5115)
        (x0r, y0r, x1r, y1r) = getRegion(x0, y0, x1, y1, rv_s, mm)
        ClipROI(x0r, y0r, x1r, y1r, df_sct_sd1, CLUSTER1, filename_seurat_clusters1_rgb_txt, "6-3")
        DrawROI(img_vv_cluster_roi, x0r, y0r, x1r, y1r)                
        DrawROI(img_vv_cluster1_roi, x0r, y0r, x1r, y1r)
        

    cv2.imwrite(filename_vv_cluster_roi, img_vv_cluster_roi)        
    cv2.imwrite(filename_vv_cluster1_roi, img_vv_cluster1_roi)

    
    print("</body>", file=fp_html)
    print("</html>", file=fp_html)
        
    fp_html.close()
    

#####

    
if __name__ == "__main__":

    main()

    
