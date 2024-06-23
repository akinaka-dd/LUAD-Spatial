#!/usr/bin/env python3

import sys
import re
import pathlib
import shutil
import numpy as np
import math
import cv2

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import pandas as pd
import statistics

from util.image import *
from util.html import *

from util.load_data import *
from util.projection import *

def main():

    binary_th = 10    

    if len(sys.argv) < 3:
        print("sys.argv")
        sys.exit()

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
    MASKDIR_p_DAT_DAPI = MASKDIR_p_DAT + "/" + "DAPI"
    MASKDIR_p_DAT_MULTI = MASKDIR_p_DAT + "/" + "MULTI"
    print("MASKDIR_p_DAT_DAPI", MASKDIR_p_DAT_DAPI)
    print("MASKDIR_p_DAT_MULTI", MASKDIR_p_DAT_MULTI)

    # PAGESDIR
    PAGESDIR = MASKDIR_p_DAT + "/" +  "PAGES"

    # ALIGNDIR
    ALIGNDIR = RESULTDIR + "/" + PREFIX + "-" + "ALIGN"
    ALIGNDIR_DAT = ALIGNDIR + "/" + "dat"            

    # NORMDIR    
    NORMDIR = RESULTDIR + "/" + PREFIX_v + "-" + "NORM"
    NORMDIR_DAT = NORMDIR + "/" + "dat"

    # EVALDIR 
    EVALDIR = RESULTDIR + "/" + PREFIX + "-" + "EVAL"
    pp = pathlib.Path(EVALDIR)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)

    # EVALDIR_DAT         
    EVALDIR_DAT = EVALDIR + "/" + "dat"
    print("EVALDIR_DAT", EVALDIR_DAT)
    pp = pathlib.Path(EVALDIR_DAT)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)
        
    #####
    # LOAD masks
    img_visium_mask_wh, img_phenocycler_mask_wh_01, img_phenocycler_mask_wh,  img_visium_mask, img_phenocycler_mask_01, img_phenocycler_mask \
        = load_masks(MASKDIR_v_DAT, PREFIX_v, MASKDIR_p_DAT_MULTI, PREFIX_p)
    
    # phenocycler boundary
    filename_phenocycler_boundary_wh_01 = \
        MASKDIR_p_DAT_MULTI + "/" + PREFIX_p + "-" + "MULTI_gray_nega_bin_cont_boundary_with_hole_01" + ".png"
    img_phenocycler_boundary_wh_01 = cv2.imread(filename_phenocycler_boundary_wh_01, cv2.IMREAD_GRAYSCALE)
    # x10
    img_phenocycler_boundary_wh = cv2.resize(img_phenocycler_boundary_wh_01, None, None, 10.0, 10.0, cv2.INTER_CUBIC)

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
    # LOAD align text file
    filename_txt = ALIGNDIR_DAT + "/" + sid + "_matrix" + ".txt"
    MATRIX, rotate_total = load_matrix(filename_txt)

    #####
    # HTML 
    filename_html = EVALDIR + "/" + PREFIX + "-EVAL" + ".html"
    print(filename_html)
    fp_html = open(filename_html, "wt")
    ww = "200px"
    ww2 = "160px"    
    print_html_header(fp_html)    

    print("<h1>", file=fp_html)
    print("Evaluation of Visium and PhenoCycler data" + " &mdash; " + sid, file=fp_html)
    print("</h1>", file=fp_html)

    fp_html.flush()
    
    #############################

    # diameters/radius of spots
    
    xo, yo = f_v2p(MATRIX, scale, 0, 0)
    xds, yds = f_v2p(MATRIX, scale, spot_diameter_fullres, 0)
    xdf, ydf = f_v2p(MATRIX, scale, fiducial_diameter_fullres, 0)    

    dp_s = math.sqrt((xds-xo)*(xds-xo) + (yds-yo)*(yds-yo))
    rp_s = int(dp_s/2)

    dp_f = math.sqrt((xdf-xo)*(xdf-xo) + (ydf-yo)*(ydf-yo))
    rp_f = int(dp_f/2)

    #####

    def CountCell(filename_seg, sid):

        pos = []

        fp = open(filename_seg, "rt")

        line = fp.readline()
        for line in fp:

            line = line.rstrip()
            vals = line.split("\t")

            if sid == "FFPE_LUAD_2_C":
                xc = float(vals[1])
                yc = float(vals[2])
                area = float(vals[4])

            else:
                xc = float(vals[2])
                yc = float(vals[3])
                area = float(vals[5])

            xc = xc/resolution_um
            yc = yc/resolution_um
            
            pos.append((xc, yc, area))

        return pos

    #####

    def PlotCor(df, df_cc, df_dd, name_v, name_p):

        # PhenoCycler

        print("PHENOCYCLER")

        count, pid, filename_PAGE_gray_nega = PAGEID[name_p]
        print(filename_PAGE_gray_nega)
        img_pp = cv2.imread(filename_PAGE_gray_nega, cv2.IMREAD_GRAYSCALE) # nega
        img_pp_01 = cv2.resize(img_pp, None, None, 0.1, 0.1, cv2.INTER_NEAREST) # nega

        _, _, filename_PAGE_gray = PAGEID_pos[name_p]
        print(filename_PAGE_gray)    
        img_pp_pos = cv2.imread(filename_PAGE_gray, cv2.IMREAD_GRAYSCALE) # pos
        img_pp_pos_01 = cv2.resize(img_pp_pos, None, None, 0.1, 0.1, cv2.INTER_NEAREST) # pos

        # clipping
        img_pp_pos_clipped_01 = cv2.bitwise_not(cv2.bitwise_and(cv2.bitwise_not(img_pp_pos_01), img_phenocycler_mask_wh_01))
        img_pp_pos_clipped = cv2.bitwise_not(cv2.bitwise_and(cv2.bitwise_not(img_pp_pos), img_phenocycler_mask_wh))

        img_pp_clipped_01 = cv2.bitwise_not(cv2.bitwise_and(cv2.bitwise_not(img_pp_pos_01), img_phenocycler_mask_wh_01))
        img_pp_clipped = cv2.bitwise_not(cv2.bitwise_and(cv2.bitwise_not(img_pp_pos), img_phenocycler_mask_wh))

        img_pp_bgra = cv2.cvtColor(img_pp, cv2.COLOR_GRAY2BGRA)
        img_pp_bgra_01 = cv2.resize(img_pp_bgra, None, None, 0.1, 0.1, cv2.INTER_NEAREST)        

        img_pp_pos_bgra = cv2.cvtColor(img_pp_pos, cv2.COLOR_GRAY2BGRA)
        img_pp_pos_bgra_01 = cv2.resize(img_pp_pos_bgra, None, None, 0.1, 0.1, cv2.INTER_NEAREST)

        # base = tissue
        img_pp_pos_bgra_base_01 = np.ones_like(img_pp_pos_bgra_01)
        img_pp_pos_bgra_base_01 = img_pp_pos_bgra_base_01 * 255
        img_pp_pos_bgra_base_01[np.where(img_phenocycler_mask_wh_01 == 255)] = [240, 240, 240, 255]

        img_pp_pos_bgra_base_plot_bin_01 = img_pp_pos_bgra_base_01.copy()
        img_pp_pos_bgra_base_plot_gray_01 = img_pp_pos_bgra_base_01.copy()

        # BINARY
        img_pp_pos_bgra_base_plot_bin_01[np.where((img_pp_pos_clipped_01 >= 0) & (img_pp_pos_clipped_01 < 255-binary_th) )] = [255, 0, 0, 255]    

        # GRAYSCALE
        img_pp_pos_bgra_01 = cv2.cvtColor(img_pp_pos_clipped_01, cv2.COLOR_GRAY2BGRA)
        img_pp_pos_bgra_01[:, :, 0] = 255  # blue

        img_pp_pos_bin_01 = np.where(img_pp_pos_clipped_01 < 255, 0, img_pp_pos_clipped_01)
        img_pp_pos_bin_bgra_01 = cv2.cvtColor(img_pp_pos_bin_01, cv2.COLOR_GRAY2BGRA)

        img_pp_pos_bgra_base_plot_gray_01 = cv2.bitwise_and(img_pp_pos_bgra_base_plot_gray_01, img_pp_pos_bin_bgra_01)
        img_pp_pos_bgra_base_plot_gray_01 = img_pp_pos_bgra_base_plot_gray_01 + img_pp_pos_bgra_01

        # boundary
        img_pp_pos_bgra_base_plot_bin_01[np.where(img_phenocycler_boundary_wh_01 == 255)] = [0, 0, 0, 255]
        img_pp_pos_bgra_base_plot_gray_01[np.where(img_phenocycler_boundary_wh_01 == 255)] = [0, 0, 0, 255]

        vals_v = []    
        vals_p = []
        vals_c = []
        vals_ec_nc = []
        vals_ec_nf = []        
        vals_nc = []
        vals_nf = []
        vals_d = []

        vals_v2 = []    
        vals_p2 = []
        
        
        cm = plt.get_cmap("rainbow")
        vmin_nc = 2000
        vmax_nc = 20000
        norm_nc = mpl.colors.Normalize(vmin=vmin_nc, vmax=vmax_nc)

        vmin_nf = 2000
        vmax_nf = 8000
        norm_nf = mpl.colors.Normalize(vmin=vmin_nf, vmax=vmax_nf)
        
        m_count_th = math.pi * rp_s * rp_s * 0.0
        print("M_COUNT_TH", m_count_th)
        
        count_green = 0
        count = 0
        for bc in df.columns.values:        

            (xv_f, yv_f, in_tissue, ci, ri) = SPOT[bc]

            xp_f, yp_f = f_v2p(MATRIX, scale, xv_f, yv_f)

            val_v = df.loc[name_v, bc]
            val_v_cc = df_cc.loc[name_v, bc]
            val_v_nc = df_nc.loc[bc, "nCount_Spatial"]
            val_v_nf = df_nf.loc[bc, "nFeature_Spatial"]
            val_v_dd = df_dd.loc[name_v, bc]
            
            v_max = 0
            v_ave = 0
            v_num = 0
            vs = []
            m_count = 0
            v_count = 0
            for dx in range(-1*rp_s, rp_s+1):
                for dy in range(-1*rp_s, rp_s+1):

                    x = xp_f + dx
                    y = yp_f + dy
                    v = img_pp[y, x]
                    m = img_phenocycler_mask_wh[y, x]

                    d = math.sqrt( (x-xp_f)*(x-xp_f) + (y-yp_f)*(y-yp_f) )

                    if v > v_max:
                        v_max = v
                    if d > rp_s:
                        continue
                    if m > 0:
                        m_count += 1
                    if v > 0:
                        vs.append(v)
                        v_count += 1

            val_p = statistics.mean(vs)
            #val_p = statistics.median(vs)
            #val_p = statistics.mean(vs_sorted[:kk])
            #val_p = max(vs)            
            
            print(count, bc, xv_f, yv_f, in_tissue, ci, ri, val_v, xp_f, yp_f, val_p, m_count)


            flag = True

            if val_v_nc < 3000:

                vals_c.append("blue")
                
            elif m_count <= m_count_th or in_tissue == 0:

                if val_v >= 2.0:
                    cv2.circle(img_pp_pos_bgra_base_plot_bin_01, (int(xp_f/10), int(yp_f/10)), int(rp_s/10), [0, 128, 0, 255], thickness=2, lineType=cv2.LINE_AA)        
                    cv2.circle(img_pp_pos_bgra_base_plot_gray_01, (int(xp_f/10), int(yp_f/10)), int(rp_s/10), [0, 128, 0, 255], thickness=1, lineType=cv2.LINE_AA)

                count_green += 1

                vals_c.append("darkgreen")

                flag = False

            elif val_v >= 2.0:
                cv2.circle(img_pp_pos_bgra_base_plot_bin_01, (int(xp_f/10), int(yp_f/10)), int(rp_s/10), [0, 0, 255, 255], thickness=2, lineType=cv2.LINE_AA)        
                cv2.circle(img_pp_pos_bgra_base_plot_gray_01, (int(xp_f/10), int(yp_f/10)), int(rp_s/10), [0, 0, 255, 255], thickness=1, lineType=cv2.LINE_AA)

                vals_c.append("red")
                
                vals_v2.append(val_v)
                vals_p2.append(val_p)

            else:

                vals_c.append("red")

                vals_v2.append(val_v)
                vals_p2.append(val_p)
                

            vals_v.append(val_v)
            vals_p.append(val_p)
            vals_nc.append(val_v_nc)
            vals_nf.append(val_v_nf)
            vals_d.append(val_v_dd)

            # NC
            v_ec_nc = val_v_nc
            if v_ec_nc > vmax_nc:
                v_ec_nc = vmax_nc
            elif v_ec_nc < vmin_nc:
                v_ec_nc = vmin_nc

            c_ec_nc = cm(norm_nc(v_ec_nc))

            if flag:
                vals_ec_nc.append(c_ec_nc)
            else:
                vals_ec_nc.append((0.5, 0.5, 0.5, 1.0))                

            # NF
            v_ec_nf = val_v_nf
            if v_ec_nf > vmax_nf:
                v_ec_nf = vmax_nf
            elif v_ec_nf < vmin_nf:
                v_ec_nf = vmin_nf

            c_ec_nf = cm(norm_nf(v_ec_nf))

            if flag:
                vals_ec_nf.append(c_ec_nf)
            else:
                vals_ec_nf.append((0.5, 0.5, 0.5, 1.0))                

            count += 1

        #####

        # High Visium spots on PhenoCycler

        print("<h2>", file=fp_html)
        print("Visium:" + name_v + " / " + "PhenoCycler:" + name_p, file=fp_html)
        print("</h2>", file=fp_html)


        # BINARY
        filename_plot_png = EVALDIR_DAT + "/" + sid + "-" + name_v + "-" + name_p + "_plot_bin_01.png"
        filename_plot_png_rel = "dat" + "/" + sid + "-" + name_v + "-" + name_p + "_plot_bin_01.png"
        img_pp_pos_bgra_base_plot_bin_01 = bg_transp(img_pp_pos_bgra_base_plot_bin_01)
        cv2.imwrite(filename_plot_png, img_pp_pos_bgra_base_plot_bin_01)
        print("<img src=" + "\"" + filename_plot_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

        # GRAYSCALE
        filename_plot_png = EVALDIR_DAT + "/" + sid + "-" + name_v + "-" + name_p + "_plot_gray_01.png"
        filename_plot_png_rel = "dat" + "/" + sid + "-" + name_v + "-" + name_p + "_plot_gray_01.png"
        img_pp_pos_bgra_base_plot_gray_01 = bg_transp(img_pp_pos_bgra_base_plot_gray_01)
        cv2.imwrite(filename_plot_png, img_pp_pos_bgra_base_plot_gray_01)
        print("<img src=" + "\"" + filename_plot_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)


        # Visium x PhenoCycler correlation (val_v_nc >= 3000 and m_count > m_count_th and in_tissue == 1 spots) 
        
        plt.rcParams['font.family'] = "Arial"        
        fig, ax = plt.subplots(figsize=(9,9))

        ax.set_box_aspect (1)

        ddv = 0.05*(max(vals_v2) - min(vals_v2))
        ax.set_xlim(min(vals_v2)-ddv, max(vals_v2)+ddv)

        ddp = 0.05*(max(vals_p2) - min(vals_p2))
        ax.set_ylim(0-ddp/2, max(vals_p2)+ddp)        

        mappable0 = ax.scatter(vals_v2, vals_p2, color=(1.0, 1.0, 1.0, 0.0), s=20, lw=1, ec="red", cmap=cm, vmin=vmin_nc, vmax=vmax_nc)

        ax.set_xlabel("Visium SCTransform residual (" + name_v + ")", fontsize=20)
        ax.set_ylabel("PhenoCycler pixel value (" + name_p + ")", fontsize=20)
        ax.tick_params(axis="x", labelsize=18)
        ax.tick_params(axis="y", labelsize=18)        

        filename_cor_in_tissue_png = EVALDIR_DAT + "/" + sid + "-" + name_v + "-" + name_p + "_cor_in_tissue.png"
        filename_cor_in_tissue_png_rel = "dat" + "/" + sid + "-" + name_v + "-" + name_p + "_cor_in_tissue.png"    
        fig.savefig(filename_cor_in_tissue_png, bbox_inches="tight", pad_inches=0.05, dpi=600, transparent=True)
        print("<img src=" + "\"" + filename_cor_in_tissue_png_rel + "\"" + " width=\"" + ww2 + "\"" + " border=\"0\"" + ">", file=fp_html)

        
        # Visium x PhenoCycler correlation (nCount_Spatial)

        plt.rcParams['font.family'] = "Arial"        
        fig, ax = plt.subplots(figsize=(9,9))

        ax.set_box_aspect (1)

        ddv = 0.05*(max(vals_v) - min(vals_v))
        ax.set_xlim(min(vals_v)-ddv, max(vals_v)+ddv)

        ddp = 0.05*(max(vals_p) - min(vals_p))
        ax.set_ylim(0-ddp/2, max(vals_p)+ddp)        

        mappable0 = ax.scatter(vals_v, vals_p, color=(1.0, 1.0, 1.0, 0.0), s=20, lw=1, ec=vals_ec_nc, cmap=cm, vmin=vmin_nc, vmax=vmax_nc)
        mappable1 = ax.scatter([-10,-10], [0,0], c=[vmin_nc, vmax_nc], s=50, cmap=cm) # dummy

        cb = fig.colorbar(mappable1, ax=ax, shrink=0.5, extend="both", extendrect=True, ticks=np.arange(2000, 20000+2500, 2500))
        cb.set_label("nCount_Spatial", loc="center")

        ax.set_xlabel("Visium SCTransform residual (" + name_v + ")", fontsize=20)
        ax.set_ylabel("PhenoCycler pixel value (" + name_p + ")", fontsize=20)
        ax.tick_params(axis="x", labelsize=18)
        ax.tick_params(axis="y", labelsize=18)        

        filename_cor_nCountSpatial_png = EVALDIR_DAT + "/" + sid + "-" + name_v + "-" + name_p + "_cor_nCountSpatial.png"
        filename_cor_nCountSpatial_png_rel = "dat" + "/" + sid + "-" + name_v + "-" + name_p + "_cor_nCountSpatial.png"    
        fig.savefig(filename_cor_nCountSpatial_png, bbox_inches="tight", pad_inches=0.05, dpi=600, transparent=True)
        print("<img src=" + "\"" + filename_cor_nCountSpatial_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

       
        # Visium x PhenoCycler correlation (nFeature_Spatial)

        plt.rcParams['font.family'] = "Arial"        
        fig, ax = plt.subplots(figsize=(9,9))

        ax.set_box_aspect (1)

        ddv = 0.05*(max(vals_v) - min(vals_v))
        ax.set_xlim(min(vals_v)-ddv, max(vals_v)+ddv)

        ddp = 0.05*(max(vals_p) - min(vals_p))
        ax.set_ylim(0-ddp/2, max(vals_p)+ddp)        

        mappable0 = ax.scatter(vals_v, vals_p, color=(1.0, 1.0, 1.0, 0.0), s=20, lw=1, ec=vals_ec_nf, cmap=cm, vmin=vmin_nf, vmax=vmax_nf)
        mappable1 = ax.scatter([-10,-10], [0,0], c=[vmin_nf, vmax_nf], s=50, cmap=cm) # dummy

        cb = fig.colorbar(mappable1, ax=ax, shrink=0.5, extend="both", extendrect=True, ticks=np.arange(2000, 8000+500, 500))
        cb.set_label("nFeature_Spatial", loc="center")

        ax.set_xlabel("Visium SCTransform residual (" + name_v + ")", fontsize=20)
        ax.set_ylabel("PhenoCycler pixel value (" + name_p + ")", fontsize=20)
        ax.tick_params(axis="x", labelsize=18)
        ax.tick_params(axis="y", labelsize=18)        

        filename_cor_nFeatureSpatial_png = EVALDIR_DAT + "/" + sid + "-" + name_v + "-" + name_p + "_cor_nFeatureSpatial.png"
        filename_cor_nFeatureSpatial_png_rel = "dat" + "/" + sid + "-" + name_v + "-" + name_p + "_cor_nFeatureSpatial.png"    
        fig.savefig(filename_cor_nFeatureSpatial_png, bbox_inches="tight", pad_inches=0.05, dpi=600, transparent=True)
        #print("<img src=" + "\"" + filename_cor_nFeatureSpatial_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
        
        # Visium (Spatial@data) x PhenoCycler correlation (nCount_Spatial)

        plt.rcParams['font.family'] = "Arial"        
        fig, ax = plt.subplots(figsize=(9,9))

        ax.set_box_aspect (1)

        ddd = 0.05*(max(vals_d) - min(vals_d))
        ax.set_xlim(min(vals_d)-ddd, max(vals_d)+ddd)

        ddp = 0.05*(max(vals_p) - min(vals_p))
        ax.set_ylim(0-ddp/2, max(vals_p)+ddp)        

        mappable0 = ax.scatter(vals_d, vals_p, color=(1.0, 1.0, 1.0, 0.0), s=20, lw=1, ec=vals_ec_nc, cmap=cm, vmin=vmin_nc, vmax=vmax_nc)
        mappable1 = ax.scatter([-10,-10], [0,0], c=[vmin_nc, vmax_nc], s=50, cmap=cm) # dummy

        cb = fig.colorbar(mappable1, ax=ax, shrink=0.5, extend="both", extendrect=True, ticks=np.arange(2000, 20000+2500, 2500))
        cb.set_label("nCount_Spatial", loc="center")

        ax.set_xlabel("Visium log-normalized value (" + name_v + ")", fontsize=20)
        ax.set_ylabel("PhenoCycler pixel value (" + name_p + ")", fontsize=20)
        ax.tick_params(axis="x", labelsize=18)
        ax.tick_params(axis="y", labelsize=18)        

        filename_cor_nCountSpatial_lognorm_png = EVALDIR_DAT + "/" + sid + "-" + name_v + "-" + name_p + "_cor_nCountSpatial_lognorm.png"
        filename_cor_nCountSpatial_lognorm_png_rel = "dat" + "/" + sid + "-" + name_v + "-" + name_p + "_cor_nCountSpatial_lognorm.png"    
        fig.savefig(filename_cor_nCountSpatial_lognorm_png, bbox_inches="tight", pad_inches=0.05, dpi=600, transparent=True)
        print("<img src=" + "\"" + filename_cor_nCountSpatial_lognorm_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
        
        fp_html.flush()


    #####

    # Visium

    # SCT@counts
    filename_sct_cc = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-counts.txt"
    
    # SCT@data
    filename_sct_dd = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-data.txt"
    
    # SCT@scale.data
    filename_sct_sd = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-scale_data.txt"
    
    # Spatial@counts
    filename_spatial_cc = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-Spatial-counts.txt"
    
    # Spatial@data
    filename_spatial_dd = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-Spatial-data.txt"

    df_sct_sd = pd.read_csv(filename_sct_sd, index_col=0)
    df_spatial_cc = pd.read_csv(filename_spatial_cc, index_col=0)
    df_spatial_dd = pd.read_csv(filename_spatial_dd, index_col=0)
    
    filename_ncount_txt = NORMDIR_DAT + "/" + PREFIX_v + "-nCount_Spatial.txt"
    filename_nfeature_txt = NORMDIR_DAT + "/" + PREFIX_v + "-nFeature_Spatial.txt"
    
    df_nc = pd.read_csv(filename_ncount_txt, index_col=0)
    df_nf = pd.read_csv(filename_nfeature_txt, index_col=0)    

    #####

    DO_COUNTCELLS = True

    if DO_COUNTCELLS:
        
        print("FILENAME_SEG", filename_seg)

        #####
        
        cellpos = CountCell(filename_seg, sid)

        name_p = "DAPI"

        count, pid, filename_PAGE_gray_nega = PAGEID[name_p]
        print(filename_PAGE_gray_nega)
        img_pp = cv2.imread(filename_PAGE_gray_nega, cv2.IMREAD_GRAYSCALE) # nega
        #img_pp_01 = cv2.resize(img_pp, None, None, 0.1, 0.1, cv2.INTER_NEAREST) # nega
        
        img_pp_cell = np.ones_like(img_pp, dtype=np.uint8)
        img_pp_cell = img_pp_cell * 255
        img_pp_cell = cv2.cvtColor(img_pp_cell, cv2.COLOR_GRAY2BGRA)

        img_pp_cell_circle = img_pp_cell.copy()


        CELL_POS = {}
        for (xc, yc, area) in cellpos:

            xi = int(xc)
            yi = int(yc)

            img_pp_cell[yi, xi] = [255, 0, 0, 255]
            
            cv2.circle(img_pp_cell_circle, (xi, yi), 3, [255, 0, 0, 255], thickness=-1, lineType=cv2.LINE_AA)

            if not xi in CELL_POS:
                CELL_POS[xi] = {}

            CELL_POS[xi][yi] = area
            

        #####

        #img_pp_cell_01 = cv2.resize(img_pp_cell, None, None, 0.1, 0.1, cv2.INTER_NEAREST)
        img_pp_cell_01 = cv2.resize(img_pp_cell_circle, None, None, 0.1, 0.1, cv2.INTER_NEAREST)        
        img_pp_cell_01 = bg_transp(img_pp_cell_01)

        
        #####

        # cell count in TSV format
        filename_pp_cell_txt = EVALDIR_DAT + "/" + sid + "_CellCount.txt"
        filename_pp_cell_txt_rel = "dat" + "/" + sid + "_CellCount.txt"

        fp_cc = open(filename_pp_cell_txt, "wt")
        print("Barcode" + "\t" + "Number of Cells", file=fp_cc)
                
        _, _, filename_PAGE_gray = PAGEID_pos[name_p]        
        
        img_pp_pos = cv2.imread(filename_PAGE_gray, cv2.IMREAD_GRAYSCALE) # pos
        
        img_pp_pos_bgra = cv2.cvtColor(img_pp_pos, cv2.COLOR_GRAY2BGRA)
        img_pp_pos_bgra_01 = cv2.resize(img_pp_pos_bgra, None, None, 0.1, 0.1, cv2.INTER_NEAREST)

        # base = tissue
        img_pp_pos_bgra_base_01 = np.ones_like(img_pp_pos_bgra_01)
        img_pp_pos_bgra_base_01 = img_pp_pos_bgra_base_01 * 255
        img_pp_pos_bgra_base_01[np.where(img_phenocycler_mask_wh_01 == 255)] = [240, 240, 240, 255]

        img_pp_pos_bgra_base = np.ones_like(img_pp_pos_bgra)
        img_pp_pos_bgra_base = img_pp_pos_bgra_base * 255
        img_pp_pos_bgra_base[np.where(img_phenocycler_mask_wh == 255)] = [240, 240, 240, 255]
        

        # base img
        img_pp_pos_bgra_base_cellplot_01 = img_pp_pos_bgra_base_01.copy()
        img_pp_pos_bgra_base_cellplot = img_pp_pos_bgra_base.copy()        

        # cell positions
        img_pp_pos_bgra_base_cellplot_01[np.where(~np.all(img_pp_cell_01[:,:,:3] == 255, axis=-1))] = [255, 0, 0, 255]
        img_pp_pos_bgra_base_cellplot[np.where(~np.all(img_pp_cell[:,:,:3] == 255, axis=-1))] = [255, 0, 0, 255]

        m_count_th = math.pi * rp_s * rp_s * 0.0
        
        
        # spot
        df = df_sct_sd
        count = 0
        vals_c = []

        hp,wp,_ = img_pp_cell.shape
        
        for bc in df.columns.values:
            
            (xv_f, yv_f, in_tissue, ci, ri) = SPOT[bc]

            if in_tissue == 0:
                continue

            xp_f, yp_f = f_v2p(MATRIX, scale, xv_f, yv_f)

            m_count = 0
            c_count = 0
            for dx in range(-1*rp_s, rp_s+1):
                for dy in range(-1*rp_s, rp_s+1):

                    #print("DX//DY", dx, dy)

                    x = xp_f + dx
                    y = yp_f + dy

                    v = img_pp_cell[y, x]
                    m = img_phenocycler_mask_wh[y, x]

                    d = math.sqrt( (x-xp_f)*(x-xp_f) + (y-yp_f)*(y-yp_f) )

                    if d > rp_s:
                        continue

                    if m > 0:
                        m_count += 1

                        if x in CELL_POS:
                            if y in CELL_POS[x]:
                                c_count += 1
                            

                        #if ~np.all(v[:3] == 255): # if not [255, 255, 255, 255]
                        #    c_count +=  1

            xp_f_01 = int(0.1*xp_f)
            yp_f_01 = int(0.1*yp_f)
            rp_s_01 = int(0.1*rp_s)

            if m_count <= m_count_th or in_tissue == 0:
                continue
            
            cv2.circle(img_pp_pos_bgra_base_cellplot_01, (xp_f_01, yp_f_01), rp_s_01, [0, 0, 255, 255], thickness=1, lineType=cv2.LINE_AA)

            #if c_count >= 80:
            #    cv2.circle(img_pp_pos_bgra_base_cellplot, (xp_f, yp_f), rp_s, [0, 0, 255, 255], thickness=1, lineType=cv2.LINE_AA)
            
            cv2.circle(img_pp_pos_bgra_base_cellplot, (xp_f, yp_f), rp_s, [0, 0, 255, 255], thickness=1, lineType=cv2.LINE_AA)            


            #print(c_count)

            print(bc + "\t" + "{:d}".format(c_count), file=fp_cc)

            #print(count, bc, c_count)

            vals_c.append(c_count)

            count += 1

        fp_cc.flush()
        fp_cc.close()
            
        # boundary
        img_pp_pos_bgra_base_cellplot_01[np.where(img_phenocycler_boundary_wh_01 == 255)] = [0, 0, 0, 255]
        img_pp_pos_bgra_base_cellplot[np.where(img_phenocycler_boundary_wh == 255)] = [0, 0, 0, 255]        
        
        img_pp_pos_bgra_base_cellplot_01 = bg_transp(img_pp_pos_bgra_base_cellplot_01)
        img_pp_pos_bgra_base_cellplot = bg_transp(img_pp_pos_bgra_base_cellplot)

        
        print("<h2>", file=fp_html)
        print("Number of cells / spot", file=fp_html)
        print("</h2>", file=fp_html)

        filename_pp_cell_png = EVALDIR_DAT + "/" + sid + "_CellCount_01.png"
        filename_pp_cell_png_rel = "dat" + "/" + sid + "_CellCount_01.png"

        cv2.imwrite(filename_pp_cell_png, img_pp_pos_bgra_base_cellplot_01)
        print("<img src=" + "\"" + filename_pp_cell_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

        # x10
        filename_pp_cell_png_10 = EVALDIR_DAT + "/" + sid + "_CellCount.png"                
        cv2.imwrite(filename_pp_cell_png_10, img_pp_pos_bgra_base_cellplot)

        
        # histgram

        plt.rcParams['font.family'] = "Arial"        
        fig, ax = plt.subplots(figsize=(9,9))

        ax.set_box_aspect (1)

        dw = 1        

        vals_c_max = max(vals_c) + 1
        vals_c_min = 0 #min(vals_c)

        bin_max = vals_c_max - vals_c_max%dw + dw
        bin_min = vals_c_min - vals_c_min%dw
        bin_num = int((bin_max - bin_min) / dw) + 1
        bins = np.linspace(bin_min, bin_max, bin_num)

        ax.hist(vals_c, bins=bins, ec="black")

        ax.set_xlabel("Number of cells / Spot", fontsize=20)
        ax.set_ylabel("Frequency of spots", fontsize=20)
        ax.tick_params(axis="x", labelsize=18)
        ax.tick_params(axis="y", labelsize=18)        

        filename_pp_cell_hist_png = EVALDIR_DAT + "/" + sid + "_CellCount_hist.png"
        filename_pp_cell_hist_png_rel = "dat" + "/" + sid + "_CellCount_hist.png"
        
        fig.savefig(filename_pp_cell_hist_png, bbox_inches="tight", pad_inches=0.05, dpi=600, transparent=True)
        print("<img src=" + "\"" + filename_pp_cell_hist_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
        
        print("<br>", file=fp_html)

        # TSV
        print("<a href=\"" + filename_pp_cell_txt_rel + "\">" + "Number of cells / spot" + "</a><br>", file=fp_html)

        fp_html.flush()
        
        
    #####
    # CORELATION
    
    name_v = "MS4A1" # "CD20"
    name_p = "CD20" # (09/P10)    
    PlotCor(df_sct_sd, df_spatial_cc, df_spatial_dd, name_v, name_p)

    name_v = "CD68" # 
    name_p = "CD68" # (08/P09)    
    PlotCor(df_sct_sd, df_spatial_cc, df_spatial_dd, name_v, name_p)

    name_v = "CDH1" # 
    name_p = "E-Cadherin" # 
    PlotCor(df_sct_sd, df_spatial_cc, df_spatial_dd, name_v, name_p)

    name_v = "TTF1" # 
    name_p = "TTF-1-NKX2-1" # 
    PlotCor(df_sct_sd, df_spatial_cc, df_spatial_dd, name_v, name_p)

    name_v = "ACTA2" # 
    name_p = "a-SMA-ACTA2" # 
    PlotCor(df_sct_sd, df_spatial_cc, df_spatial_dd, name_v, name_p)

    name_v = "CXCL13" # 
    name_p = "BCA1-CXCL13" # 
    PlotCor(df_sct_sd, df_spatial_cc, df_spatial_dd, name_v, name_p)

    name_v = "CD44" # 
    name_p = "CD44" # (00/P00)    
    PlotCor(df_sct_sd, df_spatial_cc, df_spatial_dd, name_v, name_p)
    
    name_v = "CD4" # 
    name_p = "CD4" # (00/P00)    
    PlotCor(df_sct_sd, df_spatial_cc, df_spatial_dd, name_v, name_p)

    name_v = "CTNNB1" # 
    name_p = "b-Catenin1" # (00/P00)    
    PlotCor(df_sct_sd, df_spatial_cc, df_spatial_dd, name_v, name_p)

    name_v = "ACTB" # 
    name_p = "Beta-actin" # (00/P00)    
    PlotCor(df_sct_sd, df_spatial_cc, df_spatial_dd, name_v, name_p)

    name_v = "HLA-A" # 
    name_p = "HLA-A" # (00/P00)    
    PlotCor(df_sct_sd, df_spatial_cc, df_spatial_dd, name_v, name_p)

    name_v = "LGALS3BP" # 
    name_p = "Mac2-Galectin3" # (00/P00)    
    PlotCor(df_sct_sd, df_spatial_cc, df_spatial_dd, name_v, name_p)

    name_v = "HLA-DPB1" # 
    name_p = "HLA-DPB1" # (00/P00)    
    PlotCor(df_sct_sd, df_spatial_cc, df_spatial_dd, name_v, name_p)

    #####
            
    print_html_footer(fp_html)            
        
    fp_html.close()
   

    
if __name__ == "__main__":

    main()


    
