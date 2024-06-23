#!/usr/bin/env python3

import sys
import re
import pathlib
import shutil
import math
import cv2

import matplotlib.pyplot as plt
import subprocess

from util.load_data import *
from util.html import *


def main():

    if len(sys.argv) < 2:
        print("sys.argv")
        sys.exit()

    a0 = sys.argv[0]
    pp = pathlib.Path(a0)
    SRCDIR = str(pp.parent)
        
    sid = sys.argv[1]

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

    # NORMDIR    
    NORMDIR = RESULTDIR + "/" + PREFIX_v + "-" + "NORM"
    NORMDIR_DAT = NORMDIR + "/" + "dat"

    filename_rds = NORMDIR_DAT + "/" + PREFIX_v + "-" + "SCT.rds"

    # CLSTDIR
    CLSTDIR = RESULTDIR + "/" + PREFIX_v + "-" + "CLST"
    pp = pathlib.Path(CLSTDIR)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)

    # CLSTDIR_DAT
    CLSTDIR_DAT = CLSTDIR + "/" + "dat"
    pp = pathlib.Path(CLSTDIR_DAT)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)

    #####
    # LOAD visium text files
    SCALEFACTORS, spot_diameter_fullres, tissue_hires_scalef, fiducial_diameter_fullres, SPOT, GRID_rc, GRID_cr, FID, FMP = \
        load_visium_text_files(0, MASKDIR_v_DAT)
    scale = tissue_hires_scalef

    
    #####
    command = SRCDIR + "/exec_clustering.r" + " " + sid + " " + filename_rds + " " + CLSTDIR_DAT 
    print(command)
    p = subprocess.run(command, shell=True, check=True)
    
    #####

    filename_html = CLSTDIR + "/" + PREFIX_v + "-" + "CLST" + ".html"    
    print(filename_html)
    fp_html = open(filename_html, "wt")

    ww = "200px"
    
    print_html_header(fp_html)

    print("<h1>", file=fp_html)
    print("Clustering and dimensionality reduction of Visium data" + " &mdash; " + sid, file=fp_html)
    print("</h1>", file=fp_html)

    #####
    
    print("<h2>", file=fp_html)
    print("Clustering (SpatialDimPlot)", file=fp_html)
    print("</h2>", file=fp_html)

    # clip.range
    print("<h3>", file=fp_html)
    print("SCTransform with clip.range=c(-20,50) option (slice)", file=fp_html)
    print("</h3>", file=fp_html)
    
    filename_png = CLSTDIR_DAT + "/" + PREFIX_v + "-SpatialDimPlot.png"
    filename_png_rel = "dat" + "/" + PREFIX_v + "-SpatialDimPlot.png"
    print("<img src=" + "\"" + filename_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    print("<br>", file=fp_html)
    
    filename_slice_txt = CLSTDIR_DAT + "/" + PREFIX_v + "-CLST-FindClusters-slice.txt"
    filename_slice_txt_rel = "dat" + "/" + PREFIX_v + "-CLST-FindClusters-slice.txt"
    print("<a href=\"" + filename_slice_txt_rel + "\">" + "str() of Seurat object" + "</a>", file=fp_html)
    
    print("<br>", file=fp_html)
    
    filename_seurat_clusters_txt = CLSTDIR_DAT + "/" + PREFIX_v + "-SCT_seurat_clusters.txt"
    filename_seurat_clusters_txt_rel = "dat" + "/" + PREFIX_v + "-SCT_seurat_clusters.txt"
    print("<a href=\"" + filename_seurat_clusters_txt_rel + "\">" + "seurat_clusters (SCT)" + "</a>", file=fp_html)
    
    print("<br>", file=fp_html)

    # Default options
    print("<h3>", file=fp_html)
    print("SCTransform with default options (return.only.var.genes=FALSE) (slice1)", file=fp_html)
    print("</h3>", file=fp_html)
    
    filename_png = CLSTDIR_DAT + "/" + PREFIX_v + "-SpatialDimPlot1.png"
    filename_png_rel = "dat" + "/" + PREFIX_v + "-SpatialDimPlot1.png"
    print("<img src=" + "\"" + filename_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    print("<br>", file=fp_html)
    
    filename_slice1_txt = CLSTDIR_DAT + "/" + PREFIX_v + "-CLST-FindClusters-slice1.txt"
    filename_slice1_txt_rel = "dat" + "/" + PREFIX_v + "-CLST-FindClusters-slice1.txt"
    print("<a href=\"" + filename_slice1_txt_rel + "\">" + "str() of Seurat object" + "</a>", file=fp_html)
    
    print("<br>", file=fp_html)
    
    filename_seurat_clusters1_txt = CLSTDIR_DAT + "/" + PREFIX_v + "-SCT_seurat_clusters1.txt"
    filename_seurat_clusters1_txt_rel = "dat" + "/" + PREFIX_v + "-SCT_seurat_clusters1.txt"
    print("<a href=\"" + filename_seurat_clusters1_txt_rel + "\">" + "seurat_clusters (SCT)" + "</a>", file=fp_html)
    
    print("<br>", file=fp_html)
    
    #####
    
    print("<h2>", file=fp_html)
    print("Clustering (on hires and fullres)", file=fp_html)    
    print("</h2>", file=fp_html)
   
    filename_seurat_clusters_rgb_txt = CLSTDIR_DAT + "/" + PREFIX_v + "-SCT_seurat_clusters_RGB.txt"
    
    def drawSpotsCluster(filename_seurat_clusters_txt, filename_seurat_clusters_rgb_txt):

        filename_visium_hi = MASKDIR_v_DAT + "/" + PREFIX_v + ".png"    
        filename_visium_ful = MASKDIR_v_DAT + "/" + PREFIX_v + "_ful.png"

        filename_visium_hi_spots = MASKDIR_v_DAT + "/" + PREFIX_v + "_spots.png"        

        img_visium_hi = cv2.imread(filename_visium_hi, cv2.IMREAD_UNCHANGED)
        img_visium_ful = cv2.imread(filename_visium_ful, cv2.IMREAD_UNCHANGED)

        img_visium_hi_spots = img_visium_hi.copy()

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

            print(rgb, r, g, b)

            RGB[cid] = (r, g, b)

        fp_cls = open(filename_seurat_clusters_txt)
        line = fp_cls.readline()

        CLST = {}

        for line in fp_cls:

            line = line.rstrip()
            vals = line.split(",")

            bc = vals[0]
            cid = int(vals[1])

            #print(bc, cid)

            CLST[bc] = cid

        cid_min = min(CLST.values())            
        cid_max = max(CLST.values())

        # diameters/radius of spots

        xo, yo = 0, 0
        xds, yds = spot_diameter_fullres, 0
        xdf, ydf = fiducial_diameter_fullres, 0    

        dp_s_ful = math.sqrt((xds-xo)*(xds-xo) + (yds-yo)*(yds-yo))
        rp_s_ful = int(dp_s_ful/2)

        rp_s_hi = int(rp_s_ful * tissue_hires_scalef)

        dp_f_ful = math.sqrt((xdf-xo)*(xdf-xo) + (ydf-yo)*(ydf-yo))
        rp_f_ful = int(dp_f_ful/2)

        rp_f_hi = int(rp_f_ful * tissue_hires_scalef)    

        for bc, (xv_ful, yv_ful, in_tissue, ci, ri) in SPOT.items():

            if bc in CLST:
                cid = CLST[bc]
            else:
                cid = -1

            xv_hi = int(xv_ful * tissue_hires_scalef)
            yv_hi = int(yv_ful * tissue_hires_scalef)

            if in_tissue == 1:

                cv2.circle(img_visium_hi_spots, (xv_hi, yv_hi), rp_s_hi, [0, 0, 255, 255], thickness=1)

                if not cid == -1:

                    (r, g, b) = RGB[cid]
                    cv2.circle(img_visium_hi, (xv_hi, yv_hi), rp_s_hi, [b, g, r, 255], thickness=-1)
                    cv2.circle(img_visium_ful, (xv_ful, yv_ful), rp_s_ful, [b, g, r, 255], thickness=-1)

            else:
                cv2.circle(img_visium_hi_spots, (xv_hi, yv_hi), rp_s_hi, [255, 0, 0, 255], thickness=1)                        
                cv2.circle(img_visium_hi, (xv_hi, yv_hi), rp_s_hi, [255, 0, 0, 255], thickness=1)            
                cv2.circle(img_visium_ful, (xv_ful, yv_ful), rp_s_ful, [255, 0, 0, 255], thickness=5)
                

        return img_visium_hi_spots, img_visium_hi, img_visium_ful, cid_max, RGB


    print("<h3>", file=fp_html)
    print("SCTransform with clip.range=c(-20,50) option (slice)", file=fp_html)
    print("</h3>", file=fp_html)

    img_visium_hi_spots, img_visium_hi, img_visium_ful, cid_max, RGB = drawSpotsCluster(filename_seurat_clusters_txt, filename_seurat_clusters_rgb_txt)

    # legend
    cluster_num = cid_max+1
    cmap = plt.get_cmap("rainbow", cluster_num)    
    colors = []
    for k in range(0, cluster_num):
        #c = cmap(k / cluster_num)
        (b, g, r) = RGB[k]
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
    filename_vv_cluster_legend = CLSTDIR_DAT + "/" + PREFIX_v + "-cluster_legend.png"            
    fig.savefig(filename_vv_cluster_legend, bbox_inches=bbox, pad_inches=0.05, dpi=600, transparent=True)
    plt.clf()    
    plt.close()

    filename_png_spots = CLSTDIR_DAT + "/" + PREFIX_v + "-SpatialDimPlot_hires_spots.png"
    filename_png_spots_rel = "dat" + "/" + PREFIX_v + "-SpatialDimPlot_hires_spots.png"
    cv2.imwrite(filename_png_spots, img_visium_hi_spots)
    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_png_spots_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    filename_png = CLSTDIR_DAT + "/" + PREFIX_v + "-SpatialDimPlot_hires.png"
    filename_png_rel = "dat" + "/" + PREFIX_v + "-SpatialDimPlot_hires.png"
    cv2.imwrite(filename_png, img_visium_hi)
    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    filename_png_ful = CLSTDIR_DAT + "/" + PREFIX_v + "-SpatialDimPlot_fullres.png"
    filename_png_ful_rel = "dat" + "/" + PREFIX_v + "-SpatialDimPlot_fullres.png"
    cv2.imwrite(filename_png_ful, img_visium_ful)        
    #print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_png_ful_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    filename_vv_cluster_legend = CLSTDIR_DAT + "/" + PREFIX_v + "-cluster_legend.png"    
    filename_vv_cluster_legend_rel = "dat" + "/" + PREFIX_v + "-cluster_legend.png"                
    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vv_cluster_legend_rel + "\"" + " width=\"" + "30px" + "\"" + " border=\"1\"" + ">", file=fp_html)

    
    print("<h3>", file=fp_html)
    print("SCTransform with default options (return.only.var.genes=FALSE) (slice1)", file=fp_html)
    print("</h3>", file=fp_html)

    filename_seurat_clusters1_txt = CLSTDIR_DAT + "/" + PREFIX_v + "-SCT_seurat_clusters1.txt"    
    filename_seurat_clusters1_rgb_txt = CLSTDIR_DAT + "/" + PREFIX_v + "-SCT_seurat_clusters1_RGB.txt"

    img_visium_hi_spots1, img_visium_hi1, img_visium_ful1, cid_max1, RGB1 = drawSpotsCluster(filename_seurat_clusters1_txt, filename_seurat_clusters1_rgb_txt)

    # legend
    cluster_num1 = cid_max1+1
    cmap = plt.get_cmap("rainbow", cluster_num1)    
    colors1 = []
    for k in range(0, cluster_num1):
        #c = cmap(k / cluster_num1)
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
    filename_vv_cluster_legend1 = CLSTDIR_DAT + "/" + PREFIX_v + "-cluster_legend1.png"            
    fig.savefig(filename_vv_cluster_legend1, bbox_inches=bbox, pad_inches=0.05, dpi=600, transparent=True)
    plt.clf()    
    plt.close()
    

    filename_png_spots1 = CLSTDIR_DAT + "/" + PREFIX_v + "-SpatialDimPlot_hires_spots1.png"
    filename_png_spots1_rel = "dat" + "/" + PREFIX_v + "-SpatialDimPlot_hires_spots1.png"
    cv2.imwrite(filename_png_spots1, img_visium_hi_spots1)
    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_png_spots1_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    filename_png1 = CLSTDIR_DAT + "/" + PREFIX_v + "-SpatialDimPlot_hires1.png"
    filename_png1_rel = "dat" + "/" + PREFIX_v + "-SpatialDimPlot_hires1.png"
    cv2.imwrite(filename_png1, img_visium_hi1)
    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_png1_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    filename_png_ful1 = CLSTDIR_DAT + "/" + PREFIX_v + "-SpatialDimPlot_fullres1.png"
    filename_png_ful1_rel = "dat" + "/" + PREFIX_v + "-SpatialDimPlot_fullres1.png"
    cv2.imwrite(filename_png_ful1, img_visium_ful1)        
    #print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_png_ful1_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    filename_vv_cluster_legend1 = CLSTDIR_DAT + "/" + PREFIX_v + "-cluster_legend1.png"    
    filename_vv_cluster_legend1_rel = "dat" + "/" + PREFIX_v + "-cluster_legend1.png"                
    print("<img style=\"vertical-align: middle;\" src=" + "\"" + filename_vv_cluster_legend1_rel + "\"" + " width=\"" + "30px" + "\"" + " border=\"1\"" + ">", file=fp_html)
    
    #####    

    print("<h2>", file=fp_html)
    print("Dimensionality reduction (DimPlot)", file=fp_html)
    print("</h2>", file=fp_html)
    
    filename_png = CLSTDIR_DAT + "/" + PREFIX_v + "-DimPlot_UMAP.png"
    filename_png_rel = "dat" + "/" + PREFIX_v + "-DimPlot_UMAP.png"
    print("<img src=" + "\"" + filename_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    #####
        
    print("<h2>", file=fp_html)
    print("Spatially variable features (FindMarkers)", file=fp_html)
    print("</h2>", file=fp_html)

    filename_degs_txt = CLSTDIR_DAT + "/" + PREFIX_v + "-SpatialFeaturePlot_FM_filenames.txt"

    fp_txt = open(filename_degs_txt, "rt")
    line = fp_txt.readline()

    cid_prev = -1

    for line in fp_txt:
        line = line.rstrip()
        vals = line.split("\t")

        cid = int(vals[0])
        nid = int(vals[1])
        filename = vals[2]
        gene = vals[3]
        pval = float(vals[4])
        pval_adj = float(vals[5])
        log2fc = float(vals[6])

        if cid > cid_prev:

            print("<h3>", file=fp_html)
            print("Cluster {:d}".format(cid), file=fp_html)
            print("</h3>", file=fp_html)

            cid_prev = cid

        filename_png = CLSTDIR_DAT + "/" + filename
        filename_png_rel = "dat" + "/" + filename
        
        print("<figure style=\"display:inline-block; margin:0px 0px 0px 0px;\">", file=fp_html)
        print("<figcaption>", file=fp_html)
        print("<b>{:}</b>".format(gene), file=fp_html)            
        print("</figcaption>", file=fp_html)                        
        print("<img style=\"margin:0px 0px 0px 0px;\" src=" + "\"" + filename_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
        print("<figcaption>", file=fp_html)

        print("p: {:.2e}<br>p_adj: {:.2e}<br>Log2FC: {:.2f}".format(pval, pval_adj, log2fc), file=fp_html)    
        
        print("</figcaption>", file=fp_html)                
        print("</figure>", file=fp_html)    
    
    print("</body>", file=fp_html)
    print("</html>", file=fp_html)
        
    fp_html.close()
    

#####

if __name__ == "__main__":

    main()
