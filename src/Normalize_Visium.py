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

import subprocess

from util.html import *

def main():
    
    if len(sys.argv) < 2:
        print("sys.argv")
        sys.exit()

    a0 = sys.argv[0]
    pp = pathlib.Path(a0)
    SRCDIR = str(pp.parent)
        
    dirname_spaceranger = sys.argv[1]
    sid = sys.argv[2]

    # directories
    
    WORKDIR = "."
    RESULTDIR = WORKDIR + "/" + sid

    PREFIX = sid   
    
    # Visium
    PREFIX_v = sid + "-" + "Visium"

    # NORMDIR
    NORMDIR = RESULTDIR + "/" + PREFIX_v + "-" + "NORM"
    print("NORMDIR", NORMDIR)
    pp = pathlib.Path(NORMDIR)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)
    NORMDIR_DAT = NORMDIR + "/" + "dat"
    print("NORMDIR_DAT", NORMDIR_DAT)
    pp = pathlib.Path(NORMDIR_DAT)
    if not pp.exists():
        pp.mkdir(exist_ok=True, parents=True)

    # Directory of SpaceRanger result containing "spatial" subdirectory
    SPACERANGERDIR = dirname_spaceranger
    
    command = SRCDIR + "/exec_preprocess.r" + " " + sid + " " + SPACERANGERDIR + " " + NORMDIR 
    print(command)
    p = subprocess.run(command, shell=True, check=True)

    #####
        
    filename_html = NORMDIR + "/" + PREFIX_v + "-" + "NORM" + ".html"
    print(filename_html)
    fp_html = open(filename_html, "wt")

    ww = "200px"

    print_html_header(fp_html)    

    print("<h1>", file=fp_html)
    print("QC/Filteration and normalization of Visium data" + " &mdash; " + sid, file=fp_html)
    print("</h1>", file=fp_html)

    #####

    # slice initial
    print("<h2>", file=fp_html)
    print("Original data", file=fp_html)
    print("</h2>", file=fp_html)
    
    filename_slice_txt = NORMDIR_DAT + "/" + PREFIX_v + "-000-slice.txt"
    filename_slice_txt_rel = "dat" + "/" + PREFIX_v + "-000-slice.txt"
    print("<a href=\"" + filename_slice_txt_rel + "\">" + "str() of Seurat object (initial)" + "</a>", file=fp_html)

    #####
    
    # slice after PFS
    print("<h3>", file=fp_html)
    print("PercentageFeatureSet", file=fp_html)
    print("</h3>", file=fp_html)

    print("<ul>", file=fp_html)
    print("  <li>percent_mito</li>", file=fp_html)
    print("  <li>percent_hb</li>", file=fp_html)
    print("  <li>percent_ribo</li>", file=fp_html)        
    print("</ul>", file=fp_html)
    
    print("In percentage (ranging 0 to 100 %), not rational fraction (0 to 1.0).<br>", file=fp_html)
    
    filename_slice_txt = NORMDIR_DAT + "/" + PREFIX_v + "-PFS-slice.txt"
    filename_slice_txt_rel = "dat" + "/" + PREFIX_v + "-PFS-slice.txt"
    print("<a href=\"" + filename_slice_txt_rel + "\">" + "str() of Seurat object (after PercentFeatureSet() for mitochondrial, hemoglobin, and ribosomal genes)" + "</a><br>", file=fp_html)

    # VlnPlot(1)
    print("<h3>", file=fp_html)
    print("VlnPlot (Original data)", file=fp_html)
    print("</h3>", file=fp_html)
    
    # nCount_Spatial
    filename_ncount_vln_png = NORMDIR_DAT + "/" + PREFIX_v + "-VLN-nCount_Spatial.png"
    filename_ncount_vln_png_rel = "dat" + "/" + PREFIX_v + "-VLN-nCount_Spatial.png"
    print("<img src=" + "\"" + filename_ncount_vln_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # nFeature_Spatial
    filename_nfeature_vln_png = NORMDIR_DAT + "/" + PREFIX_v + "-VLN-nFeature_Spatial.png"
    filename_nfeature_vln_png_rel = "dat" + "/" + PREFIX_v + "-VLN-nFeature_Spatial.png"    
    print("<img src=" + "\"" + filename_nfeature_vln_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # percent_mito
    filename_percent_mito_vln_png = NORMDIR_DAT + "/" + PREFIX_v + "-VLN-percent_mito.png"
    filename_percent_mito_vln_png_rel = "dat" + "/" + PREFIX_v + "-VLN-percent_mito.png"
    print("<img src=" + "\"" + filename_percent_mito_vln_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
    
    # percent_hb
    filename_percent_hb_vln_png = NORMDIR_DAT + "/" + PREFIX_v + "-VLN-percent_hb.png"
    filename_percent_hb_vln_png_rel = "dat" + "/" + PREFIX_v + "-VLN-percent_hb.png"
    print("<img src=" + "\"" + filename_percent_hb_vln_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # percent_ribo
    filename_percent_ribo_vln_png = NORMDIR_DAT + "/" + PREFIX_v + "-VLN-percent_ribo.png"
    filename_percent_ribo_vln_png_rel = "dat" + "/" + PREFIX_v + "-VLN-percent_ribo.png"    
    print("<img src=" + "\"" + filename_percent_ribo_vln_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    #####

    print("<br>", file=fp_html)

    # nCount_Spatial, txt
    filename_ncount_txt = NORMDIR_DAT + "/" + PREFIX_v + "-nCount_Spatial.txt"
    filename_ncount_txt_rel = "dat" + "/" + PREFIX_v + "-nCount_Spatial.txt"
    print("<a href=\"" + filename_ncount_txt_rel + "\">" + "nCount_Spatial" + "</a>", file=fp_html)

    print("<br>", file=fp_html)    

    # nFeature_Spatial, txt
    filename_nfeature_txt = NORMDIR_DAT + "/" + PREFIX_v + "-nFeature_Spatial.txt"
    filename_nfeature_txt_rel = "dat" + "/" + PREFIX_v + "-nFeature_Spatial.txt"
    print("<a href=\"" + filename_nfeature_txt_rel + "\">" + "nFeature_Spatial" + "</a>", file=fp_html)
    


    # FeatureScatter(1)
    print("<h3>", file=fp_html)
    print("FeatureScatter (Original data)", file=fp_html)
    print("</h3>", file=fp_html)
   
    # nCount_Spatial vs nFeature_Spatial
    filename_fss_png = NORMDIR_DAT + "/" + PREFIX_v + "-FSS-nCount_Spatial-nFeature_Spatial.png"
    filename_fss_png_rel = "dat" + "/" + PREFIX_v + "-FSS-nCount_Spatial-nFeature_Spatial.png"
    print("<img src=" + "\"" + filename_fss_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
    
    
    # SpatialFeaturePlot(1)
    print("<h3>", file=fp_html)
    print("SpatialFeaturePlot (Original data)", file=fp_html)
    print("</h3>", file=fp_html)

    # nCount_Spatial
    filename_ncount_sft_png = NORMDIR_DAT + "/" + PREFIX_v + "-SFT-nCount_Spatial.png"
    filename_ncount_sft_png_rel = "dat" + "/" + PREFIX_v + "-SFT-nCount_Spatial.png"    
    print("<img src=" + "\"" + filename_ncount_sft_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # nFeature_Spatial
    filename_nfeature_sft_png = NORMDIR_DAT + "/" + PREFIX_v + "-SFT-nFeature_Spatial.png"
    filename_nfeature_sft_png_rel = "dat" + "/" + PREFIX_v + "-SFT-nFeature_Spatial.png"
    print("<img src=" + "\"" + filename_nfeature_sft_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
    
    # percent_mito
    filename_percent_mito_sft_png = NORMDIR_DAT + "/" + PREFIX_v + "-SFT-percent_mito.png"
    filename_percent_mito_sft_png_rel = "dat" + "/" + PREFIX_v + "-SFT-percent_mito.png"
    print("<img src=" + "\"" + filename_percent_mito_sft_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # percent_hb
    filename_percent_hb_sft_png = NORMDIR_DAT + "/" + PREFIX_v + "-SFT-percent_hb.png"
    filename_percent_hb_sft_png_rel = "dat" + "/" + PREFIX_v + "-SFT-percent_hb.png"    
    print("<img src=" + "\"" + filename_percent_hb_sft_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
    
    # percent_ribo
    filename_percent_ribo_sft_png = NORMDIR_DAT + "/" + PREFIX_v + "-SFT-percent_ribo.png"
    filename_percent_ribo_sft_png_rel = "dat" + "/" + PREFIX_v + "-SFT-percent_ribo.png"    
    print("<img src=" + "\"" + filename_percent_ribo_sft_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    
    #####
    
    # slice after QC
    print("<h2>", file=fp_html)
    print("QC/Filteration", file=fp_html)
    print("</h2>", file=fp_html)

    print("slice = slice[, slice$nCount_Spatial > 500 & slice$nFeature_Spatial > 1000 & slice$nFeature_Spatial < 35000 & slice$percent_hb < 50]<br>", file=fp_html)
    
    filename_qc_slice_txt = NORMDIR_DAT + "/" + PREFIX_v + "-QC-slice.txt"
    filename_qc_slice_txt_rel = "dat" + "/" + PREFIX_v + "-QC-slice.txt"
    print("<a href=\"" + filename_qc_slice_txt_rel + "\">" + "str() of Seurat object (after QC/filteration)" + "</a><br>", file=fp_html)


    # VlnPlot(2)
    print("<h3>", file=fp_html)
    print("VlnPlot (after QC)", file=fp_html)
    print("</h3>", file=fp_html)
    
    # nCount_Spatial
    filename_qc_ncount_vln_png = NORMDIR_DAT + "/" + PREFIX_v + "-QC-VLN-nCount_Spatial.png"
    filename_qc_ncount_vln_png_rel = "dat" + "/" + PREFIX_v + "-QC-VLN-nCount_Spatial.png"
    print("<img src=" + "\"" + filename_qc_ncount_vln_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # nFeature_Spatial
    filename_qc_nfeature_vln_png = NORMDIR_DAT + "/" + PREFIX_v + "-QC-VLN-nFeature_Spatial.png"
    filename_qc_nfeature_vln_png_rel = "dat" + "/" + PREFIX_v + "-QC-VLN-nFeature_Spatial.png"    
    print("<img src=" + "\"" + filename_qc_nfeature_vln_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # percent_mito
    filename_qc_percent_mito_vln_png = NORMDIR_DAT + "/" + PREFIX_v + "-QC-VLN-percent_mito.png"
    filename_qc_percent_mito_vln_png_rel = "dat" + "/" + PREFIX_v + "-QC-VLN-percent_mito.png"    
    print("<img src=" + "\"" + filename_qc_percent_mito_vln_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # percent_hb
    filename_qc_percent_hb_vln_png = NORMDIR_DAT + "/" + PREFIX_v + "-QC-VLN-percent_hb.png"
    filename_qc_percent_hb_vln_png_rel = "dat" + "/" + PREFIX_v + "-QC-VLN-percent_hb.png"
    print("<img src=" + "\"" + filename_qc_percent_hb_vln_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
    
    # percent_ribo
    filename_qc_percent_ribo_vln_png = NORMDIR_DAT + "/" + PREFIX_v + "-QC-VLN-percent_ribo.png"
    filename_qc_percent_ribo_vln_png_rel = "dat" + "/" + PREFIX_v + "-QC-VLN-percent_ribo.png"    
    print("<img src=" + "\"" + filename_qc_percent_ribo_vln_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)


    # FeatureScatter(2)
    print("<h3>", file=fp_html)
    print("FeatureScatter (after QC)", file=fp_html)
    print("</h3>", file=fp_html)
    
    # nCount_Spatial vs nFeature_Spatial
    filename_qc_fss_png = NORMDIR_DAT + "/" + PREFIX_v + "-QC-FSS-nCount_Spatial-nFeature_Spatial.png"
    filename_qc_fss_png_rel = "dat" + "/" + PREFIX_v + "-QC-FSS-nCount_Spatial-nFeature_Spatial.png"
    print("<img src=" + "\"" + filename_qc_fss_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)


    # SpatialFeaturePlot(2)
    print("<h3>", file=fp_html)
    print("SpatialFeaturePlot (after QC)", file=fp_html)
    print("</h3>", file=fp_html)
    
    # nCount_Spatial
    filename_qc_ncount_sft_png = NORMDIR_DAT + "/" + PREFIX_v + "-QC-SFT-nCount_Spatial.png"
    filename_qc_ncount_sft_png_rel = "dat" + "/" + PREFIX_v + "-QC-SFT-nCount_Spatial.png"    
    print("<img src=" + "\"" + filename_qc_ncount_sft_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # nFeature_Spatial
    filename_qc_nfeature_sft_png = NORMDIR_DAT + "/" + PREFIX_v + "-QC-SFT-nFeature_Spatial.png"
    filename_qc_nfeature_sft_png_rel = "dat" + "/" + PREFIX_v + "-QC-SFT-nFeature_Spatial.png"
    print("<img src=" + "\"" + filename_qc_nfeature_sft_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # percent_mito
    filename_qc_percent_mito_sft_png = NORMDIR_DAT + "/" + PREFIX_v + "-QC-SFT-percent_mito.png"
    filename_qc_percent_mito_sft_png_rel = "dat" + "/" + PREFIX_v + "-QC-SFT-percent_mito.png"
    print("<img src=" + "\"" + filename_qc_percent_mito_sft_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
   
    # percent_hb
    filename_qc_percent_hb_sft_png = NORMDIR_DAT + "/" + PREFIX_v + "-QC-SFT-percent_hb.png"
    filename_qc_percent_hb_sft_png_rel = "dat" + "/" + PREFIX_v + "-QC-SFT-percent_hb.png"
    print("<img src=" + "\"" + filename_qc_percent_hb_sft_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
    
    # percent_ribo
    filename_qc_percent_ribo_sft_png = NORMDIR_DAT + "/" + PREFIX_v + "-QC-SFT-percent_ribo.png"
    filename_qc_percent_ribo_sft_png_rel = "dat" + "/" + PREFIX_v + "-QC-SFT-percent_ribo.png"
    print("<img src=" + "\"" + filename_qc_percent_ribo_sft_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
    
    
    #####

    # After normalization
    print("<h2>", file=fp_html)
    print("Normalization", file=fp_html)
    print("</h2>", file=fp_html)

    # GroupCorrelation (for Seurat v4 only)
    DOIT = False 
    if DOIT:
        print("<h3>", file=fp_html)
        print("GroupCorrelation", file=fp_html)
        print("</h3>", file=fp_html)

        # log normalization
        filename_gcp_log_png = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-LOG-GCP.png"
        filename_gcp_log_png_rel = "dat" + "/" + PREFIX_v + "-NORM-LOG-GCP.png"    
        print("<img src=" + "\"" + filename_gcp_log_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

        # SCT normalization
        filename_gcp_sct_png = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-GCP.png"
        filename_gcp_sct_png_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-GCP.png"
        print("<img src=" + "\"" + filename_gcp_sct_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)
    

    # After log normalization
    print("<h3>", file=fp_html)
    print("Log normalization", file=fp_html)
    print("</h3>", file=fp_html)

    print("log1p(p*10000) is stored in Spatial@data<br>", file=fp_html)    

    filename_slice_txt = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-LOG-slice_before_GroupCorrelation.txt"
    filename_slice_txt_rel = "dat" + "/" + PREFIX_v + "-NORM-LOG-slice_before_GroupCorrelation.txt"
    print("<a href=\"" + filename_slice_txt_rel + "\">" + "str() of Seurat object (before GroupCorrelation)" + "</a><br>", file=fp_html)

    filename_slice_txt = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-LOG-slice.txt"
    filename_slice_txt_rel = "dat" + "/" + PREFIX_v + "-NORM-LOG-slice.txt"
    print("<a href=\"" + filename_slice_txt_rel + "\">" + "str() of Seurat object (after GroupCorrelation)" + "</a><br>", file=fp_html)

    #####
    
    # After SCTransform normalization
    print("<h3>", file=fp_html)
    print("SCTransform normalization", file=fp_html)
    print("</h3>", file=fp_html)

    print("residuals (normalized values) is stored in SCT@scale.data<br>", file=fp_html)    
    print("log1p(corrected counts) is stored in SCT@data<br>", file=fp_html)
    
    filename_slice_txt = NORMDIR_DAT + "/" +  PREFIX_v + "-NORM-SCT-slice_before_GroupCorrelation.txt"
    filename_slice_txt_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-slice_before_GroupCorrelation.txt"
    print("<a href=\"" + filename_slice_txt_rel + "\">" + "str() of Seurat object (before GroupCorrelation)" + "</a><br>", file=fp_html)

    filename_slice_txt = NORMDIR_DAT + "/" +  PREFIX_v + "-NORM-SCT-slice.txt"
    filename_slice_txt_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-slice.txt"
    print("<a href=\"" + filename_slice_txt_rel + "\">" + "str() of Seurat object (after GroupCorrelation)" + "</a><br>", file=fp_html)


    # VlnPlot(3)
    print("<h3>", file=fp_html)
    print("VlnPlot (after SCT normalization)", file=fp_html)
    print("</h3>", file=fp_html)
    
    # nCount_Spatial
    filename_sct_ncount_vln_png = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-VLN-nCount_SCT.png"
    filename_sct_ncount_vln_png_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-VLN-nCount_SCT.png"
    print("<img src=" + "\"" + filename_sct_ncount_vln_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # nFeature_SCT
    filename_sct_nfeature_vln_png = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-VLN-nFeature_SCT.png"
    filename_sct_nfeature_vln_png_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-VLN-nFeature_SCT.png"
    print("<img src=" + "\"" + filename_sct_nfeature_vln_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # percent_mito
    filename_sct_percent_mito_vln_png = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-VLN-percent_mito.png"
    filename_sct_percent_mito_vln_png_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-VLN-percent_mito.png"
    print("<img src=" + "\"" + filename_sct_percent_mito_vln_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # percent_hb
    filename_sct_percent_hb_vln_png = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-VLN-percent_hb.png"
    filename_sct_percent_hb_vln_png_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-VLN-percent_hb.png"
    print("<img src=" + "\"" + filename_sct_percent_hb_vln_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # percent_ribo
    filename_sct_percent_ribo_vln_png = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-VLN-percent_ribo.png"
    filename_sct_percent_ribo_vln_png_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-VLN-percent_ribo.png"    
    print("<img src=" + "\"" + filename_sct_percent_ribo_vln_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    #####


    print("<br>", file=fp_html)

    # nCount_SCT, txt
    filename_ncount_sct_txt = NORMDIR_DAT + "/" + PREFIX_v + "-nCount_SCT.txt"
    filename_ncount_sct_txt_rel = "dat" + "/" + PREFIX_v + "-nCount_SCT.txt"
    print("<a href=\"" + filename_ncount_sct_txt_rel + "\">" + "nCount_SCT" + "</a>", file=fp_html)

    print("<br>", file=fp_html)    

    # nFeature_SCT, txt
    filename_nfeature_sct_txt = NORMDIR_DAT + "/" + PREFIX_v + "-nFeature_SCT.txt"
    filename_nfeature_sct_txt_rel = "dat" + "/" + PREFIX_v + "-nFeature_SCT.txt"
    print("<a href=\"" + filename_nfeature_sct_txt_rel + "\">" + "nFeature_SCT" + "</a>", file=fp_html)
    

    # FeatureScatter(3)
    print("<h3>", file=fp_html)
    print("FeatureScatter (after SCT normalization)", file=fp_html)
    print("</h3>", file=fp_html)

    # nCount_SCT vs nFeature_SCT
    filename_sct_fss_sct_sct_png = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-FSS-nCount_SCT-nFeature_SCT.png"
    filename_sct_fss_sct_sct_png_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-FSS-nCount_SCT-nFeature_SCT.png"    
    print("<img src=" + "\"" + filename_sct_fss_sct_sct_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # nCount_Spatial vs nFeature_SCT
    filename_sct_fss_spatial_sct_png = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-FSS-nCount_Spatial-nFeature_SCT.png"
    filename_sct_fss_spatial_sct_png_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-FSS-nCount_Spatial-nFeature_SCT.png"    
    print("<img src=" + "\"" + filename_sct_fss_spatial_sct_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # nCount_Spatial vs nCount_SCT
    filename_sct_fss_spatial_sct_png = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-FSS-nCount_Spatial-nCount_SCT.png"
    filename_sct_fss_spatial_sct_png_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-FSS-nCount_Spatial-nCount_SCT.png"    
    print("<img src=" + "\"" + filename_sct_fss_spatial_sct_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    
    # SpatialFeaturePlot(3)
    print("<h3>", file=fp_html)
    print("SpatialFeaturePlot (after SCT normalization)", file=fp_html)
    print("</h3>", file=fp_html)

    # nCount_SCT
    filename_sct_ncount_sft_png = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-SFT-nCount_SCT.png"
    filename_sct_ncount_sft_png_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-SFT-nCount_SCT.png"    
    print("<img src=" + "\"" + filename_sct_ncount_sft_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # nFeature_SCT
    filename_sct_nfeature_sft_png = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-SFT-nFeature_SCT.png"
    filename_sct_nfeature_sft_png_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-SFT-nFeature_SCT.png"
    print("<img src=" + "\"" + filename_sct_nfeature_sft_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # percent_mito
    filename_sct_percent_mito_sft_png = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-SFT-percent_mito.png"
    filename_sct_percent_mito_sft_png_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-SFT-percent_mito.png"
    print("<img src=" + "\"" + filename_sct_percent_mito_sft_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # percent_hb
    filename_sct_percent_hb_sft_png = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-SFT-percent_hb.png"
    filename_sct_percent_hb_sft_png_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-SFT-percent_hb.png"
    print("<img src=" + "\"" + filename_sct_percent_hb_sft_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

    # percent_ribo
    filename_sct_percent_ribo_sft_png = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-SFT-percent_ribo.png"
    filename_sct_percent_ribo_sft_png_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-SFT-percent_ribo.png"
    print("<img src=" + "\"" + filename_sct_percent_ribo_sft_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)


    
    print("<h2>", file=fp_html)
    print("Datasets", file=fp_html)
    print("</h2>", file=fp_html)


    ##### #####
    # CSV files

    print("<h3>", file=fp_html)
    print("CSV files", file=fp_html)
    print("</h3>", file=fp_html)


    ##### SCT

    # SCT@counts
    filename_sct_cc = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-counts.txt"
    filename_sct_cc_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-counts.txt"    
    print("<a href=\"" + filename_sct_cc_rel + "\">" + "SCT@counts" + "</a><br>", file=fp_html)

    # SCT@data
    filename_sct_dd = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-data.txt"
    filename_sct_dd_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-data.txt"    
    print("<a href=\"" + filename_sct_dd_rel + "\">" + "SCT@data" + "</a><br>", file=fp_html)

    # SCT@scale.data
    filename_sct_sd = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-SCT-scale_data.txt"
    filename_sct_sd_rel = "dat" + "/" + PREFIX_v + "-NORM-SCT-scale_data.txt"    
    print("<a href=\"" + filename_sct_sd_rel + "\">" + "SCT@scale.data" + "</a><br>", file=fp_html)
    

    ##### Spatial

    # Spatial@counts
    filename_spatial_cc = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-Spatial-counts.txt"
    filename_spatial_cc_rel = "dat" + "/" + PREFIX_v + "-NORM-Spatial-counts.txt"    
    print("<a href=\"" + filename_spatial_cc_rel + "\">" + "Spatial@counts" + "</a><br>", file=fp_html)

    # Spatial@dataq
    filename_spatial_dd = NORMDIR_DAT + "/" + PREFIX_v + "-NORM-Spatial-data.txt"
    filename_spatial_dd_rel = "dat" + "/" + PREFIX_v + "-NORM-Spatial-data.txt"    
    print("<a href=\"" + filename_spatial_dd_rel + "\">" + "Spatial@data" + "</a><br>", file=fp_html)


    print("<h2>", file=fp_html)
    print("Examples", file=fp_html)
    print("</h2>", file=fp_html)

    #gene_name = "MS4A1"

    def MakeSFP(gene_name):

        print("P:" + gene_name)

        print("<h3>", file=fp_html)
        print(gene_name, file=fp_html)
        print("</h3>", file=fp_html)

        print("SCT@counts, SCT@data, SCT@scale.data, Spatial@counts, and Spatial@data<br>", file=fp_html)

        # SCT
        data_ss = "SCT_counts"
        filename_gene_png = NORMDIR_DAT + "/" + PREFIX_v + "-SFP-" + data_ss + "-" + gene_name + ".png"
        filename_gene_png_rel = "dat" + "/" + PREFIX_v + "-SFP-" + data_ss + "-" + gene_name + ".png"
        print("<img src=" + "\"" + filename_gene_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

        data_ss = "SCT_data"
        filename_gene_png = NORMDIR_DAT + "/" + PREFIX_v + "-SFP-" + data_ss + "-" + gene_name + ".png"
        filename_gene_png_rel = "dat" + "/" + PREFIX_v + "-SFP-" + data_ss + "-" + gene_name + ".png"
        print("<img src=" + "\"" + filename_gene_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

        data_ss = "SCT_scale_data"
        filename_gene_png = NORMDIR_DAT + "/" + PREFIX_v + "-SFP-" + data_ss + "-" + gene_name + ".png"
        filename_gene_png_rel = "dat" + "/" + PREFIX_v + "-SFP-" + data_ss + "-" + gene_name + ".png"
        print("<img src=" + "\"" + filename_gene_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

        # Spatial
        data_ss = "Spatial_counts"
        filename_gene_png = NORMDIR_DAT + "/" + PREFIX_v + "-SFP-" + data_ss + "-" + gene_name + ".png"
        filename_gene_png_rel = "dat" + "/" + PREFIX_v + "-SFP-" + data_ss + "-" + gene_name + ".png"
        print("<img src=" + "\"" + filename_gene_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

        data_ss = "Spatial_data"
        filename_gene_png = NORMDIR_DAT + "/" + PREFIX_v + "-SFP-" + data_ss + "-" + gene_name + ".png"
        filename_gene_png_rel = "dat" + "/" + PREFIX_v + "-SFP-" + data_ss + "-" + gene_name + ".png"
        print("<img src=" + "\"" + filename_gene_png_rel + "\"" + " width=\"" + ww + "\"" + " border=\"0\"" + ">", file=fp_html)

        
    MakeSFP("MS4A1")
    MakeSFP("CD68")
    MakeSFP("CDH1")
    MakeSFP("TTF1")            
    MakeSFP("CD4") 
    MakeSFP("PTPRC") # CD45
    MakeSFP("APOE")
    
    #####
    
    print_html_footer(fp_html)        
        
    fp_html.close()
    

#####

if __name__ == "__main__":

    main()
    
