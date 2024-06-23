#!/usr/bin/env python3

import cv2
import ast

# LOAD masks
def load_masks(MASKDIR_v_DAT, PREFIX_v, MASKDIR_p_DAT_MULTI, PREFIX_p):

    # Visium
    filename_visium_mask = MASKDIR_v_DAT + "/" + PREFIX_v + "_tr_gray_without_circle_nega_mask.png"
    filename_visium_mask_wh = MASKDIR_v_DAT  + "/" + PREFIX_v + "_tr_gray_without_circle_nega_mask_with_hole.png"

    # PhenoCycler
    filename_phenocycler_mask_01 = MASKDIR_p_DAT_MULTI + "/" + PREFIX_p + "-" + "MULTI_gray_nega_mask_01.png"
    filename_phenocycler_mask_wh_01 = MASKDIR_p_DAT_MULTI + "/" + PREFIX_p + "-" + "MULTI_gray_nega_mask_with_hole_01.png"

    filename_phenocycler_mask = MASKDIR_p_DAT_MULTI + "/" + PREFIX_p + "-" + "MULTI_gray_nega_mask.png"
    filename_phenocycler_mask_wh = MASKDIR_p_DAT_MULTI + "/" + PREFIX_p + "-" + "MULTI_gray_nega_mask_with_hole.png"

    # with hole
    img_visium_mask_wh = cv2.imread(filename_visium_mask_wh, cv2.IMREAD_GRAYSCALE)
    img_phenocycler_mask_wh_01 = cv2.imread(filename_phenocycler_mask_wh_01, cv2.IMREAD_GRAYSCALE)
    img_phenocycler_mask_wh = cv2.imread(filename_phenocycler_mask_wh, cv2.IMREAD_GRAYSCALE)

    # without hole
    img_visium_mask = cv2.imread(filename_visium_mask, cv2.IMREAD_GRAYSCALE)
    img_phenocycler_mask_01 = cv2.imread(filename_phenocycler_mask_01, cv2.IMREAD_GRAYSCALE)
    img_phenocycler_mask = cv2.imread(filename_phenocycler_mask, cv2.IMREAD_GRAYSCALE)

    return \
        img_visium_mask_wh, img_phenocycler_mask_wh_01, img_phenocycler_mask_wh, \
        img_visium_mask, img_phenocycler_mask_01, img_phenocycler_mask

# LOAD PhenoCycler resolution
def load_phenocycler_resolution(filename_reso):

    fp_reso = open(filename_reso, "rt")

    line = fp_reso.readline()
    line = line.rstrip()
    vals = line.split("\t")
    resolution_nm = float(vals[1])
    resolution_um = resolution_nm / 1000

    print("RESOLUTION_NM", resolution_nm)
    print("RESOLUTION_UM", resolution_um)

    return resolution_nm, resolution_um
    
# LOAD PhenoCycler pages
def load_phenocycler_pages(PAGESDIR, PREFIX_p, filename_rgb):

    fp_rgb = open(filename_rgb, "rt")

    IMAGES = []
    PAGEID = {}
    PAGEID_pos = {}    

    count = 0
    line = fp_rgb.readline()
    for line in fp_rgb:
        line = line.rstrip()
        vals = line.split("\t")

        pid = vals[0]
        name = vals[1]
        name_org = vals[2]

        v = (pid, name, name_org)

        IMAGES.append(v)

        PAGEDIR = PAGESDIR + "/" + "PAGE{:s}".format(pid) + "-" + name
        filename_PAGE_gray_nega = PAGEDIR + "/" + PREFIX_p + "-" + "PAGE{:s}".format(pid) + "_gray_nega" + "-" + name + ".png"
        filename_PAGE_gray = PAGEDIR + "/" + PREFIX_p + "-" + "PAGE{:s}".format(pid) + "_gray" + "-" + name + ".png"        

        PAGEID[name] = (count, pid, filename_PAGE_gray_nega)
        PAGEID_pos[name] = (count, pid, filename_PAGE_gray)

        count += 1

    return IMAGES, PAGEID, PAGEID_pos

# LOAD visium text files
def load_visium_text_files(DEGREE, MASKDIR_v_DAT):

    filename_json = MASKDIR_v_DAT + "/" + "scalefactors_json.json"
    filename_pos = MASKDIR_v_DAT + "/" + "tissue_positions_list.csv"
    filename_ft = MASKDIR_v_DAT + "/" + "features.tsv"
    filename_fmp = MASKDIR_v_DAT + "/" + "fiducial_positions_list.txt"

    fp_json = open(filename_json, "rt")
    line = fp_json.readline()
    SCALEFACTORS = ast.literal_eval(line)
    spot_diameter_fullres = SCALEFACTORS["spot_diameter_fullres"]
    tissue_hires_scalef = SCALEFACTORS["tissue_hires_scalef"]
    fiducial_diameter_fullres = SCALEFACTORS["fiducial_diameter_fullres"]
    scale = tissue_hires_scalef
    fp_json.close()

    fp_pos = open(filename_pos, "rt")

    SPOT = {}
    GRID_rc = {}
    GRID_cr = {}

    for line in fp_pos:

        line = line.rstrip()
        vals = line.split(",")

        bc = vals[0]
        in_tissue = int(vals[1])
        ri = int(vals[2])
        ci = int(vals[3])
        yv_f = int(vals[4])
        xv_f = int(vals[5])

        xx = xv_f
        yy = yv_f

        if DEGREE == cv2.ROTATE_90_COUNTERCLOCKWISE: # direction of the y-axis
            xx = yv_f # x*cos(-90) - y*sin(-90) = 
            yy = int(hv_h / scale) - xv_f # H + (x*sin(-90) + y*cos(-90))

        SPOT[bc] = (xx, yy, in_tissue, ci, ri)

        if not ri in GRID_rc:
            GRID_rc[ri] = {}
        GRID_rc[ri][ci] = (xx, yy)

        if not ci in GRID_cr:
            GRID_cr[ci] = {}
        GRID_cr[ci][ri] = (xx, yy)

    FID = {}
    fp_ft = open(filename_ft, "rt")
    for line in fp_ft:
        line = line.rstrip()
        vals = line.split()

        fid = vals[0]
        gene = vals[1]

        FID[gene] = fid

    FMP = [] # fiducial marker position
    fp_fmp = open(filename_fmp, "rt")
    for line in fp_fmp:
        line = line.rstrip()
        vals = line.split("\t")

        xc = int(vals[0])
        yc = int(vals[1])

        FMP.append((xc, yc))

    return SCALEFACTORS, spot_diameter_fullres, tissue_hires_scalef, fiducial_diameter_fullres, SPOT, GRID_rc, GRID_cr, FID, FMP


# LOAD visium cluster text file
def load_visium_clusters(filename_seurat_clusters_txt):

    fp_cls = open(filename_seurat_clusters_txt, "rt")

    CLUSTER = {}

    count = 0
    line = fp_cls.readline()
    for line in fp_cls:
        line = line.rstrip()
        vals = line.split(",")

        bc = vals[0]
        cluster_id = int(vals[1])

        #print(bc, cluster_id)
        CLUSTER[bc] = cluster_id

        count += 1

    return CLUSTER

