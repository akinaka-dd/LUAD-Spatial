#!/usr/bin/env python3

import numpy as np
import cv2
import math


purple = [153,89,119]
pink = [132,70,202]
lightred = [132,70,202]
orangered = [0, 69, 255]
white = [255,255,255]
black = [0,0,0]
blue = [255,0,0]    
green = [0,255,0]
magenta = [255,0,255]
darkmagenta = [128,0,128]    
lightblue =  [202,70,70]
darkgreen = [0,160, 0]
royalblue = [225, 105, 65]
lightgray = [210,210,210]

alpha = 255

binary_th = 10

def bg_transp(img):

    ss_num = len(img.shape)    
    if ss_num == 2: # grayscale
        img_bgra = cv2.cvtColor(img, cv2.COLOR_GRAY2BGRA)
        img_bgr = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)            
    elif ss_num == 3: # BGR or BGRA

        h, w, c = img.shape

        if c == 3:
            img_bgra = cv2.cvtColor(img, cv2.COLOR_BGR2BGRA)
            img_bgr = img            
        elif c == 4:
            img_bgra = img
            img_bgr = cv2.cvtColor(img, cv2.COLOR_BGRA2BGR)
            
    else:
        print(ss_num, "BG_TRANSP OTHERWISE")
        sys.exit()

    img_bgra[:, :, 3] = np.where(np.all(img_bgr == 255, axis=-1), 0, 255)

    return img_bgra

def COUNT(img_visium_mask_with_hole, img_codex_ww_opt):

    hh,ww = img_codex_ww_opt.shape

    img_with_hole_gray = np.zeros_like(img_visium_mask_with_hole)
    hh0,ww0 = img_with_hole_gray.shape

    hhx = hh
    wwx = ww
    if hh0 < hh:
        flag_clipped = True
        hhx = hh0
    if ww0 < ww:
        flag_clipped = True
        wwx = ww0

    bin_th_vis, img_vis_bin = cv2.threshold(img_visium_mask_with_hole, binary_th, 255, cv2.THRESH_BINARY)
    bin_th_cod, img_cod_bin = cv2.threshold(img_codex_ww_opt, binary_th, 255, cv2.THRESH_BINARY)

    img_with_hole_and = cv2.bitwise_and(img_vis_bin[0:hhx,0:wwx], img_cod_bin[0:hhx,0:wwx])
    img_with_hole_or = cv2.bitwise_or(img_vis_bin[0:hhx,0:wwx], img_cod_bin[0:hhx,0:wwx])

    count_and = cv2.countNonZero(img_with_hole_and)
    count_or = cv2.countNonZero(img_with_hole_or)

    count_vis = cv2.countNonZero(img_vis_bin[0:hhx,0:wwx])
    count_cod = cv2.countNonZero(img_cod_bin[0:hhx,0:wwx])

    return count_and, count_or, count_vis, count_cod


## ZOOM

def ZOOM(img_mask, img_mask_with_hole, h_max, w_max, ratio):

    img_mask_resized = cv2.resize(img_mask, dsize=None, fx=ratio, fy=ratio, interpolation=cv2.INTER_NEAREST)
    img_mask_with_hole_resized = cv2.resize(img_mask_with_hole, dsize=None, fx=ratio, fy=ratio, interpolation=cv2.INTER_NEAREST)
    
    h_p, w_p = img_mask.shape

    img_mask_resized_r = np.zeros((h_max, w_max), dtype=np.uint8)    
    h_pr, w_pr = img_mask_resized.shape

    h_vv = h_pr
    if h_vv > h_max:
        h_vv = h_max
    
    w_vv = w_pr
    if w_vv > w_max:
        w_vv = w_max
    
    img_mask_resized_r[:h_vv, :w_vv] = img_mask_resized[:h_vv, :w_vv]

    img_mask_with_hole_resized_r = np.zeros((h_max, w_max), dtype=np.uint8)        
    h_pr, w_pr = img_mask_with_hole_resized.shape    
    img_mask_with_hole_resized_r[:h_vv, :w_vv] = img_mask_with_hole_resized[:h_vv, :w_vv]

    return img_mask_resized_r, img_mask_with_hole_resized_r


## SHIFT

def SHIFT(img_mask_resized, img_mask_with_hole_resized, dx, dy):

    rows,cols = img_mask_resized.shape
    M = np.float32([[1,0,dx],[0,1,dy]])
    img_mask_resized_shifted = cv2.warpAffine(img_mask_resized, M, (cols,rows))
    img_mask_with_hole_resized_shifted = cv2.warpAffine(img_mask_with_hole_resized, M, (cols,rows))

    return img_mask_resized_shifted, img_mask_with_hole_resized_shifted


def CENTERING(img_mask, img_mask_with_hole, xc, yc):

    contours, hierachy = cv2.findContours(img_mask, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    contour_max_area = max(contours, key=lambda x: cv2.contourArea(x))

    mom = cv2.moments(contour_max_area)
    xc_r = int(mom["m10"]/mom["m00"])
    yc_r = int(mom["m01"]/mom["m00"])

    dx = xc - xc_r
    dy = yc - yc_r

    img_cc, img_wh_cc = \
        SHIFT(img_mask, img_mask_with_hole, dx, dy)

    return img_cc, img_wh_cc, dx, dy


## ROTATE

def ROTATE(img_mask_resized_shifted, img_mask_with_hole_resized_shifted, xc, yc, h_max, w_max, degree):

    if degree % 180 == 0:

        #h_v, w_v = img_mask_resized_shifted.shape
        sp = img_mask_resized_shifted.shape
        h_v = sp[0]
        w_v = sp[1]
        
        h_x = min(h_v, h_max)
        w_x = min(w_v, w_max)
        
        img_mask_resized_shifted_rotated = np.zeros((h_max, w_max), dtype=np.uint8)
        img_mask_resized_shifted_rotated[:h_x, :w_x] = img_mask_resized_shifted[:h_x, :w_x]

        img_mask_with_hole_resized_shifted_rotated = np.zeros((h_max, w_max), dtype=np.uint8)
        img_mask_with_hole_resized_shifted_rotated[:h_x, :w_x] = img_mask_with_hole_resized_shifted[:h_x, :w_x]

        return img_mask_resized_shifted_rotated, img_mask_with_hole_resized_shifted_rotated

    R = cv2.getRotationMatrix2D((xc,yc), degree, 1.0)
    img_mask_resized_shifted_rotated = cv2.warpAffine(img_mask_resized_shifted, R, (w_max, h_max), flags=cv2.INTER_NEAREST)
    img_mask_with_hole_resized_shifted_rotated = cv2.warpAffine(img_mask_with_hole_resized_shifted, R, (w_max, h_max), flags=cv2.INTER_NEAREST)

    return img_mask_resized_shifted_rotated, img_mask_with_hole_resized_shifted_rotated


def rotate_img(img, dg):

    h, w, _ = img.shape

    dg_rad = np.pi * (dg / 180)

    w_rotated = int(h*np.absolute(np.sin(dg_rad)) + w*np.absolute(np.cos(dg_rad)))
    h_rotated = int(h*np.absolute(np.cos(dg_rad)) + w*np.absolute(np.sin(dg_rad)))

    M = cv2.getRotationMatrix2D((w/2,h/2), dg, 1.0)

    M[0][2] = M[0][2] - w/2 + w_rotated/2
    M[1][2] = M[1][2] - h/2 + h_rotated/2

    img_rotated = cv2.warpAffine(img, M, (w_rotated,h_rotated), flags=cv2.INTER_NEAREST)

    return img_rotated


