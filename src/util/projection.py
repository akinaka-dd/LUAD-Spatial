#!/usr/bin/env python3

import math

def load_matrix(filename_txt):

    fp_txt = open(filename_txt, "rt")

    rotate_total = 0    
    MATRIX = []
    for line in fp_txt:
        line = line.rstrip()
        vals = line.split("\t")

        opc = vals[0]
        ops = [] 

        if opc == "ZOOM":
            r = float(vals[1])
            ops = [r]

        elif opc == "SHIFT":
            dx = int(vals[1])
            dy = int(vals[2])
            ops = [dx, dy]

        elif opc == "ROTATE":
            dg = float(vals[1])
            xc = int(vals[2])
            yc = int(vals[3])
            ops = [dg, xc, yc]
            rotate_total += dg            

        MATRIX.append((opc, ops))

    return MATRIX, rotate_total

def f_v2p(M, scale, xv_f, yv_f): 

    # to hires resolution
    xv_h = xv_f * scale
    yv_h = yv_f * scale

    xt_h = xv_h
    yt_h = yv_h

    for transform in reversed(M):

        opc, ops = transform

        if opc == "ZOOM":
            ratio = ops[0]
            xt_h = xt_h / ratio
            yt_h = yt_h / ratio

        elif opc == "SHIFT":
            dx = ops[0]
            dy = ops[1]
            xt_h = xt_h - dx
            yt_h = yt_h - dy

        elif opc == "ROTATE":
            dg = ops[0] 
            xc = ops[1]
            yc = ops[2]
            xt_h = xt_h - xc
            yt_h = yt_h - yc
            sin = math.sin(math.radians(+1*dg))
            cos = math.cos(math.radians(+1*dg))
            xt_h2 = cos * xt_h - sin * yt_h
            yt_h2 = sin * xt_h + cos * yt_h
            xt_h = xt_h2 + xc
            yt_h = yt_h2 + yc

    xp_f = int(xt_h)
    yp_f = int(yt_h)

    return xp_f, yp_f


def f_p2v(M, scale, xp_f, yp_f): 

    xt_h = xp_f
    yt_h = yp_f

    for transform in M:

        opc, ops = transform

        if opc == "ZOOM":
            ratio = ops[0]
            xt_h = xt_h * ratio
            yt_h = yt_h * ratio

        elif opc == "SHIFT":
            dx = ops[0]
            dy = ops[1]
            xt_h = xt_h + dx
            yt_h = yt_h + dy

        elif opc == "ROTATE":
            dg = ops[0]
            xc = ops[1]
            yc = ops[2]
            xt_h = xt_h - xc
            yt_h = yt_h - yc
            sin = math.sin(math.radians(-1*dg))
            cos = math.cos(math.radians(-1*dg))
            xt_h2 = cos * xt_h - sin * yt_h
            yt_h2 = sin * xt_h + cos * yt_h
            xt_h = xt_h2 + xc
            yt_h = yt_h2 + yc

    xv_f = int(xt_h/scale)
    yv_f = int(yt_h/scale)

    return xv_f, yv_f

