# LUAD-Spatial
Spatial omics analysis of human LUAD by Visium and PhenoCycler

DOIT_LUAD3B.py

DOIT_GetDATA_LUAD3B.py

DOIT_Visium_LUAD3B.py

<table>
<tr>
<td><img src="img/visium/FFPE_LUAD_3_B-Visium_00.png"></td>
<td><img src="img/visium/FFPE_LUAD_3_B-Visium_tr_gray_with_circle_bgr.png"></td>
<td><img src="img/visium/FFPE_LUAD_3_B-Visium_tr_gray_without_circle_nega_bin.png"></td>
<td><img src="img/visium/FFPE_LUAD_3_B-Visium_tr_gray_without_circle_nega_bin_bgr_cont.png"></td>
<td><img src="img/visium/FFPE_LUAD_3_B-Visium_tr_gray_without_circle_nega_mask_with_hole_bgr.png"></td>
</tr>
<tr>
<td>H&E stain image</td>
<td>Fiducial marker removal</td>
<td>Binarization</td>
<td>Contour extraction</td>
<td>Mask image</td>
</tr>
</table>


DOIT_PhenoCycler_LUAD3B.py

<table>
<tr>
<td><img src="img/phenocycler/FFPE_LUAD_3_B-PhenoCycler-MULTI_bgr_bb_01.png"></td>
<td><img src="img/phenocycler/FFPE_LUAD_3_B-PhenoCycler-MULTI_gray_nega_01.png"></td>
<td><img src="img/phenocycler/FFPE_LUAD_3_B-PhenoCycler-MULTI_gray_nega_bin_01.png"></td>
<td><img src="img/phenocycler/FFPE_LUAD_3_B-PhenoCycler-MULTI_gray_nega_bin_dilate_erode_dilate_cont_01.png"></td>
<td><img src="img/phenocycler/FFPE_LUAD_3_B-PhenoCycler-MULTI_gray_nega_mask_with_hole_bgr_01.png"></td>
</tr>
<tr>
<td>Fluorescence image</td>
<td>Grayscaling</td>
<td>Binarization</td>
<td>Contour extraction</td>
<td>Mask image</td>
</tr>
</table>


DOIT_Normalize_LUAD3B.py

DOIT_Align_LUAD3B.py

<table>
<tr>
<td><img src="img/align/FFPE_LUAD_3_B_init_wh.png" ></t>
<td><img src="img/align/FFPE_LUAD_3_B_zoom_wh.png"></td>
<td><img src="img/align/FFPE_LUAD_3_B_shift_wh.png"></td>
<td><img src="img/align/FFPE_LUAD_3_B_rotate_wh.png"></td>
<td><img src="img/align/FFPE_LUAD_3_B_opt_wh.png"></td>
</tr>
<tr>
<td></td>
<td>Scaling</td>
<td>Translation</td>
<td>Rotation</td>
<td>Optimization</td>
</tr>
</table>

<video controls width="150" muted="false" src="https://github.com/akinaka-dd/LUAD-Spatial/assets/173239467/52cc85bb-f3b8-482d-98ca-164f7bfcb809"></video>
</div>

DOIT_Eval_LUAD3B.py

Visium (MS4A1) and PhenoCycler (CD20)</td>

<table>
<tr>
<td><img width="200px" src="img/eval/FFPE_LUAD_3_B-MS4A1-CD20_plot_bin_01.png"></td>
<td><img width="200px" src="img/eval/FFPE_LUAD_3_B-MS4A1-CD20_cor_in_tissue.png"></td>
</tr>
<tr>
<td>Spatial correlation</td>
<td>Quantitative correlation</td>
</tr>
</table>

DOIT_Clustering_LUAD3B.py

DOIT_GetROI_LUAD3B.py

<table>
<tr>
<td><img width="200px" src="img/roi/FFPE_LUAD_3_B-PhenoCycler-Visium-map_ext_01.png"></td>
<td><img width="200px" src="img/roi/FFPE_LUAD_3_B-PhenoCycler-Visium-map_ext_rotated_01.png"></td>
</tr>
<tr>
<td>Visium on PhenoCycler</td>
<td>PhenoCycler on Visium</td>
</tr>
</table>


<table>
<tr>
<td><img width="200px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-map0_01.png"></td>
<td><img width="200px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-map.png"></td>
</tr>
<tr>
<td>A region in the PhenoCycler coordinate</td>
<td>Projection onto the Visium coordinate</td>
</tr>
</table>



