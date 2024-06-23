# LUAD-Spatial
Spatial omics analysis of human LUAD by Visium and PhenoCycler

## DOIT_LUAD3B.py
- Executes the following DOIT-xxxx scripts (1 to 7).
- Results are stored in FFPE_LUAD_3_B directory in the current directory.

``` 4d
% path_to_the_src_directory/src/DOIT_LUAD3B.py
```

### 1. DOIT_GetDATA_LUAD3B.py

### 2. DOIT_Visium_LUAD3B.py

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


### 3. DOIT_PhenoCycler_LUAD3B.py

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


### 4. DOIT_Normalize_LUAD3B.py

Normalization.

### 5. DOIT_Align_LUAD3B.py

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

### 6. DOIT_Eval_LUAD3B.py

Visium (MS4A1) and PhenoCycler (CD20)</td>

<table>
<tr>
<td><img width="200px" src="img/eval/FFPE_LUAD_3_B-MS4A1-CD20_plot_bin_01.png"></td>
<td><![Alt Text](URL )  width="200px" src="img/eval/FFPE_LUAD_3_B-MS4A1-CD20_cor_in_tissue.png"></td>
</tr>
<tr>
<td>Spatial correlation</td>
<td>Quantitative correlation</td>
</tr>
</table>
	
	DOIT_Clustering_LUAD3B.py

### 7. DOIT_GetROI_LUAD3B.py

<table>
<tr>
<td><img width="250px" src="img/roi/FFPE_LUAD_3_B-PhenoCycler-Visium-map_ext_01.png"></td>
<td><img width="250px" src="img/roi/FFPE_LUAD_3_B-PhenoCycler-Visium-map_ext_rotated_01.png"></td>
</tr>
<tr>
<td>Visium image aligned to PhenoCycler image</td>
<td>PhenoCycler image aligned to Visium image</td>
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

#### Example of ROI (1)

<table>
<tr>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05-celltype_spots.png">
<td><img height="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05-celltype_legend.png">
<td><img width="300px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05-celltype_barchart.png">
</tr>
</table>


<table>
<tr>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05_heatmap_spots.png">
<td><img height="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05-heatmap_legend.png">
</tr>
</table>

<table>
<tr>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05_PAGE01-DAPI_rgba_spots.png">
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05_PAGE05-E-Cadherin_rgba_spots.png">
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05_PAGE03-CD4_rgba_spots.png">
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05_PAGE17-CD8_rgba_spots.png">
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05_PAGE10-CD20_rgba_spots.png">
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05_PAGE09-CD68_rgba_spots.png">
<td><img width="50px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05_colorbar_rgba.png"> 
</tr>
<tr>
<td>DAPI</td>
<td>E-cadherin</td>
<td>CD4</td>
<td>CD8</td>
<td>CD20</td>
<td>CD68</td>
</tr>
</table>

<table>
<tr>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05_visium.png"></td>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05_visium_spots-CDH1.png"></td>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05_visium_spots-CD4.png"></td>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05_visium_spots-CD8A.png"></td>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05_visium_spots-MS4A1.png"></td>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05_visium_spots-CD68.png"></td>
<td><img width="50px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI05_visium_spots_colorbar.png"></td>
</tr>
<tr>
<td>H&E stain</td>
<td>CDH1</td>
<td>CD4</td>
<td>CD8A</td>
<td>MS4A1</td>
<td>CD68</td>
</tr>
</table>


#### Example of ROI (2)

<table>
<tr>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00-celltype_spots.png">
<td><img height="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00-celltype_legend.png">
<td><img width="300px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00-celltype_barchart.png">
</tr>
</table>


<table>
<tr>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00_heatmap_spots.png">
<td><img height="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00-heatmap_legend.png">
</tr>
</table>

<table>
<tr>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00_PAGE01-DAPI_rgba_spots.png">
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00_PAGE05-E-Cadherin_rgba_spots.png">
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00_PAGE03-CD4_rgba_spots.png">
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00_PAGE17-CD8_rgba_spots.png">
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00_PAGE10-CD20_rgba_spots.png">
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00_PAGE09-CD68_rgba_spots.png">
<td><img width="50px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00_colorbar_rgba.png"> 
</tr>
<tr>
<td>DAPI</td>
<td>E-cadherin</td>
<td>CD4</td>
<td>CD8</td>
<td>CD20</td>
<td>CD68</td>
</tr>
</table>


<table>
<tr>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00_visium.png"></td>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00_visium_spots-CDH1.png"></td>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00_visium_spots-CD4.png"></td>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00_visium_spots-CD8A.png"></td>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00_visium_spots-MS4A1.png"></td>
<td><img width="150px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00_visium_spots-CD68.png"></td>
<td><img width="50px" src="img/roi/FFPE_LUAD_3_B-Visium-PhenoCycler_MULTI-ROI00_visium_spots_colorbar.png"></td>
</tr>
<tr>
<td>H&E stain</td>
<td>CDH1</td>
<td>CD4</td>
<td>CD8A</td>
<td>MS4A1</td>
<td>CD68</td>
</tr>
</table>


