# LUAD-Spatial
**Multimodal spatial omics analysis of human LUAD by Visium and PhenoCycler**

## DOIT_LUAD3B.py
- A Python script integrates spatial omics data by 10x Genomics Visium and Akoya Biosciences PhenoCycler.
  - Supposes Visium (not CytAssist) and PhenoCycler-Open (not Fusion) in the example of this repository.
- Imvokes the following DOIT-xxxx scripts (1 to 8).
  - Tested with Python 3.12.
  - Requires packages such as numpy, pandas, cv2, tifffile, depending on environment.
  - Internally invokes R scripts.
	- Tested with R 4.3.2 and Seurat v4.4.0 (and will be compatible with v5).
- Results are stored in "FFPE_LUAD_3_B" directory in the current working directory.

``` 4d
% path_to_the_src_directory/src/DOIT_LUAD3B.py
```

### 1. DOIT_GetDATA_LUAD3B.py
- Downloads and expands the datafile (LUAD3B_dat.tgz) in "dat" directory in the current working directory. 

### 2. DOIT_Visium_LUAD3B.py
- Generates mask images of Visium tissue region using H&E stain images in the output of SpaceRanger ("dat/Visium_FFPE_LUAD_3_B/spatial").
- Results are browsable from "FFPE_LUAD_3_B/FFPE_LUAD_3_B-Visium-MASK/FFPE_LUAD_3_B-Visium.html".

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
- Generate mask images of PhenoCycler tissue region using fluorescence intensity in the QPTIFF data by CODEX Processor ("dat/PhenoCycler_LUAD_3_B/Experiment_LUAD_3_B.qptiff")
- Results are browsable from "FFPE_LUAD_3_B/FFPE_LUAD_3_B-PhenoCycler-MASK/FFPE_LUAD_3_B-PhenoCycler.html".

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
- Normalizes Visium expression levels by SCTransform normalization (version 1) provided by Seurat.
  - Internally invokes an R script (exec_preprocess.r).
- Results are browsable from "FFPE_LUAD_3_B/FFPE_LUAD_3_B-Visium-NORM/FFPE_LUAD_3_B-Visium-NORM.html".

### 5. DOIT_Align_LUAD3B.py
- Aligns mask images of Visium (red) and PhenoCycler (blue) respectively obtained in step 2 and 3 above.
  - Mask images were scaled, translated, and rotated (affine transformations) to maximize their overlap.
  - Overlap is evaluated by their intersection over union (IoU), i.e., IoU = (Visium ∩ PhenoCycler) / (Visium ∪ PhenoCycler). 
  - By exhaustive grid search in the neighborhood of the solution obtained above, IoU was optimized to generate the final result. 
- Results are browsable from "FFPE_LUAD_3_B/FFPE_LUAD_3_B-ALIGN/FFPE_LUAD_3_B-Visium-PhenoCycler-ALIGN.html".

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
- Evalutes correlation between Visium gene expression level and PhenoCycler antibody fluorescence intensity.
  - In circular regions corresponding to Visium spots average of fluorescence intensities are calculated. 
  - Spatial and quantitative correlations are evaluated.
- Results are browsable from "FFPE_LUAD_3_B/FFPE_LUAD_3_B-EVAL/FFPE_LUAD_3_B-EVAL.html".
  
### Example of comparison
- MS4A1 gene expression by Visium and CD20 antibody fluorescence intensity by PhenoCycler.

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
	
### 7. DOIT_Clustering_LUAD3B.py
- Carries out clustering of Visium spots by Seurat.
  - Internally invokes an R script (exec_clustering.r).
  - Carries out dimensionality reduction and extracts spatially variable features.
- Results are browsable from "FFPE_LUAD_3_B/FFPE_LUAD_3_B-Visium-CLST/FFPE_LUAD_3_B-Visium-CLST.html".

### 8. DOIT_GetROI_LUAD3B.py
- Extracts an aligned/integrated region of interest (ROI) image.
- Results are browsable from "FFPE_LUAD_3_B/FFPE_LUAD_3_B-ROI/FFPE_LUAD_3_B-ROI.html".

#### Mutual coordinate mapping

<table>
<tr>
<td><img width="250px" src="img/roi/FFPE_LUAD_3_B-PhenoCycler-Visium-map_ext_01.png"></td>
<td><img width="250px" src="img/roi/FFPE_LUAD_3_B-PhenoCycler-Visium-map_ext_rotated_01.png"></td>
</tr>
<tr>
<td>PhenoCycler image coordinate (Visium on PhenoCycler)</td>
<td>Visium image coordinate (PhenoCycler on Visium)</td>
</tr>
</table>

#### Example of projection of a ROI in PhenoCycler image onto Visium image
#### 
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


