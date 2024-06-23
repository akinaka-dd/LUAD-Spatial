#!/usr/bin/env Rscript

library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

library(stringr)

library(data.table)

library("parallel")
options(mc.cores = detectCores(logical = FALSE))


main = function() {

    args = commandArgs(trailingOnly=T)
    print(args)

    sid = args[1]
    SPACERANGERDIR = args[2]
    NORMDIR = args[3]

    PREFIX = paste0(sid, "-", "Visium")

    NORMDIR_DAT = paste0(NORMDIR, "/", "dat")

    print(paste("SID:", sid))    
    print(paste("SPACERANGERDIR:", SPACERANGERDIR))
    print(paste("NORMDIR:", NORMDIR))
    
    name = sid

    # LOAD
    dirname = SPACERANGERDIR
    slice = Load10X_Spatial(data.dir = dirname)

    print(slice)

    # Contents of the slice object    
    filename_slice_txt = paste0(NORMDIR_DAT, "/", PREFIX, "-000-slice.txt")
    print(filename_slice_txt)    
    sink(filename_slice_txt)
    print(str(slice))  # str(slice, max.level=2)
    sink()
    
    # percent_mito & percent_hb & percent_ribo
    print("percent_mito & percent_hb & percent_ribo")

    slice = PercentageFeatureSet(slice, "^MT-", col.name = "percent_mito") # human 
    slice = PercentageFeatureSet(slice, "^HB[^(P)]", col.name = "percent_hb") # human
    slice = PercentageFeatureSet(slice, "^RP[SL]", col.name = "percent_ribo") # human 

    # Contents of the slice object
    filename_slice_txt = paste0(NORMDIR_DAT, "/", PREFIX, "-PFS-slice.txt")
    print(filename_slice_txt)    
    sink(filename_slice_txt)
    print(str(slice)) 
    sink()

    # number of features and samples (genes and spots)
    slice_dim = dim(x=slice)
    print(slice_dim)

    # number of features
    slice_nrow = nrow(x = slice)
    print(slice_nrow)

    # number of samples
    slice_ncol = ncol(x = slice)
    print(slice_ncol)


    print(head(x=rownames(slice), n=5))
    print(tail(x=colnames(slice), n=5))
    

    #####
    # themes
    
    # theme for VlnPlot
    theme_v = theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.x=element_blank(), 
        axis.title.y=element_text(size=20),              
        plot.title=element_text(size=20),
        aspect.ratio=1              
        )
    
    # theme for SpatialFeaturePlot
    theme_f = theme(legend.position="right",
          plot.title=element_text(size=20, face="bold", hjust=0.5),
          plot.margin=margin(t=0,r=0,b=0,l=0),              
          legend.title=element_text(size=8, hjust=0.5),
          legend.text=element_text(size=8),
          aspect.ratio=1              
          )
    
    # theme for FeatureScatter
    theme_s = theme(#legend.position="right",
          plot.title=element_text(size=20, face="bold", hjust=0.5),
          plot.margin=margin(t=0,r=0,b=0,l=0),              
          legend.title=element_text(size=8, hjust=0.5),
          legend.text=element_text(size=8),
          aspect.ratio=1              
          )

    # VlnPlot(1)
    # nCount_Spatial
    filename_ncount_vln_png = paste0(NORMDIR_DAT, "/", PREFIX, "-VLN-nCount_Spatial.png")
    print(filename_ncount_vln_png)    
    plot = VlnPlot(slice, features = "nCount_Spatial", pt.size=0.5) + labs(title=name, y="nCount_Spatial") + NoLegend() + theme_v
    ggsave(filename_ncount_vln_png)

    # txt CSV
    #slice@meta.data$nCount_Spatial
    df_nc = select(slice@meta.data, nCount_Spatial)
    print(head(df_nc, 5))
    filename_ncount_txt = paste0(NORMDIR_DAT, "/", PREFIX, "-nCount_Spatial.txt")
    print(filename_ncount_txt)
    
    write.table(df_nc, file=filename_ncount_txt, sep=",", na="", row.names=T, col.names=NA, quote=F)
                
    # nFeature_Spatial
    filename_nfeature_vln_png = paste0(NORMDIR_DAT, "/", PREFIX, "-VLN-nFeature_Spatial.png")
    print(filename_nfeature_vln_png)    
    plot = VlnPlot(slice, features = "nFeature_Spatial", pt.size=0.5) + labs(title=name, y="nFeature_Spatial") + NoLegend() + theme_v
    ggsave(filename_nfeature_vln_png)

    # txt
    df_nf = select(slice@meta.data, nFeature_Spatial)
    print(head(df_nf, 5))
    filename_nfeature_txt = paste0(NORMDIR_DAT, "/", PREFIX, "-nFeature_Spatial.txt")
    print(filename_nfeature_txt)

    write.table(df_nf, file=filename_nfeature_txt, sep=",", na="", row.names=T, col.names=NA, quote=F)
    

    # percent_mito
    filename_percent_mito_vln_png = paste0(NORMDIR_DAT, "/", PREFIX, "-VLN-percent_mito.png")
    print(filename_percent_mito_vln_png)    
    plot = VlnPlot(slice, features = "percent_mito", pt.size=0.5) + labs(title=name, y="percent_mito") + NoLegend() + theme_v    
    ggsave(filename_percent_mito_vln_png)

    # percent_hb
    filename_percent_hb_vln_png = paste0(NORMDIR_DAT, "/", PREFIX, "-VLN-percent_hb.png")
    print(filename_percent_hb_vln_png)    
    plot = VlnPlot(slice, features = "percent_hb", pt.size=0.5) + labs(title=name, y="percent_hb") + NoLegend() + theme_v    
    ggsave(filename_percent_hb_vln_png)

    # percent_ribo
    filename_percent_ribo_vln_png = paste0(NORMDIR_DAT, "/", PREFIX, "-VLN-percent_ribo.png")
    print(filename_percent_ribo_vln_png)    
    plot = VlnPlot(slice, features = "percent_ribo", pt.size=0.5) + labs(title=name, y="percent_ribo") + NoLegend() + theme_v    
    ggsave(filename_percent_ribo_vln_png)

    # FeatureScatter(1)
    # nCount_Spatial vs nFeature_Spatial
    filename_fss_png = paste0(NORMDIR_DAT, "/", PREFIX, "-FSS-nCount_Spatial-nFeature_Spatial.png")
    print(filename_fss_png)    
    plot = FeatureScatter(slice, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial", pt.size=0.5) + labs(title=name) + NoLegend() + theme_s  
    ggsave(filename_fss_png)


    # SpatialFeaturePlot(1)
    # nCount_Spatial
    filename_ncount_sft_png = paste0(NORMDIR_DAT, "/", PREFIX, "-SFT-nCount_Spatial.png")
    print(filename_ncount_sft_png)    
    plot = SpatialFeaturePlot(slice, features = "nCount_Spatial", pt.size.factor=1.3, crop=FALSE) + labs(title=name) + theme_f   
    ggsave(filename_ncount_sft_png)

    # nFeature_Spatial
    filename_nfeature_sft_png = paste0(NORMDIR_DAT, "/", PREFIX, "-SFT-nFeature_Spatial.png")
    print(filename_nfeature_sft_png)    
    plot = SpatialFeaturePlot(slice, features = "nFeature_Spatial", pt.size.factor=1.3, crop=FALSE) + labs(title=name) + theme_f
    ggsave(filename_nfeature_sft_png)

    # percent_mito
    filename_percent_mito_sft_png = paste0(NORMDIR_DAT, "/", PREFIX, "-SFT-percent_mito.png")
    print(filename_percent_mito_sft_png)    
    plot = SpatialFeaturePlot(slice, features = "percent_mito", pt.size.factor=1.2, crop=FALSE) + labs(title=name) + theme_f
    ggsave(filename_percent_mito_sft_png)

    # percent_hb
    filename_percent_hb_sft_png = paste0(NORMDIR_DAT, "/", PREFIX, "-SFT-percent_hb.png")
    print(filename_percent_hb_sft_png)    
    plot = SpatialFeaturePlot(slice, features = "percent_hb", pt.size.factor=1.2, crop=FALSE) + labs(title=name) + theme_f    
    ggsave(filename_percent_hb_sft_png)

    # percent_ribo
    filename_percent_ribo_sft_png = paste0(NORMDIR_DAT, "/", PREFIX, "-SFT-percent_ribo.png")
    print(filename_percent_ribo_sft_png)    
    plot = SpatialFeaturePlot(slice, features = "percent_ribo", pt.size.factor=1.2, crop=FALSE) + labs(title=name) + theme_f    
    ggsave(filename_percent_ribo_sft_png)

    ###############################
    # QC/Filteration
    print("QC/Filteration(1)")

    slice_wo_qc = slice
    
    slice = slice[, slice$nCount_Spatial > 500 & slice$nFeature_Spatial > 1000 & slice$nFeature_Spatial < 35000 & slice$percent_hb < 50]
    
    # Contents of the slice object
    filename_qc_slice_txt = paste0(NORMDIR_DAT, "/", PREFIX, "-QC-slice.txt")
    print(filename_qc_slice_txt)    
    sink(filename_qc_slice_txt)
    print(str(slice)) 
    sink()
    

    ###############################

    # VlnPlot(2)    
    # nCount_Spatial
    filename_qc_ncount_vln_png = paste0(NORMDIR_DAT, "/", PREFIX, "-QC-VLN-nCount_Spatial.png")
    print(filename_qc_ncount_vln_png)    
    plot = VlnPlot(slice, features = "nCount_Spatial", pt.size=0.5) + labs(title=name, y="nCount_Spatial") + NoLegend() + theme_v    
    ggsave(filename_qc_ncount_vln_png)

    # nFeature_Spatial
    filename_qc_nfeature_vln_png = paste0(NORMDIR_DAT, "/", PREFIX, "-QC-VLN-nFeature_Spatial.png")
    print(filename_qc_nfeature_vln_png)    
    plot = VlnPlot(slice, features = "nFeature_Spatial", pt.size=0.5) + labs(title=name, y="nFeature_Spatial") + NoLegend() + theme_v
    
    ggsave(filename_qc_nfeature_vln_png)

    # percent_mito
    filename_qc_percent_mito_vln_png = paste0(NORMDIR_DAT, "/", PREFIX, "-QC-VLN-percent_mito.png")
    print(filename_qc_percent_mito_vln_png)    
    plot = VlnPlot(slice, features = "percent_mito", pt.size=0.5) + labs(title=name, y="percent_mito") + NoLegend() + theme_v    
    ggsave(filename_qc_percent_mito_vln_png)

    # percent_hb
    filename_qc_percent_hb_vln_png = paste0(NORMDIR_DAT, "/", PREFIX, "-QC-VLN-percent_hb.png")
    print(filename_qc_percent_hb_vln_png)    
    plot = VlnPlot(slice, features = "percent_hb", pt.size=0.5) + labs(title=name, y="percent_hb") + NoLegend() + theme_v    
    ggsave(filename_qc_percent_hb_vln_png)

    # percent_ribo
    filename_qc_percent_ribo_vln_png = paste0(NORMDIR_DAT, "/", PREFIX, "-QC-VLN-percent_ribo.png")
    print(filename_qc_percent_ribo_vln_png)    
    plot = VlnPlot(slice, features = "percent_ribo", pt.size=0.5) + labs(title=name, y="percent_ribo") + NoLegend() + theme_v    
    ggsave(filename_qc_percent_ribo_vln_png)
    
    # FeatureScatter(2)
    # nCount_Spatial vs nFeature_Spatial
    filename_qc_fss_png = paste0(NORMDIR_DAT, "/", PREFIX, "-QC-FSS-nCount_Spatial-nFeature_Spatial.png")
    print(filename_qc_fss_png)    
    plot = FeatureScatter(slice, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial", pt.size=0.5) + labs(title=name) + NoLegend() + theme_s    
    ggsave(filename_qc_fss_png)

    # SpatialFeaturePlot(2)    
    # nCount_Spatial
    filename_qc_ncount_sft_png = paste0(NORMDIR_DAT, "/", PREFIX, "-QC-SFT-nCount_Spatial.png")
    print(filename_qc_ncount_sft_png)    
    plot = SpatialFeaturePlot(slice, features = "nCount_Spatial", pt.size.factor=1.3, crop=FALSE) + labs(title=name) + theme_f       
    ggsave(filename_qc_ncount_sft_png)

    # nFeature_Spatial
    filename_qc_nfeature_sft_png = paste0(NORMDIR_DAT, "/", PREFIX, "-QC-SFT-nFeature_Spatial.png")
    print(filename_qc_nfeature_sft_png)    
    plot = SpatialFeaturePlot(slice, features = "nFeature_Spatial", pt.size.factor=1.3, crop=FALSE) + labs(title=name) + theme_f    
    ggsave(filename_qc_nfeature_sft_png)

    # percent_mito
    filename_qc_percent_mito_sft_png = paste0(NORMDIR_DAT, "/", PREFIX, "-QC-SFT-percent_mito.png")
    print(filename_qc_percent_mito_sft_png)    
    plot = SpatialFeaturePlot(slice, features = "percent_mito", pt.size.factor=1.2, crop=FALSE) + labs(title=name) + theme_f
    ggsave(filename_qc_percent_mito_sft_png)

    # percent_hb
    filename_qc_percent_hb_sft_png = paste0(NORMDIR_DAT, "/", PREFIX, "-QC-SFT-percent_hb.png")
    print(filename_qc_percent_hb_sft_png)    
    plot = SpatialFeaturePlot(slice, features = "percent_hb", pt.size.factor=1.2, crop=FALSE) + labs(title=name) + theme_f    
    ggsave(filename_qc_percent_hb_sft_png)

    # percent_ribo
    filename_qc_percent_ribo_sft_png = paste0(NORMDIR_DAT, "/", PREFIX, "-QC-SFT-percent_ribo.png")
    print(filename_qc_percent_ribo_sft_png)    
    plot = SpatialFeaturePlot(slice, features = "percent_ribo", pt.size.factor=1.2, crop=FALSE) + labs(title=name) + theme_f    
    ggsave(filename_qc_percent_ribo_sft_png)
    
    slice = slice_wo_qc

    ###############################
    # Normalization
    print("Normalization")

    #####

    print("Log normalization")
    slice = NormalizeData(slice, verbose = FALSE, assay = "Spatial")

    # GroupCorrelation (for Seurat v4 only)
    if (FALSE) {

        filename_slice_txt = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-LOG-slice_before_GroupCorrelation.txt")
        print(filename_slice_txt)    
        sink(filename_slice_txt)
        print(str(slice)) 
        sink()

        print("GroupCorrelation")
        # LOG
        slice = GroupCorrelation(slice, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
    }

    # slice.txt
    filename_slice_txt = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-LOG-slice.txt")
    print(filename_slice_txt)    
    sink(filename_slice_txt)
    print(str(slice)) 
    sink()

    #####

    print("SCTransform normalization")

    slice = SCTransform(slice, assay = "Spatial", return.only.var.genes = FALSE, verbose = TRUE, clip.range=c(-20, 50))
        
    # GroupCorrelation (for Seurat v4 only)
    if (FALSE) {

        filename_slice_txt = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-slice_before_GroupCorrelation.txt")
        print(filename_slice_txt)    
        sink(filename_slice_txt)
        print(str(slice)) 
        sink()
        
        print("GroupCorrelation")
        # SCT
        slice = GroupCorrelation(slice, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)
    }

    # slice.txt
    filename_slice_txt = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-slice.txt")
    print(filename_slice_txt)    
    sink(filename_slice_txt)
    print(str(slice)) 
    sink()

    #####
    #####
    
    # GroupCorrelationPlot (for Seurat v4 only)

    if (FALSE) {
        print("GroupCorrelationPlot")
    
        # Plot LOG
        filename_gcp_log_png = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-LOG-GCP.png")
        print(filename_gcp_log_png)    
        p_log = GroupCorrelationPlot(slice, assay = "Spatial", cor = "nCount_Spatial_cor") +
            labs(title=paste0(name, " / ", "Log normalization"), x="Group") +
            theme(plot.title = element_text(size=20, hjust = 0.5),
                  axis.text.x=element_text(size=10),
                  axis.text.y=element_text(size=20),
                  axis.title.x=element_text(size=20),
                  axis.title.y=element_text(size=20),
                  aspect.ratio=1                            
                  )
        ggsave(filename_gcp_log_png)

        # Plot SCT
        filename_gcp_sct_png = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-GCP.png")
        print(filename_gcp_sct_png)        
        p_sct = GroupCorrelationPlot(slice, assay = "SCT", cor = "nCount_Spatial_cor") +
            labs(title=paste0(name, " / ", "SCT normalization"), x="Group") +            
            theme(plot.title = element_text(size=20, hjust = 0.5),
                  axis.text.x=element_text(size=10),
                  axis.text.y=element_text(size=20),
                  axis.title.x=element_text(size=20),
                  axis.title.y=element_text(size=20),
                  aspect.ratio=1                            
                  )
        ggsave(filename_gcp_sct_png)

    }
    
    ###############################

    # Contents of the slice object

    filename_slice_txt = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-slice.txt")
    print(filename_slice_txt)    
    sink(filename_slice_txt)
    print(str(slice)) 
    sink()

    # VlnPlot(3)    
    # nCount_SCT
    filename_sct_ncount_vln_png = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-VLN-nCount_SCT.png")
    print(filename_sct_ncount_vln_png)    
    plot = VlnPlot(slice, features = "nCount_SCT", pt.size=0.5) + labs(title=name, y="nCount_SCT") + NoLegend() + theme_v    
    ggsave(filename_sct_ncount_vln_png)

    # txt CSV
    df_nc_sct = select(slice@meta.data, nCount_SCT)
    print(head(df_nc_sct, 5))
    filename_ncount_sct_txt = paste0(NORMDIR_DAT, "/", PREFIX, "-nCount_SCT.txt")
    print(filename_ncount_sct_txt)
    
    write.table(df_nc_sct, file=filename_ncount_sct_txt, sep=",", na="", row.names=T, col.names=NA, quote=F)

    # nFeature_SCT
    filename_sct_nfeature_vln_png = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-VLN-nFeature_SCT.png")
    print(filename_sct_nfeature_vln_png)    
    plot = VlnPlot(slice, features = "nFeature_SCT", pt.size=0.5) + labs(title=name, y="nFeature_SCT") + NoLegend() + theme_v        
    ggsave(filename_sct_nfeature_vln_png)

    # txt
    df_nf_sct = select(slice@meta.data, nFeature_SCT)
    print(head(df_nf_sct, 5))
    filename_nfeature_sct_txt = paste0(NORMDIR_DAT, "/", PREFIX, "-nFeature_SCT.txt")
    print(filename_nfeature_sct_txt)

    write.table(df_nf, file=filename_nfeature_sct_txt, sep=",", na="", row.names=T, col.names=NA, quote=F)

    # percent_mito
    filename_sct_percent_mito_vln_png = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-VLN-percent_mito.png")
    print(filename_sct_percent_mito_vln_png)    
    plot = VlnPlot(slice, features = "percent_mito", pt.size=0.5) + labs(title=name, y="percent_mito") + NoLegend() + theme_v    
    ggsave(filename_sct_percent_mito_vln_png)

    # percent_hb
    filename_sct_percent_hb_vln_png = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-VLN-percent_hb.png")
    print(filename_sct_percent_hb_vln_png)    
    plot = VlnPlot(slice, features = "percent_hb", pt.size=0.5) + labs(title=name, y="percent_hb") + NoLegend() + theme_v    
    ggsave(filename_sct_percent_hb_vln_png)

    # percent_ribo
    filename_sct_percent_ribo_vln_png = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-VLN-percent_ribo.png")
    print(filename_sct_percent_ribo_vln_png)    
    plot = VlnPlot(slice, features = "percent_ribo", pt.size=0.5) + labs(title=name, y="percent_ribo") + NoLegend() + theme_v    
    ggsave(filename_sct_percent_ribo_vln_png)


    # FeatureScatter(3)
    # nCount_SCT vs nFeature_SCT
    filename_sct_fss_sct_sct_png = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-FSS-nCount_SCT-nFeature_SCT.png")
    print(filename_sct_fss_sct_sct_png)    
    plot = FeatureScatter(slice, feature1 = "nCount_SCT", feature2 = "nFeature_SCT", pt.size=0.5) + labs(title=name) + NoLegend() + theme_s    
    ggsave(filename_sct_fss_sct_sct_png)

    # nCount_SPatial vs nFeature_SCT
    filename_sct_fss_spatial_sct_png = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-FSS-nCount_Spatial-nFeature_SCT.png")
    print(filename_sct_fss_spatial_sct_png)    
    plot = FeatureScatter(slice, feature1 = "nCount_Spatial", feature2 = "nFeature_SCT", pt.size=0.5) + labs(title=name) + NoLegend() + theme_s    
    ggsave(filename_sct_fss_spatial_sct_png)

    # nCount_SPatial vs nCount_SCT
    filename_sct_fss_spatial_sct_png = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-FSS-nCount_Spatial-nCount_SCT.png")
    print(filename_sct_fss_spatial_sct_png)    
    plot = FeatureScatter(slice, feature1 = "nCount_Spatial", feature2 = "nCount_SCT", pt.size=0.5) + labs(title=name) + NoLegend() + theme_s    
    ggsave(filename_sct_fss_spatial_sct_png)

    
    # SpatialFeaturePlot(3)
    # nCount_SCT
    filename_sct_ncount_sft_png = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-SFT-nCount_SCT.png")
    print(filename_sct_ncount_sft_png)    
    plot = SpatialFeaturePlot(slice, features = "nCount_SCT", pt.size.factor=1.3, crop=FALSE) + labs(title=name) + theme_f   
    ggsave(filename_sct_ncount_sft_png)

    # nFeature_SCT
    filename_sct_nfeature_sft_png = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-SFT-nFeature_SCT.png")
    print(filename_sct_nfeature_sft_png)    
    plot = SpatialFeaturePlot(slice, features = "nFeature_SCT", pt.size.factor=1.3, crop=FALSE) + labs(title=name) + theme_f
    ggsave(filename_sct_nfeature_sft_png)

    # percent_mito
    filename_sct_percent_mito_sft_png = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-SFT-percent_mito.png")
    print(filename_sct_percent_mito_sft_png)    
    plot = SpatialFeaturePlot(slice, features = "percent_mito", pt.size.factor=1.2, crop=FALSE) + labs(title=name) + theme_f
    ggsave(filename_sct_percent_mito_sft_png)

    # percent_hb
    filename_sct_percent_hb_sft_png = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-SFT-percent_hb.png")
    print(filename_sct_percent_hb_sft_png)    
    plot = SpatialFeaturePlot(slice, features = "percent_hb", pt.size.factor=1.2, crop=FALSE) + labs(title=name) + theme_f    
    ggsave(filename_sct_percent_hb_sft_png)

    # percent_ribo
    filename_sct_percent_ribo_sft_png = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-SFT-percent_ribo.png")
    print(filename_sct_percent_ribo_sft_png)    
    plot = SpatialFeaturePlot(slice, features = "percent_ribo", pt.size.factor=1.2, crop=FALSE) + labs(title=name) + theme_f    
    ggsave(filename_sct_percent_ribo_sft_png)

    ##### #####
    # CSV files

    ##### SCT
    
    filename_sct_cc = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-counts.txt")    
    cc = as.data.frame(as.matrix(slice@assays$SCT@counts))
    fwrite(x=cc, row.names=TRUE, col.names=TRUE, file=filename_sct_cc)

    filename_sct_dd = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-data.txt")
    dd = as.data.frame(as.matrix(slice@assays$SCT@data))
    fwrite(x=dd, row.names=TRUE, col.names=TRUE, file=filename_sct_dd)

    filename_sct_sd = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-SCT-scale_data.txt")
    sd = as.data.frame(as.matrix(slice@assays$SCT@scale.data))
    fwrite(x=sd, row.names=TRUE, col.names=TRUE, file=filename_sct_sd)

    ##### Spatial

    filename_spatial_cc = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-Spatial-counts.txt")    
    cc_sp = as.data.frame(as.matrix(slice@assays$Spatial@counts))
    fwrite(x=cc_sp, row.names=TRUE, col.names=TRUE, file=filename_spatial_cc)

    filename_spatial_dd = paste0(NORMDIR_DAT, "/", PREFIX, "-NORM-Spatial-data.txt")
    dd_sp = as.data.frame(as.matrix(slice@assays$Spatial@data))
    fwrite(x=dd_sp, row.names=TRUE, col.names=TRUE, file=filename_spatial_dd)

    

    #####

    # Example

    default_assay = DefaultAssay(slice)
    print(default_assay)

    MakeSFP = function(gene_name) {

        print(paste0("R:", gene_name))

        data_ss = "SCT_counts"
        title_ss = paste0(name, " / ", "SCT@counts")
        filename_gene = paste0(NORMDIR_DAT, "/", PREFIX, "-SFP-", data_ss, "-", gene_name, ".png")          
        plot = SpatialFeaturePlot(slice, slot="counts", features=gene_name, pt.size.factor=1.3, crop=FALSE) + labs(title=title_ss) + theme_f
        ggsave(filename_gene)

        data_ss = "SCT_data"
        title_ss = paste0(name, " / ", "SCT@data")
        filename_gene = paste0(NORMDIR_DAT, "/", PREFIX, "-SFP-", data_ss, "-", gene_name, ".png")        
        plot = SpatialFeaturePlot(slice, slot="data", features=gene_name, pt.size.factor=1.3, crop=FALSE) + labs(title=title_ss) + theme_f   
        ggsave(filename_gene)

        data_ss = "SCT_scale_data"
        title_ss = paste0(name, " / ", "SCT@scale.data")
        filename_gene = paste0(NORMDIR_DAT, "/", PREFIX, "-SFP-", data_ss, "-", gene_name, ".png")            
        plot = SpatialFeaturePlot(slice, slot="scale.data", features=gene_name, pt.size.factor=1.3, crop=FALSE) + labs(title=title_ss) + theme_f   
        ggsave(filename_gene)
    
        #####

        DefaultAssay(slice) = "Spatial"
        default_assay = DefaultAssay(slice)
        print(default_assay)
    
        data_ss = "Spatial_counts"
        title_ss = paste0(name, " / ", "Spatial@counts")
        filename_gene = paste0(NORMDIR_DAT, "/", PREFIX, "-SFP-", data_ss, "-", gene_name, ".png")            
        plot = SpatialFeaturePlot(slice, slot="counts", features=gene_name, pt.size.factor=1.3, crop=FALSE) + labs(title=title_ss) + theme_f   
        ggsave(filename_gene)

        data_ss = "Spatial_data"
        title_ss = paste0(name, " / ", "Spatial@data")
        filename_gene = paste0(NORMDIR_DAT, "/", PREFIX, "-SFP-", data_ss, "-", gene_name, ".png")            
        plot = SpatialFeaturePlot(slice, slot="data", features=gene_name, pt.size.factor=1.3, crop=FALSE) + labs(title=title_ss) + theme_f   
        ggsave(filename_gene)

        DefaultAssay(slice) = "SCT"
        default_assay = DefaultAssay(slice)
        print(default_assay)

    }

    MakeSFP("MS4A1")
    MakeSFP("CD68")
    MakeSFP("CDH1")
    MakeSFP("TTF1")            
    MakeSFP("CD4")        
    MakeSFP("PTPRC") # CD45
    MakeSFP("APOE")        

    #####
    # Save RDS file

    filename_rds =paste0(NORMDIR_DAT, "/", PREFIX, "-", "SCT", ".rds")
    print(filename_rds)
    saveRDS(slice, filename_rds)
    
}

#####

main()







