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

library("scales")

main = function() {

    args = commandArgs(trailingOnly=T)
    print(args)

    sid = args[1]
    filename_rds = args[2]
    CLSTDIR = args[3]

    PREFIX = paste0(sid, "-", "Visium")

    
    # slice1 with default options
    
    dirname_sr = paste0("./dat/Visium_", sid)
    slice1 = Load10X_Spatial(data.dir = dirname_sr)
    slice1 = SCTransform(slice1, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
    slice1 = RunPCA(slice1, assay="SCT", verbose=FALSE)
    slice1 = FindNeighbors(slice1, reduction="pca", dims=1:30)
    slice1 = FindClusters(slice1, verbose=FALSE)
    slice1 = RunUMAP(slice1, reduction="pca", dims=1:30)

    # SpatialDImPlot
    filename_plot1 = paste0(CLSTDIR, "/", PREFIX, "-SpatialDimPlot1.png" )
    print(filename_plot1)
    plot1 = SpatialDimPlot(slice1, crop=FALSE, pt.size.factor=1.0, group.by="ident", label=TRUE, label.size=3) + labs(title=sid) + 
        theme(axis.text.x=element_blank(), 
              axis.text.y=element_blank(), 
              axis.title.x=element_blank(), 
              axis.title.y=element_blank(), 
              plot.title=element_text(size=8),
              legend.title=element_blank(),
              legend.background=element_rect(fill=NA,colour=NA),
              legend.key=element_rect(fill=NA,colour=NA),
              )
    
    ggsave(filename_plot1)
    

    # Cluster IDs of spots txt
    df_clusters1 = select(slice1@meta.data, seurat_clusters)
    print(head(df_clusters1, 5))
    filename_seurat_clusters1_txt = paste0(CLSTDIR, "/", PREFIX, "-SCT_seurat_clusters1.txt")
    print(filename_seurat_clusters1_txt)
    write.table(df_clusters1, file=filename_seurat_clusters1_txt, sep=",", na="", row.names=T, col.names=NA, quote=F)
    
    # SCT@scale.data        
    filename_sct_sd1 = paste0(CLSTDIR, "/", PREFIX, "-NORM-SCT-scale_data1.txt")
    sd1 = as.data.frame(as.matrix(slice1@assays$SCT@scale.data))
    fwrite(x=sd1, row.names=TRUE, col.names=TRUE, file=filename_sct_sd1)
    
    seurat_clusters1_levels = levels(slice1@meta.data$seurat_clusters)
    cluster1_num = length(seurat_clusters1_levels)        
    
    filename_seurat_clusters1_rgb_txt = paste0(CLSTDIR, "/", PREFIX, "-SCT_seurat_clusters1_RGB.txt")
    print(hue_pal()(cluster1_num))
    
    write.table(hue_pal()(cluster1_num), file=filename_seurat_clusters1_rgb_txt, sep=",", na="", row.names=T, col.names=NA, quote=F)
    
    
    print(paste0("readRDS: ", filename_rds))
    slice = readRDS(filename_rds)
    
    # Dimensionality reduction
    print("Dimensionality reduction")
    slice = RunPCA(slice, assay="SCT", verbose=FALSE)
    slice = FindNeighbors(slice, reduction="pca", dims=1:30)
    slice = FindClusters(slice, resolution=1.0, verbose=FALSE)
    slice = RunUMAP(slice, reduction="pca", dims=1:30)

    # Contents of the slice object
    filename_slice_txt = paste0(CLSTDIR, "/", PREFIX, "-CLST-FindClusters-slice.txt")
    print(filename_slice_txt)    
    sink(filename_slice_txt)
    print(str(slice)) 
    sink()

    # Cluster IDs of spots txt
    df_clusters = select(slice@meta.data, seurat_clusters)
    
    print(head(df_clusters, 5))
    filename_seurat_clusters_txt = paste0(CLSTDIR, "/", PREFIX, "-SCT_seurat_clusters.txt")
    print(filename_seurat_clusters_txt)
    
    write.table(df_clusters, file=filename_seurat_clusters_txt, sep=",", na="", row.names=T, col.names=NA, quote=F)

    #####

    # factor type
    seurat_clusters_levels = levels(slice@meta.data$seurat_clusters)
    print("LEVELS")
    print(seurat_clusters_levels)

    cluster_num = length(seurat_clusters_levels)

    print("CLUSTER_NUM")
    print(cluster_num)

    filename_seurat_clusters_rgb_txt = paste0(CLSTDIR, "/", PREFIX, "-SCT_seurat_clusters_RGB.txt")
    print(hue_pal()(cluster_num))

    write.table(hue_pal()(cluster_num), file=filename_seurat_clusters_rgb_txt, sep=",", na="", row.names=T, col.names=NA, quote=F)

    #show_col(hue_pal()(cluster_num))    
    
    filename_plot = paste0(CLSTDIR, "/", PREFIX, "-DimPlot_UMAP.png" )
    print(filename_plot)
    plot = DimPlot(slice, reduction="umap", group.by="ident", label=TRUE) + labs(title=sid) + 
        theme(axis.text.x=element_text(size=10),
              axis.text.y=element_text(size=10),
              axis.title.x=element_text(size=10),
              axis.title.y=element_text(size=10),
              plot.title=element_text(size=8),
              legend.title=element_blank(),
              legend.background=element_rect(fill=NA,colour=NA),
              legend.key=element_rect(fill=NA,colour=NA)                            
              )
    ggsave(filename_plot)


    # Plot Clusters
    
    filename_plot = paste0(CLSTDIR, "/", PREFIX, "-SpatialDimPlot.png" )
    print(filename_plot)
    plot = SpatialDimPlot(slice, crop=FALSE, pt.size.factor=1.0, group.by="ident", label=TRUE, label.size=3) + labs(title=sid) + 
        theme(axis.text.x=element_blank(), 
              axis.text.y=element_blank(), 
              axis.title.x=element_blank(), 
              axis.title.y=element_blank(), 
              plot.title=element_text(size=8),
              legend.title=element_blank(),
              legend.background=element_rect(fill=NA,colour=NA),
              legend.key=element_rect(fill=NA,colour=NA),
              )
    ggsave(filename_plot)    


    # DEGs (FindMarkers)
    print("DEGs (FindMarkers)")

    filename_degs_txt = paste0(CLSTDIR, "/", PREFIX, "-SpatialFeaturePlot_FM_filenames.txt" )

    print("Filename\tgene\tp_val\tp_val_adj\tavg_log2FC")

    filename_v = c()
    name_v = c()
    pval_v = c()
    pval_adj_v = c()
    log2fc_v = c()

    k_v = c()
    n_v = c()

    for (k in 0:(cluster_num-1)) {

        de_markers <- FindMarkers(slice, ident.1=k)

        print(head(de_markers))
        

        print(rownames(de_markers)[1])

        rownames = rownames(de_markers)
        colnames = colnames(de_markers)

        for (n in 1:10) {
            
            name = rownames[n]
            p_val = de_markers[n, ]$p_val
            avg_log2FC = de_markers[n, ]$avg_log2FC
            p_val_adj = de_markers[n, ]$p_val_adj            

            print(paste("P_VAL:", p_val))
            print(paste("P_VAL_ADJ:", p_val_adj))
            print(paste("AVG_LOG2FC:", avg_log2FC))

            filename_plot = paste0(CLSTDIR, "/", PREFIX, "-SpatialFeaturePlot_FM", "_", sprintf("%02d", k), "-", sprintf("%02d", n), "_", name, ".png" )
            filename_plot_bb = paste0(PREFIX, "-SpatialFeaturePlot_FM", "_", sprintf("%02d", k), "-", sprintf("%02d", n), "_", name, ".png" )
            print(filename_plot)

            plot = SpatialFeaturePlot(slice, crop=FALSE, pt.size.factor=1.0, features=name, alpha=c(0.1, 1)) + 
                theme(legend.position="right", 
                      plot.title=element_text(size=8),            
                      legend.title=element_text(size=5),
                      legend.text=element_text(size=5),
                      legend.justification="center"
                      )
            ggsave(filename_plot)

            print(paste(filename_plot, name, p_val, p_val_adj, avg_log2FC))

            k_v = append(k_v, k)
            n_v = append(n_v, n)
            
            filename_v = append(filename_v, filename_plot_bb)
            name_v = append(name_v, name)
            pval_v = append(pval_v, p_val)
            pval_adj_v = append(pval_adj_v, p_val_adj)
            log2fc_v = append(log2fc_v, avg_log2FC)

        }

        dd = data.frame(cluster=k_v, n=n_v, filename=filename_v, gene=name_v, p_val=pval_v, p_val_adj=pval_adj_v, avg_log2FC=log2fc_v)
        print(dd)

        write.table(dd, filename_degs_txt, sep="\t", row.names=FALSE, quote=FALSE)

    }
    
}


main()

