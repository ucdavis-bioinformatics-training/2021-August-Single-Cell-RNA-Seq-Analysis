---
title: "Introduction to Single Cell RNAseq Part 4"
author: "UCD Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---

Last Updated: March 23 2021, 5pm

# Part 4: PCA and choice in number of PCS

## Load libraries

```r
library(Seurat)
library(ggplot2)
```

## Load the Seurat object

```r
load(file="pre_sample_corrected.RData")
experiment.aggregate
```

<div class='r_output'> An object of class Seurat 
 36601 features across 4000 samples within 1 assay 
 Active assay: RNA (36601 features, 3783 variable features)
</div>

## Now doing so for 'real'

ScaleData - Scales and centers genes in the dataset. If variables are provided in vars.to.regress, they are individually regressed against each gene, and the resulting residuals are then scaled and centered unless otherwise specified. Here we regress out cell cycle results S.Score and G2M.Score, percentage mitochondria (percent.mito) and the number of features (nFeature_RNA).


```r
experiment.aggregate <- ScaleData(
  object = experiment.aggregate,
  vars.to.regress = c("batchid", "S.Score", "G2M.Score", "percent.mito", "nFeature_RNA"))
```

## Dimensionality reduction with PCA

Next we perform PCA (principal components analysis) on the scaled data.  


```r
?RunPCA
```


```r
experiment.aggregate <- RunPCA(object = experiment.aggregate, npcs=100)
```

Seurat then provides a number of ways to visualize the PCA results

Visualize PCA loadings

```r
VizDimLoadings(experiment.aggregate, dims = 1, ncol = 1) + theme_minimal(base_size = 8)
```

![](scRNA_Workshop-PART4_files/figure-html/viz_pca-1.png)<!-- -->

```r
VizDimLoadings(experiment.aggregate, dims = 2, ncol = 1) + theme_minimal(base_size = 8)
```

![](scRNA_Workshop-PART4_files/figure-html/viz_pca-2.png)<!-- -->

Principal components plot

```r
DimPlot(object = experiment.aggregate, reduction = "pca")
```

![](scRNA_Workshop-PART4_files/figure-html/plot_pca-1.png)<!-- -->

Draws a heatmap focusing on a principal component. Both cells and genes are sorted by their principal component scores. Allows for nice visualization of sources of heterogeneity in the dataset.


```r
DimHeatmap(object = experiment.aggregate, dims = 1:6, cells = 500, balanced = TRUE)
```

![](scRNA_Workshop-PART4_files/figure-html/heatmap_pca-1.png)<!-- -->

```r
DimHeatmap(object = experiment.aggregate, dims = 7:12, cells = 500, balanced = TRUE)
```

![](scRNA_Workshop-PART4_files/figure-html/heatmap_pca-2.png)<!-- -->

#### Questions

1. Go back to the original data (rerun the load RData section) and then try modifying the ScaleData vars.to.regres, remove some variables, try adding in orig.ident? See how choices effect the pca plot

### Selecting which PCs to use
To overcome the extensive technical noise in any single gene, Seurat clusters cells based on their PCA scores, with each PC essentially representing a metagene that combines information across a correlated gene set. Determining how many PCs to include downstream is therefore an important step.

ElbowPlot plots the standard deviations (or approximate singular values if running PCAFast) of the principle components for easy identification of an elbow in the graph. This elbow often corresponds well with the significant PCs and is much faster to run.  This is the traditional approach to selecting principal components.


```r
ElbowPlot(experiment.aggregate, ndims = 0:100)
```

![](scRNA_Workshop-PART4_files/figure-html/elbow-1.png)<!-- -->

The JackStraw function randomly permutes a subset of data, and calculates projected PCA scores for these 'random' genes, then compares the PCA scores for the 'random' genes with the observed PCA scores to determine statistical signifance. End result is a p-value for each gene's association with each principal component. We identify significant PCs as those who have a strong enrichment of low p-value genes.


```r
experiment.aggregate <- JackStraw(object = experiment.aggregate, dims = 100)
```


```r
experiment.aggregate <- ScoreJackStraw(experiment.aggregate, dims = 1:100)
JackStrawPlot(object = experiment.aggregate, dims = 1:100) + theme(legend.position="bottom")
```

![](scRNA_Workshop-PART4_files/figure-html/plot_jackstraw-1.png)<!-- -->

## Finally, lets save the filtered and normalized data

```r
save(experiment.aggregate, file="pca_sample_corrected.RData")
```

## Get the next Rmd file

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2021-August-Single-Cell-RNA-Seq-Analysis/master/data_analysis/scRNA_Workshop-PART5.Rmd", "scRNA_Workshop-PART5.Rmd")
```

## Session Information

```r
sessionInfo()
```

<div class='r_output'> R version 4.0.3 (2020-10-10)
 Platform: x86_64-apple-darwin17.0 (64-bit)
 Running under: macOS Big Sur 10.16
 
 Matrix products: default
 BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
 LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
 
 locale:
 [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
 
 attached base packages:
 [1] stats     graphics  grDevices utils     datasets  methods   base     
 
 other attached packages:
 [1] ggplot2_3.3.3      SeuratObject_4.0.0 Seurat_4.0.1      
 
 loaded via a namespace (and not attached):
   [1] Rtsne_0.15            colorspace_2.0-0      deldir_0.2-10        
   [4] ellipsis_0.3.1        ggridges_0.5.3        spatstat.data_2.1-0  
   [7] leiden_0.3.7          listenv_0.8.0         farver_2.1.0         
  [10] ggrepel_0.9.1         fansi_0.4.2           codetools_0.2-18     
  [13] splines_4.0.3         knitr_1.31            polyclip_1.10-0      
  [16] jsonlite_1.7.2        ica_1.0-2             cluster_2.1.1        
  [19] png_0.1-7             uwot_0.1.10           shiny_1.6.0          
  [22] sctransform_0.3.2     spatstat.sparse_2.0-0 compiler_4.0.3       
  [25] httr_1.4.2            assertthat_0.2.1      Matrix_1.3-2         
  [28] fastmap_1.1.0         lazyeval_0.2.2        later_1.1.0.1        
  [31] htmltools_0.5.1.1     tools_4.0.3           igraph_1.2.6         
  [34] gtable_0.3.0          glue_1.4.2            RANN_2.6.1           
  [37] reshape2_1.4.4        dplyr_1.0.5           Rcpp_1.0.6           
  [40] scattermore_0.7       jquerylib_0.1.3       vctrs_0.3.6          
  [43] nlme_3.1-152          lmtest_0.9-38         xfun_0.22            
  [46] stringr_1.4.0         globals_0.14.0        mime_0.10            
  [49] miniUI_0.1.1.1        lifecycle_1.0.0       irlba_2.3.3          
  [52] goftest_1.2-2         future_1.21.0         MASS_7.3-53.1        
  [55] zoo_1.8-9             scales_1.1.1          spatstat.core_2.0-0  
  [58] promises_1.2.0.1      spatstat.utils_2.1-0  parallel_4.0.3       
  [61] RColorBrewer_1.1-2    yaml_2.2.1            reticulate_1.18      
  [64] pbapply_1.4-3         gridExtra_2.3         sass_0.3.1           
  [67] rpart_4.1-15          stringi_1.5.3         highr_0.8            
  [70] rlang_0.4.10          pkgconfig_2.0.3       matrixStats_0.58.0   
  [73] evaluate_0.14         lattice_0.20-41       ROCR_1.0-11          
  [76] purrr_0.3.4           tensor_1.5            patchwork_1.1.1      
  [79] htmlwidgets_1.5.3     labeling_0.4.2        cowplot_1.1.1        
  [82] tidyselect_1.1.0      parallelly_1.24.0     RcppAnnoy_0.0.18     
  [85] plyr_1.8.6            magrittr_2.0.1        R6_2.5.0             
  [88] generics_0.1.0        DBI_1.1.1             pillar_1.5.1         
  [91] withr_2.4.1           mgcv_1.8-34           fitdistrplus_1.1-3   
  [94] survival_3.2-10       abind_1.4-5           tibble_3.1.0         
  [97] future.apply_1.7.0    crayon_1.4.1          KernSmooth_2.23-18   
 [100] utf8_1.2.1            spatstat.geom_2.0-1   plotly_4.9.3         
 [103] rmarkdown_2.7         grid_4.0.3            data.table_1.14.0    
 [106] digest_0.6.27         xtable_1.8-4          tidyr_1.1.3          
 [109] httpuv_1.5.5          munsell_0.5.0         viridisLite_0.3.0    
 [112] bslib_0.2.4
</div>