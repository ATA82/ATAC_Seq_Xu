# Bulk ATAC-Seq analysis on CD13+/CD24+ sorted cells from human adult kidney tissue.

## General:
6 ATAC-Seq libraries generated from samples extracted from the kidney of one female renal cell carcinoma patient:
1. **FS3.CD13_R1**: CD13+ sorted cells from human adult kidney tissue.
2. **FS3.CD13_R2**: CD13+ sorted cells from human adult kidney tissue.
3. **FS3.CD13_R3**: CD13+ sorted cells from human adult kidney tissue.
4. **FS3.CD24_R1**: CD24+ sorted cells from human adult kidney tissue.
5. **FS3.CD24_R2**: CD24+ sorted cells from human adult kidney tissue.
6. **FS3.CD24_R3**: CD24+ sorted cells from human adult kidney tissue.

## Description of Processed Data:
Available at https://figshare.com/s/728705bc42446275044d

_Extended Peak Table:_ 
A table in an extended narrow peak format containing all peaks overlapping 1KB region aroud TSS.
_Count Table:_
contains for each of the 6 samples counts of the reads mapping on a region overlapping the selected peaks.

Processed gene expression values from the single-cell RNA-seq are https://doi.org/10.6084/m9.figshare.11786238

## Method Description
The ATAC-Seq data analysis was performed as follows: Generation of fastq files and adapter removal were completed with the Illumina software bcl2fastq (v2.20) [https://support.illumina.com/sequencing/sequencing_software/
bcl2fastq-conversion-software.html]. Library complexity was evaluated with Preseq (version 2.0) [1]. Sequence quality control was performed with FastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc). Additionally, MultiQC (version 1.7) [2] was used to aggregate the FastQC results across samples (See Supplementary ??). Subsequently, all reads were aligned to the human reference genome (GRCh38/hg38) using bwa (version 0.7.17) [3] with default parameters. A clear overrepresentation of mitochondrial reads, especially in CD13, was observed using Rsamtools [4] (version 2.0.0). All mitochondrial reads were filtered out using samtools [5] (version 0.1.19) and followed up by additional sequence quality control. Afterwards, peak calling was performed using common software MACS2 (version 2.2.5) [6] with the following parameter (–keep-dup all) and otherwise default parameters for paired-end reads. BigWig track files for visualisation of peaks in genome browsers were generated using the rtracklayer package (version 1.44.2) [7] and normalized with 1 Million as scaling factor. After generating a merged peak file from all samples, peaks annotation was performed using ChIPseeker (version 1.20.0) [8] with standard parameters. For downstream analysis only peaks overlaping a 1KB region around a TSS were selected. The featureCounts software from the Rsubread package (version 1.34.7) [9] was used to count reads mapping on a region overlaping the selected peaks. Analysis of differential chromatine accessibility was performed using DESeq2 (version 1.24.0) [10] with the “local” fit type and 10-4 as minimal dispersion value for gene dipersion estimation. For PCA-analysis and heatmap generation raw counts were normalized using the variance stabilization transformation approach of DESeq2 to take account for library size and biases introduced by low counts genes. The piano R package (version 2.0.2) [11] was used to analyze the Gene Set Enrichment of GO-terms. The GO-terms of the biological processes were selected from the Molecular Signatures Database (MsigDB) (version 6.0) [12]. The gene sets p-values were computed based on the consensus approach of the piano package applied on the fisher, stouffer and reporter tests with a consensus cutoff of 50 using the mean-method. The null distribution was selected as method for significance assessment of gene sets in all three tests and the Benjamini-Hochberg correction method was applied for multiple-tests p-value adjustment. 

    (1) Daley, T. & Smith, A. D. Predicting the molecular complexity of sequencing libraries. Nat. Methods 10, 325–327 (2013).
    (2) Ewels, P., Magnusson, M., Lundin, S. & Kaller, M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics 32, 3047–3048 (2016).
    (3) Li, H. & Durbin, R. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics 25, 1754–1760 (2009).
    (4) Lawrence, M. & Morgan, M. Scalable genomics with R and Bioconductor. Stat. Sci. 29, 214–226 (2014).
    (5) Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078–2079 (2009).
    (6) Feng, J., Liu, T., Qin, B., Zhang, Y. & Liu, X. S. Identifying ChIP–seq enrichment using MACS. Nat. Protoc. 7, 1728–1740 (2012).
    (7) Lawrence, M., Gentleman, R. & Carey, V. rtracklayer: an R package for interfacing with genome browsers. Bioinformatics 25, 1841–1842 (2009).
    (8) Yu, G., Wang, L. G. & He, Q. Y. ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. Bioinformatics 31, 2382–2383 (2015).
    (9) Liao, Y., Smyth, G. K. & Shi, W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics 30, 923–930 (2014).
    (10) Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 15, 550 (2014).
    (11) Väremo, L., Nielsen, J. & Nookaew, I. Enriching the gene set analysis of genome-wide data by incorporating directionality of gene expression and combining statistical hypotheses and methods. Nucleic Acids Res. 41, 4378–4391 (2013).
    (12) Liberzon, A. et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst. 1, 417–425 (2015).

## Usage
The R markdown script called Xu_Tubuloid.Rmd contains the main pipeline for the ATAC-Seq Analysis in this paper. Addionally generating the extended merged peak file was performed using a java class included in this repository and the command mergeBed from bedtools. This is also described in the Rmd file. 

The R markdown script contains also the code chunks used to generate the ATAC-Seq plots integrated in the figures of this paper. Only the track figures were created separately with IGV.

The analysis is based on 6 ATAC-Seq libraries based on Samples from one female renal cell carcinoma patient:
FS3.CD13_R1: CD13+ sorted cells from human adult kidney tissue.
FS3.CD13_R2: CD13+ sorted cells from human adult kidney tissue.
FS3.CD13_R3: CD13+ sorted cells from human adult kidney tissue.
FS3.CD24_R1: CD24+ sorted cells from human adult kidney tissue.
FS3.CD24_R2: CD24+ sorted cells from human adult kidney tissue.
FS3.CD24_R3: CD24+ sorted cells from human adult kidney tissue.

[Optionally] Request access to the raw sequencing 10x data (FastQ files) to Rafael Kramann following the manuscript details. Then preprocess the data following your preferences. The processed data with peaks tables of each sample alongside of the merge version the reads counts table of the peaks overlaping 1KB aroud TSS region are available at https://figshare.com/s/728705bc42446275044d.

## Environment
    
Package Name | Source | Version |
------------ | ------------- | ------------- |
blc2fastq | https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html | 2.20 |
Preseq | http://smithlabresearch.org/software/preseq/ | 2.0 |
FastQC | http://www.bioinformatics.babraham.ac.uk/projects/fastqc) |0.11.8|
MultiQC | https://multiqc.info/ | 1.7 |
bwa | bio-bwa.sourceforge.net | 0.7.17 |
Rsamtools | bioconductor | 2.0.0 |
samtools | http://www.htslib.org/download/ | 0.1.19 |
macs2 | https://github.com/taoliu/MACS | 2.2.5 |
rtracklayer | bioconductor | 1.44.2 |
ChIPseeker | bioconductor | 1.20.0 |
Rsubread | bioconductor | 1.34.7 |
DESeq2 | bioconductor | 1.24.0 |
piano | bioconductor | 2.0.2 |


Detailed environment for reproducibility

Herein it is described the sessionInfo() from R:
```
R version 3.6.0 (2019-04-26)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=de_DE.UTF-8       LC_NUMERIC=C               LC_TIME=de_DE.UTF-8        LC_COLLATE=de_DE.UTF-8     LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=de_DE.UTF-8    LC_PAPER=de_DE.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=de_DE.UTF-8
[12] LC_IDENTIFICATION=C       

attached base packages:
 [1] grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] knitr_1.24                              snowfall_1.84-6.1                       snow_0.4-3                              icesTAF_3.3-1                           piano_2.0.2                             pander_0.6.3                            dplyr_0.8.3                            
 [8] gdata_2.18.0                            Hmisc_4.2-0                             Formula_1.2-3                           survival_2.44-1.1                       edgeR_3.26.8                            limma_3.40.6                            ggplotify_0.0.4                        
[15] pheatmap_1.0.12                         TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.6 GenomicFeatures_1.36.4                  org.Hs.eg.db_3.8.2                      AnnotationDbi_1.46.0                    clusterProfiler_3.12.0                 
[22] ChIPseeker_1.20.0                       scales_1.0.0                            ATACseqQC_1.8.5                         lattice_0.20-38                         gridExtra_2.3                           Rsamtools_2.0.0                         DESeq2_1.24.0                          
[29] SummarizedExperiment_1.14.1             DelayedArray_0.10.0                     BiocParallel_1.18.1                     matrixStats_0.54.0                      Biobase_2.44.0                          ggrepel_0.8.1                           ggpubr_0.2.3                           
[36] magrittr_1.5                            ggplot2_3.2.1                           drc_3.0-1                               MASS_7.3-51.4                           keras_2.2.5.0                           tensorflow_2.0.0                        RColorBrewer_1.1-2                     
[43] fields_10.0                             maps_3.3.0                              spam_2.5-1                              dotCall64_1.0-0                         plotrix_3.7-6                           kohonen_3.0.10                          BSgenome_1.52.0                        
[50] Biostrings_2.52.0                       XVector_0.24.0                          rtracklayer_1.44.2                      GenomicRanges_1.36.0                    GenomeInfoDb_1.20.0                     IRanges_2.18.1                          S4Vectors_0.22.0                       
[57] RSQLite_2.1.2                           jsonlite_1.6                            shiny_1.3.2                             BiocGenerics_0.30.0                     seqplots_1.22.2                        

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.1                tidyr_0.8.3                   acepack_1.4.1                 bit64_0.9-7                   multcomp_1.4-10               data.table_1.12.2             rpart_4.1-15                  RCurl_1.95-4.12               AnnotationFilter_1.8.0       
 [10] generics_0.0.2                cowplot_1.0.0                 lambda.r_1.2.3                TH.data_1.0-10                europepmc_0.3                 bit_1.1-14                    enrichplot_1.4.0              xml2_1.2.2                    httpuv_1.5.1                 
 [19] assertthat_0.2.1              viridis_0.5.1                 xfun_0.8                      hms_0.5.0                     evaluate_0.14                 promises_1.0.1                progress_1.2.2                caTools_1.17.1.2              dbplyr_1.4.2                 
 [28] readxl_1.3.1                  igraph_1.2.4.1                DBI_1.0.0                     geneplotter_1.62.0            htmlwidgets_1.3               futile.logger_1.4.3           purrr_0.3.2                   backports_1.1.4               GenomicScores_1.8.1          
 [37] grImport2_0.1-5               annotate_1.62.0               gridBase_0.4-7                biomaRt_2.40.3                vctrs_0.2.0                   ensembldb_2.8.1               abind_1.4-5                   withr_2.1.2                   ggforce_0.3.1                
 [46] Gviz_1.28.3                   triebeard_0.3.0               checkmate_1.9.4               rGADEM_2.32.0                 GenomicAlignments_1.20.1      prettyunits_1.0.2             cluster_2.1.0                 DOSE_3.10.2                   seqLogo_1.50.0               
 [55] lazyeval_0.2.2                crayon_1.3.4                  relations_0.6-9               genefilter_1.66.0             slam_0.1-46                   pkgconfig_2.0.2               labeling_0.3                  tweenr_1.0.1                  ProtGenerics_1.16.0          
 [64] nnet_7.3-12                   rlang_0.4.0                   sandwich_2.5-1                seqinr_3.6-1                  BiocFileCache_1.8.0           MotIV_1.40.0                  AnnotationHub_2.16.1          dichromat_2.0-0               VennDiagram_1.6.20           
 [73] randomForest_4.6-14           cellranger_1.1.0              polyclip_1.10-0               graph_1.62.0                  Matrix_1.2-17                 urltools_1.7.3                carData_3.0-2                 boot_1.3-23                   zoo_1.8-6                    
 [82] Rsubread_1.34.7               base64enc_0.1-3               whisker_0.4                   ggridges_0.5.1                png_0.1-7                     viridisLite_0.3.0             shinydashboard_0.7.1          bitops_1.0-6                  visNetwork_2.0.8             
 [91] KernSmooth_2.23-15            blob_1.2.0                    stringr_1.4.0                 qvalue_2.16.0                 regioneR_1.16.5               jpeg_0.1-8.1                  gridGraphics_0.4-1            ggsignif_0.6.0                memoise_1.1.0                
[100] plyr_1.8.4                    gplots_3.0.1.1                bibtex_0.4.2                  zlibbioc_1.30.0               compiler_3.6.0                ade4_1.7-13                   htmlTable_1.13.2              formatR_1.7                   tidyselect_0.2.5             
[109] stringi_1.4.3                 forcats_0.4.0                 yaml_2.2.0                    GOSemSim_2.10.0               locfit_1.5-9.1                ChIPpeakAnno_3.18.2           latticeExtra_0.6-28           VariantAnnotation_1.30.1      polynom_1.4-0                
[118] fastmatch_1.1-0               tools_3.6.0                   rio_0.5.16                    rstudioapi_0.10               foreign_0.8-72                idr_1.2                       farver_1.1.0                  ggraph_1.0.2                  digest_0.6.20                
[127] rvcheck_0.1.3                 BiocManager_1.30.4            Rcpp_1.0.2                    car_3.0-4                     later_0.8.0                   motifStack_1.28.0             httr_1.4.1                    biovizBase_1.32.0             colorspace_1.4-1             
[136] XML_3.98-1.20                 reticulate_1.13               splines_3.6.0                 RBGL_1.60.0                   multtest_2.40.0               preseqR_4.0.0                 xtable_1.8-4                  futile.options_1.0.1          marray_1.62.0                
[145] UpSetR_1.4.0                  zeallot_0.1.0                 R6_2.4.0                      sets_1.0-18                   pillar_1.4.2                  htmltools_0.3.6               mime_0.7                      glue_1.3.1                    DT_0.8                       
[154] interactiveDisplayBase_1.22.0 class_7.3-15                  codetools_0.2-16              fgsea_1.10.0                  mvtnorm_1.0-11                tibble_2.1.3                  curl_4.0                      tfruns_1.4                    gtools_3.8.1                 
[163] shinyjs_1.0                   zip_2.0.4                     GO.db_3.8.2                   openxlsx_4.1.2                rmarkdown_1.14                munsell_0.5.0                 DO.db_2.9                     GenomeInfoDbData_1.2.1        haven_2.1.1                  
[172] reshape2_1.4.3                gtable_0.3.0    

```
