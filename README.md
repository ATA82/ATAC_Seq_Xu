# ATAC_Seq_Xu

## General:
6 ATAC-Seq samples from one female renal cell carcinoma patient:
FS3.CD13_R1: CD13+ sorted cells from human adult kidney tissue.
FS3.CD13_R2: CD13+ sorted cells from human adult kidney tissue.
FS3.CD13_R3: CD13+ sorted cells from human adult kidney tissue.
FS3.CD24_R1: CD24+ sorted cells from human adult kidney tissue.
FS3.CD24_R2: CD24+ sorted cells from human adult kidney tissue.
FS3.CD24_R3: CD24+ sorted cells from human adult kidney tissue.

## Description of Processed Data:
Available at https://figshare.com/s/728705bc42446275044d \n
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
