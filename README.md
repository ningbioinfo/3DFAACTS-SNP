# 3DFAACTS-SNP

Here is the source code for the 3DFAACTS-SNP workflow from manuscript "**3DFAACTS-SNP: Using regulatory T cell-specific epigenomics data to uncover candidate mechanisms of Type-1 Diabetes (T1D) risk**", containing three major scripts:

1. _T1D_3DFAACTS-SNP_workflow.R_: to process the workflow and generate 3DFAACTS SNPs.
2. _HiCheatmapPlotting.R_: to plot the Hi-C matrix.
3. _GvizPlotting.R_: to plot multiple types of genomics data for visualisation of data integration.

## R packages dependencies:

- For running _T1D_3DFAACTS-SNP_workflow.R_:

  - tidyverse
  - GenomicRanges
  - org.Hs.eg.db
  - TxDb.Hsapiens.UCSC.hg19.knownGene
  - data.table
  - rtracklayer

- For running _HiCheatmapPlotting.R_:

  - Gviz
  - GenomicInteractions
  - coMET

### There are two python scripts are used to generate Hi-C data and ATAC-seq data that described in the manuscript, they require the following dependencies:

- pandas

### Data availability:

Most of the data that requires for the pipeline and scripts can be found in this repository. For ATAC-seq and Hi-C data, they are archived in ENA with accession id PRJEB39882.
