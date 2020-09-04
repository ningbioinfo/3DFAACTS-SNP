# This script is used for plotting Hi-C normalised matrix for the manuscript "3DFAACTS-SNP: Using regulatory T cell-specific epigenomics data to uncover candidate mechanisms of Type-1 Diabetes (T1D) risk".
# Author: Ning Liu

# In this script, we aim at plotting the Hi-C matrix in Figure S2 from Additional file 1 of the manuscript.

source("./hicheatmap.R")

chr <- "chr2"
start <- 204122714
end <- 204812714
hmstart <- 203922714
hmend <- 205092714

hic_m <- from_densematrix_to_hmplotting(paste0("./HiC_test_chr2.matrix"))
hicHeatmap(hic_m, chr, hmstart, hmend, zrange = c(0,3), max_y = 11, flip = F)
#plotTADs2HM(paste0("TADs/armatus_HiC_merged_40k_",chr,"_",chr,"_40k.tad.consensus.txt"), hmstart, hmend, 40000, tad_col = "green")
plotTADs2HM("./TAD_test.tsv", chr, hmstart, hmend, res = 40000)
plot_sub_triangle(start,end,40000)
