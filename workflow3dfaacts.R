workflow_3DFAACTSsnp <- function(snp,atac,hic,pro,enh,TF=0){
  require(GenomicRanges)
  require(tidyverse)
  
  snps.atac <- subsetByOverlaps(snp, atac)
  
  hic.t1d <- interaction_intersecting_elements(interaction = hic, 
                                               gr_elements = snps.atac)
  
  snps.ah <- y
  
  ov.pro <- findOverlaps(snps.ah, pro)
  
  ov.enh <- findOverlaps(snps.ah, enh)
  
  snps.ah.pro <- snps.ah[ov.pro@from]
  
  snps.ah.enh <- snps.ah[ov.enh@from]
  
  snps.ah.enh.pro <- c(snps.ah.pro,snps.ah.enh) %>% 
    sort() %>% 
    unique()
  
  snps.ahpe.foxp3 <- subsetByOverlaps(snps.ah.enh.pro, TF)
  
  return(snps.ahpe.foxp3)
}