# This script is the source code describing the 3DFAACTS-SNP workflow from the manuscript "3DFAACTS-SNP: Using regulatory T cell-specific epigenomics data to uncover candidate mechanisms of Type-1 Diabetes (T1D) risk".
# Author: Ning Liu
# GitHub repo: https://github.com/ningbioinfostruggling/3DFAACTS-SNP

library(tidyverse)
library(GenomicRanges)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(data.table)
library(rtracklayer)
setwd('~/Documents/T1DPaper/PublicCodes/')

########## load data
t1d.snps <- fread('./T1D_finemappedSNPs_Onengut2015.txt') %>%
  as_tibble() %>%
  magrittr::set_colnames(c("chr","rsid","start","MAF")) %>%
  .[complete.cases(.),] %>%
  mutate(end = start) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  sort() %>%
  unique()

atac_treg <- import.bed("Treg_all_ATAC_summits_collapsed_1kb.bed")

foxp3 <- import.bed("./GSE20995_Treg.ab.HsTiling2.0.FDR0.5.hg19.bed")

treg_hic <- read_delim('./Treg_HiC_cis.bedpe', 
                       delim = '\t', col_names = F) %>%
  magrittr::set_colnames(c("chr1","start1","end1",
                           "chr2","start2","end2","count")) %>%
  filter(chr1==chr2)

hg19_genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene) %>%
  as_tibble() %>%
  left_join(., as.data.frame(org.Hs.egSYMBOL)) %>%
  dplyr::select(-gene_id) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  sort()

txid <- keys(TxDb.Hsapiens.UCSC.hg19.knownGene, "TXID")

df <- select(TxDb.Hsapiens.UCSC.hg19.knownGene, txid, "GENEID", "TXID")
df <- cbind(df,select(org.Hs.eg.db, df$GENEID, c("SYMBOL"))$SYMBOL) %>%
  .[complete.cases(.),] %>%
  magrittr::set_colnames(c("tx_id","ENTREZID","symbol")) %>%
  mutate(tx_id = as.integer(tx_id)) %>%
  as_tibble()

promoters <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene) %>% 
  as_tibble() %>% 
  filter(tx_id %in% df$tx_id) %>%
  left_join(df) %>%
  dplyr::select(-c(tx_name)) %>%
  mutate(tss = ifelse(strand == "+", start, end )) %>% 
  mutate(pro_start = ifelse(strand == "+", tss-2000, tss), pro_end = ifelse(strand == "+", tss, tss+2000)) %>%
  dplyr::select(-c(start,end,strand,tss,width, ENTREZID,tx_id)) %>% 
  magrittr::set_colnames(c("seqnames","symbol","start","end")) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

enhancers.per <- read_delim("./FANTOM5_permissive_Enh.bed", delim = "\t", 
                            col_names = F) %>%
  .[1:3] %>%
  magrittr::set_colnames(c("chr","start","end")) %>%
  mutate(name = paste0("Enh_",seq_len(nrow(.)))) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  sort()

enhancers.tcell <- read_delim("./FANTOM5_Tcell_specific_Enh.bed", delim = "\t", col_names = F) %>%
  .[1:3] %>%
  magrittr::set_colnames(c("chr","start","end")) %>%
  mutate(name = paste0("TEnh_",seq_len(nrow(.)))) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  sort()

enhancers <- c(subsetByOverlaps(enhancers.per, enhancers.tcell, invert = T), enhancers.tcell) %>% sort()

chromhmm <- read_delim("./E044_15_coreMarks_dense.bed", delim = "\t", col_names = F) %>%
  .[c(1:4,9)] %>%
  magrittr::set_colnames(c("chr","start","end","states","col"))

chromhmm_promoter <- chromhmm %>%
  filter(grepl("Tss",states) | states == "11_BivFlnk") %>%
  dplyr::select(-col) %>%
  mutate(states = str_remove(states, "[0-9]*_")) %>%
  split(., f = .$states) %>%
  lapply(function(x) x %>% mutate(states = paste0(states,"_chmm_",seq(nrow(x))))) %>%
  bind_rows() %>%
  magrittr::set_colnames(c("chr","start","end","symbol")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  sort()

chromhmm_enhancer <- chromhmm %>%
  filter(grepl("Enh",states)) %>%
  dplyr::select(-col) %>%
  mutate(states = str_remove(states, "[0-9]*_")) %>%
  split(., f = .$states) %>%
  lapply(function(x) x %>% mutate(states = paste0(states,"_chmm_",seq(nrow(x))))) %>%
  bind_rows() %>%
  magrittr::set_colnames(c("chr","start","end","name")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  sort()

promoters <- c(promoters, chromhmm_promoter) %>% 
  sort() %>%
  as_tibble() %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

enhancers <- c(enhancers, chromhmm_enhancer) %>% 
  sort() %>% 
  as_tibble() %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

########## 3DFAACTS-SNP workflow

source('./iie.R')
source('./workflow3dfaacts.R')

t1d_3dfaacts_snp <- workflow_3DFAACTSsnp(t1d.snps, atac_treg, 
                                         treg_hic, promoters, 
                                         enhancers, foxp3)

