#!/usr/bin/Rscripts
library(dplyr)
library(data.table)
library(GenomicRanges)
library(parallel)
library(ggplot2)
library(tidyr)
library(reshape2)  

library(ggplot2)
library(tidyr)
library(pheatmap)
library(stringr)
library(magrittr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(reshape2)  
library("ggsci")
library("scales")

theme <- theme(panel.background = element_blank(),panel.border = element_rect(fill = NA),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      strip.background=element_blank(),axis.text.x=element_text(colour = "black"),
      axis.text.y = element_text(colour = "black"),axis.ticks = element_line(colour = "black"),
      plot.margin = unit(c(5,1,1,1),"line"))

#### validate aymmetric co-occupancy peaks in scNanoSeq-CUT&Tag ####

library(dplyr)
library(data.table)
library(GenomicRanges)
library(parallel)
library(vcfR)

###########################################################################################
###########################################################################################
# note: A/B is not pat/total, but alt/total ("1"/total)
vcf <- read.vcfR("./ENCODE_GM12878_H3K4me3_GRCh38_GIAB_highconf.het.vcf.gz")
dp <- extract.info(vcf,element = 'DP', as.numeric = T)
ad <- extract.info(vcf,element = 'AO', as.numeric = T)
ratio <- ad/dp

# generate a GRange object for ATAC-seq SNP
vcf_ref <- read.vcfR('./HG001_GRCh38_GIAB_highconf.het.vcf.gz')
gt_ref <-
   extract.gt(vcf_ref)[paste0(getCHROM(vcf), "_", getPOS(vcf)), ]

snp <-
   GRanges(
      seqnames = getCHROM(vcf),
      ranges = IRanges(start = getPOS(vcf), width = 1),
      dp = dp,
      ratio = ratio,
      gt_ref = gt_ref
   )
seqlevelsStyle(snp) <- 'UCSC'
snp <- snp[!is.na(snp$ratio)]
snp <- snp[snp$dp > 15]

# visualize
hist(snp$ratio, breaks = 100)
plot(snp$dp, snp$ratio, cex = .1)


# load long atac-seq biased peak call set
if (T) {
   long.asym <- readRDS('../H3K4me3/20231103_GM12878_H3K4me3_TGS_peaks_biased_20frag.Rds')
   
   df <-
      fread('../H3K4me3/GM12878_H3K4me3_sc_seacr_top0.05.peaks.stringent.17cell_support.bed',
            header = F,
            sep = '\t') %>% as.data.frame()
   names(df) <-
      c('chr', 'start', 'end')
   peak.gr <- GRanges(seqnames = df$chr,
                      ranges = IRanges(df$start, df$end))
   seqlevelsStyle(peak.gr) <- 'UCSC'
}

# validation of bias by NGS ATAC-seq SNP ratio
if (T) {
   # maternal or paternal specific asymmetric peaks
   mode <- 'mat'
   mode <- 'pat'
   
   long.cut.off <- 0.05
   ngs.cut.off <- 0.05
   
   if (mode == 'mat') {
      sig <- long.asym$fdr <= long.cut.off & long.asym$ratio <= .1
      peak.index <-
         long.asym$peak[sig] %>% unique
      minor.gt <- '0|1'
      
      # the count of autosome ASPs
      table(long.asym[sig,]$chr=='chrX')
   }
   
   if (mode == 'pat') {
      sig <- long.asym$fdr < long.cut.off & long.asym$ratio >= .9
      peak.index <-
         long.asym$peak[sig] %>% unique
      minor.gt <- '1|0'
      
      # the count of autosome ASPs
      table(long.asym[sig,]$chr=='chrX')
   } 
   
   hit <- findOverlapPairs(peak.gr[peak.index], snp)
   
   binom.p <- c()
   for (i in seq_along(hit@second)) {
      x <- hit@second[i]
      binom.p <- c(binom.p,
                   binom.test(round(x$ratio * x$dp), x$dp, p = 0.5)$p.value)
   }
   binom.p <- p.adjust(binom.p, method = 'fdr')
   
   ##### STATISTICS ####
   # validatable peaks
   hit@first %>% unique() %>% length
   # validation SNPs
   hit@second %>% unique() %>% length
   
   criteria <- (binom.p < ngs.cut.off)
   hit.peak <- hit@first[criteria]
   hit.snp <- hit@second[criteria]
   
   # gross precision (by peak)
   (hit.peak %>% unique %>% length) / (hit@first %>% unique %>% length)
   # gross precision (by SNP, if any SNP in the peak supports)
   (hit.snp %>% unique %>% length) / (hit@second %>% unique %>% length)
   
   # validate consistency of bias direction between long and short atac
   table(hit.snp$gt_ref == minor.gt,
         hit.snp$ratio > .5)
   valid <- !xor(hit.snp$gt_ref == minor.gt,
                 hit.snp$ratio > .5)
   # consistent precision in significant subset (by SNP)
   sum(valid) / length(valid)
   
   # overall precision (by peak)
   (hit.peak[valid] %>% unique %>% length) / 
      (hit@first %>% unique %>% length)
   
   # overall precision (by SNP)
   (hit.snp[valid] %>% unique %>% length) / 
      (hit@second %>% unique %>% length)
###################################################################################################################################   
   ###dont consider Pvalue,just consider ratio
   criteria <- (binom.p <100)
   hit.peak <- hit@first[criteria]
   hit.snp <- hit@second[criteria]
   # validate consistency of bias direction between long and short atac
   table(hit.snp$gt_ref == minor.gt,
         hit.snp$ratio > .5)
   valid <- !xor(hit.snp$gt_ref == minor.gt,
                 hit.snp$ratio > .5)
   # overall precision (by peak)
   (hit.peak[valid] %>% unique %>% length) / 
      (hit@first %>% unique %>% length)
   
   # overall precision (by SNP)
   (hit.snp[valid] %>% unique %>% length) / 
      (hit@second %>% unique %>% length)
   
}

