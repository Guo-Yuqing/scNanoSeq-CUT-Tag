#!/usr/bin/env Rscripts
###
 # @Descripttion:scNanoSeq-CUT&Tag co-occupancy
 # @version:
 # @Date: 2023-6-24 22:35:11
 # @LastEditors:
 # @LastEditTime: 2023-6-24 22:35:20
###
argv <- commandArgs(TRUE)

cell.type <- argv[1]
modification <- argv[2]

setwd(paste0("co_occupancy/",cell.type ,"/coa_ks_test"))
######################################################################################################################################################################
library(data.table)
library(GenomicRanges)
library(dplyr)
library(genomation)
library(parallel)
library(data.table)
library(GenomicRanges)
library(dplyr)
library(InteractionSet)
library(GenomicInteractions)
library(genomation)

library(ggplot2)
library(PerformanceAnalytics)
library(MASS)
library(viridis)

library("ggsci")
library("scales")
######################################################################################################################################################################
##########################################################################################################################
###########################################################################################################################
# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# select pairs of co-accessible peaks, 
# make interaction objects, 
# and write them in bedpe files.
if (T) {
    
    peaks <-
        readBed(paste0(cell.type,"_",modification,"_peak.bed")) 

    gi.clean<-NULL
    coa.all<-NULL
    
    for (chrom in c(1:22,'X')){
      # get subset of peaks by chromosome
      peaks.chr <- subset(peaks,seqnames==paste0('chr',chrom))
      peak.test<-readRDS(paste0('scNanoCuttag_',cell.type,"_",modification,'_chr',chrom,'.co-occupancy.ks.Rds'))
      peak.test$region=paste0(peak.test$start,":",peak.test$end)

      # extract co-accessible peak pairs(reads support directly & reads length>1kb &  peak length <1kb)
      peak_pair <- read.table(paste0("../peakpair/chr",chrom,"_peakpairs.bed"),header = F,sep = '\t')
      peak_pair$left_region <- paste0(peak_pair$V2+1,":",peak_pair$V3)
      peak_pair$right_region <- paste0(peak_pair$V8+1,":",peak_pair$V9)
      
      peaks.chr.data <- as.data.frame(peaks.chr)
      peaks.chr.data$region <- paste0(peaks.chr.data$start,":",peaks.chr.data$end)
      peak_pair <- peak_pair[peak_pair$left_region %in% peaks.chr.data$region & peak_pair$right_region %in% peaks.chr.data$region,]
      
      for (i in 1:nrow(peak_pair)){
        
        peak_pair[i,"left_judge"] <- peak.test[peak.test$region==peak_pair[i,"left_region"],'left.p.val.adj']
        peak_pair[i,"right_judge"] <- peak.test[peak.test$region==peak_pair[i,"right_region"],'right.p.val.adj']
        peak_pair[i,"left_index"] <- as.numeric(rownames(peaks.chr.data[peaks.chr.data$region==peak_pair[i,"left_region"],]))
        peak_pair[i,"right_index"] <- as.numeric(rownames(peaks.chr.data[peaks.chr.data$region==peak_pair[i,"right_region"],]))
        
      }
      
      peak_pair <- peak_pair[peak_pair$left_judge<0.01 & peak_pair$right_judge<0.01  & peak_pair$V3-peak_pair$V2<=1500 &peak_pair$V9-peak_pair$V8<=1500 ,]
      peak_pair <- na.omit(peak_pair)


      sig.left.index <- peak_pair$left_index
      sig.right.index <- peak_pair$right_index
      
      # generate a GenomicInteraction object
      gi <- GInteractions(sig.left.index, sig.right.index, peaks.chr)
      
    
    export.bedpe(gi,
                 paste0("../coa_bedpe/",cell.type,"_",modification,'_co-occupancy.bedpe'))
    

    }
  }





