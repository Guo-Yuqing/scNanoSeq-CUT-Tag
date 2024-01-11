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
      
      peak_pair <- peak_pair[peak_pair$left_judge<0.01 & peak_pair$right_judge<0.01  & peak_pair$V3-peak_pair$V2<=1000 &peak_pair$V9-peak_pair$V8<=1000 ,]
      peak_pair <- na.omit(peak_pair)


      sig.left.index <- peak_pair$left_index
      sig.right.index <- peak_pair$right_index
      
      # generate a GenomicInteraction object
      gi <- GInteractions(sig.left.index, sig.right.index, peaks.chr)
      
      coa <- data.frame(
        span = pairdist(gi),
        left.end = peak.test[sig.left.index,]$left.length,
        right.end = peak.test[sig.right.index,]$right.length,
        left.control = peak.test[sig.left.index, ]$left.control,
        right.control = peak.test[sig.right.index, ]$right.control
      ) %>%
        mutate(
          left.spring = left.end / left.control,
          right.spring = right.end / right.control,
          res = span / left.end
        )

      if (is.null(gi.clean)){
        gi.clean<-gi[coa$res<2]
      }else{
        gi.clean<-c(gi.clean,gi[coa$res<2])
      }
      
      if (is.null(coa.all)){
        coa.all<-coa
      }else{
        coa.all<-rbind(coa.all,coa)
      }
    }
    
    export.bedpe(gi.clean,
                 paste0("../coa_bedpe/",cell.type,"_",modification,'_co-occupancy.bedpe'))
    
    # density plot of span ~ linkage
    {
      set.seed(1)
      
      df <- data.frame(
        log10.span = log10(coa.all$span),
        log10.left.end = log10(coa.all$left.end),
        log10.res = log10(coa.all$res)
      ) 
      
      df$density <-
        get_density(df$log10.span,
                    df$log10.left.end,
                    n = 100)
      
      pdf(paste0("../coa_bedpe/",cell.type,"_",modification,'.span_relationship.pdf'),
          width = 7, height = 7)
      
      # color by density
      p <- ggplot(df,aes(x=log10.span, y=log10.left.end, color = density)) +
        geom_point() +
        stat_density2d(h = c(0.5,0.05), color='grey90') +
        geom_hline(yintercept = median(log10(coa.all$right.control))) +
        geom_abline(slope = 1,intercept = -log10(2))+
        scale_color_viridis()+
        coord_fixed(ratio=diff(range(df$log10.span))/diff(range(df$log10.left.end)))+
        theme_bw()
      
      plot(p)
      
      df <-
        rbind(
          data.frame(dist = peaks %>% start %>% diff, group = 'nanoCuttag_peak_distance'),
          data.frame(dist = coa.all$left.end, group = 'median_peak_read_length'),
          data.frame(dist = coa.all$span[coa.all$res<2], group = 'co-occupancy_span \n (res < 2-fold)')
        )
      
      p <- ggplot(data = df, aes(log10(dist), fill = group)) +
        geom_density(alpha = .5) +
        geom_vline(xintercept = median(log10(coa.all$left.control)))+
        theme_bw()
      
      plot(p)
      
      dev.off()
    }
  }





