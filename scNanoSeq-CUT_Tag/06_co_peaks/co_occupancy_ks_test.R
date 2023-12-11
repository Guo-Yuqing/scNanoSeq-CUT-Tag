#!/usr/bin/env Rscripts
###
 # @Descripttion:scNanoSeq-CUT&Tag co-occupancy
 # @version:
 # @LastEditTime: 2023-6-24 22:35:20
###
argv <- commandArgs(TRUE)

cell.type <- argv[1]
modification <- argv[2]
chrom <- argv[3]
cores = 5

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
print("Read table")
# nanocuttag fragments
frags.chr <- fread(paste0("../flank_frag/",cell.type,"_",modification,"_",chrom,"_sort.gz"))
colnames(frags.chr) <- c('chr', 'start', 'end', 'read', 'mapQ', 'side')

# nanocuttag peaks
peaks <-
  readBed(paste0(cell.type,"_",modification,"_peak.bed")) 
# get subset of peaks by chromosome
peaks.chr <- subset(peaks,seqnames==chrom)

# get subset of fragments by chromosome and filter supplemntary alignemnts

tb <- table(frags.chr$read)
frags.chr <- frags.chr %>% subset(! read %in% names(tb)[tb>2])
frags.chr <- frags.chr[order(frags.chr$read,frags.chr$chr,frags.chr$start),]
rm(tb);rm(frags);gc()

# calculate read length
rl.chr <- subset(frags.chr, side == '+')$start - 
  subset(frags.chr, side == '-')$start
names(rl.chr) <- frags.chr$read %>% unique

# the worker function for co-occupancy hypothesis test
print("KS_test")
worker <- function(i) {
  print(i)
  peak <- resize(peaks.chr[i], width = 1e3, fix = 'center')
  control <- resize(peaks.chr[i], width = 1e5, fix = 'center')
  fr <- subset(frags.chr, start >= start(peak) & end <= end(peak))
  
  fr.control <-
    subset(frags.chr, start >= start(control) & end <= end(control))
  result <- list()
  for (sd in c('-', '+')) {
    index <- subset(fr, side == sd)$read
    index.control <- subset(fr.control, side == sd)$read
    if (length(index) != 0) {
      rl.chr.peak <- rl.chr[index]
      rl.chr.peak <- rl.chr.peak[rl.chr.peak>1000]
      rl.chr.control <- rl.chr[index.control]
      rl.chr.control <- rl.chr.control[rl.chr.control>1000]
      if (length(rl.chr.peak)>0){
      pval <- ks.test(rl.chr.peak, rl.chr.control)$p.value
      result[[sd]] <-
        c(pval, median(rl.chr.peak), median(rl.chr.control))
      } else{
        result[[sd]] <- c(NA, NA, NA) 
    }
      }else{
      result[[sd]] <- c(NA, NA, NA)
    }
  }
  result
}


result <- mclapply(seq_along(peaks.chr), worker, mc.cores = cores)

# foramt hypothesis test results
peak.test <- data.frame(
  chr = seqnames(peaks.chr),
  start = start(peaks.chr),
  end = end(peaks.chr),
  left.p.val = sapply(result, function(x) x[['-']][1]),
  left.length = sapply(result, function(x) x[['-']][2]),
  left.control= sapply(result, function(x) x[['-']][3]),
  right.p.val = sapply(result, function(x) x[['+']][1]),
  right.length = sapply(result, function(x) x[['+']][2]),
  right.control = sapply(result, function(x) x[['+']][3])
) %>%
  mutate(left.p.val.adj=p.adjust(left.p.val,method = 'fdr'),
         right.p.val.adj=p.adjust(right.p.val,method = 'fdr'))

saveRDS(peak.test, paste0('scNanoCuttag_',cell.type,"_",modification,'_',chrom,'.co-occupancy.ks.Rds')) 



