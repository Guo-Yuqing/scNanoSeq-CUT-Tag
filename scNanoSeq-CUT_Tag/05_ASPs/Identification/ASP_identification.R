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


###########################################################################################
###########################################################################################
testHaplotypeAsymmetric <- function(peak,frag,mc=NULL) {
  frag.cut.off=20
  
  worker <- function(i) {
    result <- list()
    isec.frag <- frag[unique(hits@from[hits@to == i])]
    
    for (ctg in unique(isec.frag$contig)) {
      isec.ctg.frag <- isec.frag[isec.frag$contig == ctg]
      count.ps1 <-
        isec.ctg.frag$cell[isec.ctg.frag$phase == '1'] %>%
        unique %>%
        length
      count.ps2 <-
        isec.ctg.frag$cell[isec.ctg.frag$phase == '2'] %>%
        unique %>%
        length
      
      if (count.ps1 + count.ps2 >= frag.cut.off){
        binom.p <-
          binom.test(count.ps1, count.ps1 + count.ps2, p = 0.5)$p.value
        result$peak <- c(result$peak, i)
        result$contig <- c(result$contig, ctg)
        result$p.value <- c(result$p.value, binom.p)
        result$ratio <-
          c(result$ratio, count.ps1 / (count.ps1 + count.ps2))
        result$count <- c(result$count, count.ps1 + count.ps2)
      }
    }
    result
  }
  
  hits <- findOverlaps(frag, peak)
  result <- mclapply(unique(hits@to), worker, mc.cores = mc)
  
  result.df <- data.frame(
    peak = result %>% sapply(function(x)
      x$peak) %>% unlist,
    contig = result %>% sapply(function(x)
      x$contig) %>% unlist,
    p.value = result %>% sapply(function(x)
      x$p.value) %>% unlist,
    ratio = result %>% sapply(function(x)
      x$ratio) %>% unlist,
    count = result %>% sapply(function(x)
      x$count) %>% unlist
  )
  result.df$fdr <- result.df$p.value %>% p.adjust(method = 'fdr')
  
  result.df<-result.df[order(result.df$count, decreasing = T), ]
  result.df$chr <-peak[result.df$peak] %>% seqnames %>% as.character
  result.df
}


loadFrags<-function(frag.file){
  df <-
    fread(frag.file,
          header = F,
          sep = '\t') %>% as.data.frame()
  names(df) <-
    c('chr', 'start', 'end', 'read', 'cell',  'contig', 'phase')
  frag <- GRanges(
    seqnames = df$chr,
    ranges = IRanges(df$start, df$end),
    read = df$read,
    cell = df$cell,
    contig = df$contig,
    phase = df$phase
  )
  seqlevelsStyle(frag) <- 'UCSC'
  frag
}

#####################################################################################################

df <-fread('./GM12878_H3K4me3_sc_seacr_top0.05.peaks.stringent.17cell_support.bed',
        header = F,
        sep = '\t') %>% as.data.frame()
names(df) <-
  c('chr', 'start', 'end','Peak_name','score','signal')
peak <- GRanges(seqnames = df$chr,
                ranges = IRanges(df$start, df$end))
seqlevelsStyle(peak) <- 'UCSC'


frags <-loadFrags('./GM12878_H3K4me3_PATMAT_flank100.gz')
result <- testHaplotypeAsymmetric(peak, frags, mc = 1)

saveRDS(result, 'GM12878_H3K4me3_TGS_peaks_biased_20frag.Rds')

table(result$contig)

#####################################################################################################
#####################################################################################################
tmp <- result.pl[result.pl$fdr <= cut.off,]
tmp$per <- round(tmp$ratio,1)
table(tmp$per)

cut.off <- 5e-2

pdf('newfilter_GM12878_H3K4me3_TGS_peaks_fdr_5e-2_20frag.pdf',width = 5,height = 5)

result[result$ratio<=0.9 & result$ratio>=0.1,]$fdr=1
result.pl<-result[order(result$fdr,decreasing = T),]

# color paternal and maternal specific peaks
col <- c('grey70','deepskyblue2', 'firebrick1')[
  (result.pl$fdr <= cut.off) + (result.pl$fdr <= cut.off & result.pl$ratio < .5)  + 1
] %>% adjustcolor(alpha.f = .5)


plot(
  result.pl$count,
  result.pl$ratio,
  pch = 20,
  col = col,
  main = 'Allele Ratio ~ Fragment Count',
  xlab = 'Haplotyped Reads per Peaks',
  ylab = 'Paternal/(Paternal+Maternal)',
  xlim = c(0, max(result.pl$count))
)

# highlight chrX
# put chrX points on the upper layer
result.pl<-rbind(subset(result.pl,chr!='chrX'),subset(result.pl,chr=='chrX'))
sig.chrx<-(result.pl$fdr <= cut.off & result.pl$chr == 'chrX')
col <- c('grey70', 'deepskyblue2', 'firebrick1')[
  sig.chrx + (sig.chrx & result.pl$ratio < .5) + 1] %>%
  adjustcolor(alpha.f = .5)
cex <- c(1,1.5)[(result.pl$fdr <= cut.off & result.pl$chr=='chrX')+1]

plot(
  result.pl$count,
  result.pl$ratio,
  cex = cex,
  col = col,
  pch = 20,
  xlab = 'Haplotyped Reads per Peaks',
  ylab = 'Paternal/(Paternal+Maternal)',
  xlim = c(0, max(result.pl$count)),
  main = 'Allele Ratio ~ Fragment Count (chrX only)'
)

dev.off()



pdf(
  'newfilter_GM12878_H3K4me3_biased_ratio_by_chromosome_fdr_5e-2_20frag.pdf',
  width = 4,
  height = 4
)
tb <- table(result$chr, result$fdr <= cut.off)
df <- data.frame(biased.peak.ratio = tb[, 2] / rowSums(tb),
                 chr = rownames(tb))
df$chr<-factor(df$chr,paste0('chr',c(1:22,'X')))
df$col<-'black'
df$col[df$chr=='chrX']<-'red'


ggplot(df, aes(x = chr, y = biased.peak.ratio)) +
  geom_bar(stat = "identity", fill = df$col)  +
#  scale_y_continuous(limits=c(0, 0.2), expand = c(0, 0)) +
  xlab('Chromosome') +
  ylab('ASPs/All peaks') +
  ggtitle('GM12878 Biased Peak Distribution') +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = 'black', angle=90),
    axis.text.y = element_text(color = 'black'),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", size = 1)
  )
dev.off()


####################################################################################################################
cut.off <- 5e-2
result.pl<-result[order(result$fdr,decreasing = T),]
tmp <- result.pl[result.pl$fdr <= cut.off ,]
tmp$per <- round(tmp$ratio,1)
table(tmp$per)

tmp$peak_name <- paste0('SEACRpeak_',tmp$peak)

write.table(tmp,'GM12878_TGS_APS_peak_fdr_5e-2_20frag.bed',sep = '\t',quote = F,row.names = F)

