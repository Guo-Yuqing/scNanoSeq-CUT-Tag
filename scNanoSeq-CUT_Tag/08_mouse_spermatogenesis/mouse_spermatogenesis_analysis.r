#!/usr/bin/env Rscripts
###
 # @Descripttion:
 # @version:
 # @Author: Yuqing Guo
 # @Date: 2023-5
 # @LastEditors: Yuqing Guo
 # @LastEditTime: 2023-10
###

setwd('mouse_Testis/raw_fragments')
library(GenomeInfoDb)
library(ArchR)
addArchRGenome("mm10")
library(Rsamtools)
library(ggsci)
library(colorspace)
################################################################################################################################################
inputFiles <- scanTabix("B6D2F1_Testis_H3K4me3.bed.gz")
names(inputFiles) <- "B6D2F1_H3K4me3_Testis"
ArrowFiles <- createArrowFiles(
  inputFiles = "B6D2F1_Testis_H3K4me3.bed.gz",
  sampleNames = names(inputFiles),
  minTSS = 1,
  minFrags = 2000,
  addTileMat = TRUE,TileMatParams=list(tileSize = 5000,blacklist = NULL),
  addGeneScoreMat = TRUE ,GeneScoreMatParams = list(blacklist = NULL)
)

dir.create("out_dir")
proj1 <- ArchRProject(
  ArrowFiles = "./B6D2F1_H3K4me3_Testis.arrow",
  outputDirectory = "out_dir",
  copyArrows = TRUE 
)
################################################################################################################################################
tmp <- gsub('B6D2F1_H3K4me3_Testis#','',proj1$cellNames)
write.table(tmp,file = 'B6D2F1_H3K4me3_Testis_filter_cellname.txt',sep = '\t',quote = F,row.names = F,col.names = F)

batch <- as.data.frame(sapply(strsplit(tmp,'_'),'[',4))
colnames(batch) <- 'batch'
batch$diff <- 'Batch2'
batch[batch$batch=='T14',]$diff <- 'Batch1'
table(batch$diff)
proj1$batch <- batch$diff
proj1$Mouse <- sapply(strsplit(tmp,'_'),'[',1)

tmp2 <- as.data.frame(tmp)
tmp2$batch <- sapply(strsplit(tmp,'_'),'[',4)
tmp2$barcode <- sapply(strsplit(tmp,'_'),'[',5)
tmp2$mouse <- gsub('B6D2F1','',sapply(strsplit(tmp,'_'),'[',1))
tmp2[tmp2$batch=='T14' & tmp2$barcode %in% c( "Bc77", "Bc79" ,"Bc80" ,"Bc81" ,"Bc82" ,"Bc83" ,
                                                "Bc86" ,"Bc87", "Bc88" ,"Bc89"),"mouse"] <- 'M1'

tmp2[tmp2$batch=='T14' & tmp2$barcode %in% c( "Bc90", "Bc91" ,"Bc92" ,"Bc93" ,"Bc94" ,
                                             "Bc95" ,"Bc96"),"mouse"] <- 'M2'
proj1$mouse <- tmp2$mouse
################################################################################################################################################

proj1 <- addIterativeLSI(
    ArchRProj = proj1,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 2,
    dimsToUse = 1:30,force = T
)

proj1 <- addClusters(
    input = proj1,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    force = T
)

proj1 <- addUMAP(
    ArchRProj = proj1,
    reducedDims = "IterativeLSI",
    name = "UMAP", force = T,verbose =F
)

################################################################################################################################################
p1 <- plotEmbedding(ArchRProj = proj1,
                    colorBy = "cellColData",
                    name = "mouse",
                    embedding = "UMAP",
                    size = 1.5)
p2 <- plotEmbedding(ArchRProj = proj1,
                    colorBy = "cellColData",
                    name = "Clusters",
                    embedding = "UMAP",
                    size = 1.5)
p3 <- plotEmbedding(ArchRProj = proj1,
                    colorBy = "cellColData",
                    name = "batch",
                    embedding = "UMAP",
                    size = 1.5)

plotPDF(p1,p2,p3 name = "B6D2F1_H3K4me3_Testis_UMAP.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

################################################################################################################################################
proj1 <- addImputeWeights(proj1)

p <- plotEmbedding(
  ArchRProj = proj1, size = 1.5,
  colorBy = "GeneScoreMatrix", 
  name = c("DMRT1","NANOS3",'KIT','DAZL',#SPG
            'RAD51','DICER1','H2AFZ','SYCE1','MEIOB','SYCP3',#LZ
            'CETN1','POU5F2','ADAM3','FA2H','OVOL2','NME8',#PD
            "TDRD1",'TXNDC2','TNP1',"PRM1","KLF17", 'DYRK4', #Sperm
            'SOX9',"WT1"#sertoli
          ), 
  embedding = "UMAP",plotAs = "points",
  imputeWeights = getImputeWeights(proj1)
)

plotPDF(plotList = p, name = "B6D2F1_H3K4me3_testis_marker.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

################################################################################################################################################

##identify celltype according to gene score

################################################################################################################################################

color <- c('red','#29A7E1',
'#0E6EB8','#122B88','#E5007F','#FF7F0EFF','#FFD700',
'#008B45FF','grey'
)
names(color) <- proj1@cellColData$celltype %>% unique %>% 
  as.character %>% gtools::mixedsort(.)

p1 <- plotEmbedding(ArchRProj = proj1,
                    colorBy = "cellColData",
                    name = "celltype",
                    embedding = "UMAP",
                    size = 1.5,pal=color,baseSize = 10)

p2 <- plotEmbedding(ArchRProj = proj1,
                    colorBy = "cellColData",
                    name = "Clusters",
                    embedding = "UMAP",
                    size =1.5,baseSize = 10)
p3 <- plotEmbedding(ArchRProj = proj1,
                    colorBy = "cellColData",
                    name = "batch",
                    embedding = "UMAP",
                    size =1.5,baseSize = 10)

p4 <- plotEmbedding(ArchRProj = proj1,
                    colorBy = "cellColData",
                    name = "mouse",
                    embedding = "UMAP",
                    size = 1.5,baseSize = 10)

plotPDF(p1,p2,p3,p4, name = "B6D2F1_H3K4me3_testis_UMAP.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

################################################################################################################################################
trajectory <- c("1-SPG","2-LZ","3-PD", "4-SPC","5-Sperm1","6-Sperm2","7-Sperm3")

proj1 <- addTrajectory(
    ArchRProj = proj1, 
    name = "Sperm", 
    groupBy = "celltype",
    trajectory = trajectory, 
    embedding = "UMAP", 
    force = TRUE
)

gene_to_plot <- c("DMRT1","NANOS3",'KIT','DAZL',#SPG
            'RAD51','DICER1','H2AFZ','SYCE1','MEIOB','SYCP3',#LZ
            'CETN1','POU5F2','ADAM3','FA2H','OVOL2','NME8',#PD
            "TDRD1",'TXNDC2','TNP1',"PRM1","KLF17", 'DYRK4', #Sperm
            'SOX9',"WT1"#sertoli
            )

plot_list  <-  list() 
for (i in 1:length(gene_to_plot)) { 
  p = plotTrajectory(proj1, 
                     trajectory = "Sperm", 
                     colorBy = "GeneScoreMatrix", 
                     name = gene_to_plot[i], continuousSet = "zissou",
                     size = 1.5,quantHex =2,plotAs = 'points',baseSize = 10)
  plot_list[[i]] = p 
} 


################################################################################################################################################
getGroupBW(
    ArchRProj = proj1,
    groupBy = 'celltype',
    normMethod = "nFrags",
    tileSize = 100,
    maxCells = 1e4
  )
################################################################################################################################################
idxPass <- which(!(proj1$celltype %in% c("Sertoli" ,"Undefined")) )
cellsPass <- proj1$cellNames[idxPass]
proj_germline<- proj1[cellsPass, ]

markersGS <- getMarkerFeatures(
    ArchRProj = proj_germline, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "celltype2",
    useGroups = c( "1-SPG","2-LZ","3-PD","4-SPC" ,"5-Sperm1" ,"6-Sperm2","7-Sperm3"),
    bgdGroups = c( "1-SPG","2-LZ","3-PD","4-SPC" ,"5-Sperm1" ,"6-Sperm2","7-Sperm3"),
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")
write.table(markerList$`1-SPG`,"H3K4me3_SPG_diffgene_FDR_0.05_Log2FC_1.xls",sep = '\t',col.names = T,row.names = F,quote = F)
write.table(markerList$`2-LZ`,"H3K4me3_LZ_diffgene_FDR_0.05_Log2FC_1.xls",sep = '\t',col.names = T,row.names = F,quote = F)
write.table(markerList$`3-PD`,"H3K4me3_PD_diffgene_FDR_0.05_Log2FC_1.xls",sep = '\t',col.names = T,row.names = F,quote = F)
write.table(markerList$`4-SPC`,"H3K4me3_SPC_diffgene_FDR_0.05_Log2FC_1.xls",sep = '\t',col.names = T,row.names = F,quote = F)
write.table(markerList$`5-Sperm1`,"H3K4me3_Sperm1_diffgene_FDR_0.05_Log2FC_1.xls",sep = '\t',col.names = T,row.names = F,quote = F)
write.table(markerList$`6-Sperm2`,"H3K4me3_Sperm2_diffgene_FDR_0.05_Log2FC_1.xls",sep = '\t',col.names = T,row.names = F,quote = F)
write.table(markerList$`7-Sperm3`,"H3K4me3_Sperm3_diffgene_FDR_0.05_Log2FC_1.xls",sep = '\t',col.names = T,row.names = F,quote = F)

saveRDS(markersGS,'testis_markersGS.rds')
################################################################################################################################################

repeat_file <- read.table('/date/guoyuqing/scnanoseqcuttag/project/mouse_Testis/peak/repeatA_repeatT_for_archr/trans_rmsk_l1base_mouse.L1Md_A_T_genomic.bed',
                       header = F,sep='\t')
colnames(repeat_file) <- c('Chr','Str','End','strand','id','symbol')

repeat_file_gr <- GRanges(
      seqnames = Rle(repeat_file$Chr),
      ranges = IRanges(repeat_file$Str, end = repeat_file$End),
      strand = Rle(repeat_file$strand),
      id = repeat_file$id,
      symbol=repeat_file$symbol)

proj1 <- addGeneScoreMatrix(
  input = proj,
  genes = repeat_file_gr,
  geneModel = "exp(-abs(x)/5000) + exp(-1)",
  matrixName = "RepeatScoreMatrix",
  extendUpstream = c(1000, 1e+05),
  extendDownstream = c(1000, 1e+05),
  geneUpstream = 5000,
  geneDownstream = 0,
  useGeneBoundaries = TRUE,
  useTSS = FALSE,
  extendTSS = FALSE,
  tileSize = 5000,
  ceiling = 4,
  geneScaleFactor = 5,
  scaleTo = 10000,
  excludeChr = c("chrY", "chrM"),
  blacklist = NULL,
  threads = getArchRThreads(),
  parallelParam = NULL,
  subThreading = TRUE,
  force = TRUE,
  logFile = createLogFile("addRepeatScoreMatrix")
)

proj1 <- addImputeWeights(proj1)
plotEmbedding(
  ArchRProj = proj1, size = 3,
  colorBy = "RepeatScoreMatrix", 
      name = c("chr6_43625720_43631165_Tpk1_L1Md_T"
  ), 
  embedding = "UMAP",plotAs = "points",continuousSet = "comet",
  base=8,imputeWeights = getImputeWeights(proj1)
)


saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-Proj1", load = FALSE)




