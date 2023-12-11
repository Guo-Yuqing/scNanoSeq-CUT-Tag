#!/usr/bin/env Rscripts
###
 # @Descripttion:
 # @version:
 # @Author: Yuqing Guo
 # @Date: 2022-12-03 22:35:11
 # @LastEditors: Yuqing Guo
 # @LastEditTime: 2023-6-27 22:35:20
###
argv <- commandArgs(TRUE)
antibody <- argv[1]


library(GenomeInfoDb)
library(ArchR)
library(Rsamtools)
addArchRGenome("hg38")

inputFiles=scanTabix(paste0(antibody,"_6celllines.bed.gz"))
names(inputFiles)=paste0("LOOP_",antibody)
ArrowFiles <- createArrowFiles(
  inputFiles = paste0(antibody,"_6celllines.bed.gz"),
  sampleNames = names(inputFiles),
  minTSS = 0.1, 
  minFrags = 2000,
  addTileMat = TRUE,TileMatParams=list(tileSize = 5000,blacklist = NULL),
  addGeneScoreMat = TRUE ,GeneScoreMatParams = list(blacklist = NULL)
)

dir.create("out_dir")
proj <- ArchRProject(
  ArrowFiles = paste0("LOOP_",antibody,".arrow"),
  outputDirectory = "out_dir",
  copyArrows = TRUE 
)
################################################################################################################################################
###if the modification are repressive markers:H3K27me3;H3K9me3

tmp=gsub(paste0('LOOP_',antibody,'#'),'',proj$cellNames)
celltype=sapply(strsplit(tmp,'_'),'[',1)
proj$celltype=celltype

write.table(tmp,file = 'celllines_filter_TSS0.1_cellname.txt',sep = '\t',quote = F,row.names = F,col.names = F)
proj1=proj

################################################################################################################################################
###if the modification are activate markers:H3K4me3;H3K27ac;H3K36me3;CTCF;RAD21

 idxPass <- which(proj$TSSEnrichment>=1)
 cellsPass <- proj$cellNames[idxPass]
 proj1<- proj[cellsPass, ]

tmp=gsub(paste0('LOOP_',antibody,'#'),'',proj1$cellNames)
celltype=sapply(strsplit(tmp,'_'),'[',1)
proj1$celltype=celltype

write.table(tmp,file = 'celllines_filter_TSS1_cellname.txt',sep = '\t',quote = F,row.names = F,col.names = F)

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
    maxClusters = 6,force = T
)

proj1 <- addUMAP(
    ArchRProj = proj1, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", dimsToUse=1:6,
    force = T,verbose =F
)

################################################################################################################################################
color = c("#E41A1C","#0072B5FF",
          "#4DAF4A",
          "#B45DC2",
          "#FF7F00",
          "#FF3DA8")
names(color) <- proj1@cellColData$celltype %>% unique %>% 
  as.character %>% gtools::mixedsort(.)

p1 <- plotEmbedding(ArchRProj = proj1, 
                    colorBy = "cellColData",
                    name = "celltype", 
                    embedding = "UMAP",
                   size = 2)

p2<- plotEmbedding(ArchRProj = proj1,
                    colorBy = "cellColData",
                    name = "Clusters",
                    embedding = "UMAP",
                   size = 2)

plotPDF(p1,p2, name = "6celllines_UMAP.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-Proj1", load = FALSE)  
################################################################################################################################################

proj1 <- addImputeWeights(proj1)

p=plotEmbedding(
  ArchRProj = proj1, size = 3,
  colorBy = "GeneScoreMatrix", 
  name = c( "GAPDH","BCR","GATA1","PRAME","HBG1","BGLT3","HBG2", #K562  
            "CD19","MS4A1","CD79A","CD79B","HLA-DRB1", #GM12878 & HG002
            "PTPRC","CD3D", "CD3E","CD4","CD8A","CD8B", #H9
            "HOXB5","HOXA5",#293T
            "THY1","COL5A1","COL12A1","DCN" #HFF1
  ), 
  embedding = "UMAP",plotAs = "points",
  imputeWeights = getImputeWeights(proj1)
)

plotPDF(plotList = p, name = "6celllines_marker.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

################################################################################################################################################
getGroupBW(
    ArchRProj = proj1,
    groupBy = "celltype",
    normMethod = "nFrags",
    tileSize = 100,
    maxCells = 1e4
  )
proj1$subtype <- paste0(proj1$celltype,"_",proj1$batch)
getGroupBW(
    ArchRProj = proj1,
    groupBy = 'subtype',
    normMethod = "nFrags",
    tileSize = 100,
    maxCells = 1e4
  )
################################################################################################################################################
markersGS <- getMarkerFeatures(
    ArchRProj = proj1,
    useMatrix = "GeneScoreMatrix",
    groupBy = "celltype",
    useGroups = c( "293T","GM12878","HG002", "H9" ,"HFF1", "K562"),
    bgdGroups = c(  "293T","GM12878","HG002", "H9" ,"HFF1", "K562"),
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")
write.table(markerList$GM12878,"GM12878_diffgene_FDR_0.05_Log2FC_1.xls",sep = '\t',col.names = T,row.names = F,quote = F)
write.table(markerList$HG002,"HG002_diffgene_FDR_0.05_Log2FC_1.xls",sep = '\t',col.names = T,row.names = F,quote = F)
write.table(markerList$H9,"H9_diffgene_FDR_0.05_Log2FC_1.xls",sep = '\t',col.names = T,row.names = F,quote = F)
write.table(markerList$HFF1,"HFF1_diffgene_FDR_0.05_Log2FC_1.xls",sep = '\t',col.names = T,row.names = F,quote = F)
write.table(markerList$`293T`,"293T_diffgene_FDR_0.05_Log2FC_1.xls",sep = '\t',col.names = T,row.names = F,quote = F)
write.table(markerList$K562,"K562_diffgene_FDR_0.05_Log2FC_1.xls",sep = '\t',col.names = T,row.names = F,quote = F)

saveRDS(markersGS,'marker_gene.rds')



