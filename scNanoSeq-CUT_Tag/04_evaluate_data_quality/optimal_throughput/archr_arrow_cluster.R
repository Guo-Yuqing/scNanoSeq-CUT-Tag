#!/bin/env Rscripts
args <- commandArgs(T)

prefix <- args[1]
rep <- args[2]
cell_num <- args[3]

###################################################################################
library(GenomeInfoDb)
library(ArchR)
library(Rsamtools)
library(reshape2) 
library(tidyr)
addArchRGenome("hg38")
###################################################################################
input_name <- paste0(prefix,"_6celllines.Rep",rep,"_readsnum",cell_num,"K.bed.gz")
inputFiles <- scanTabix(input_name)
names(inputFiles) <- paste0("LOOP_",prefix)

ArrowFiles <- createArrowFiles(
  inputFiles = input_name,
  sampleNames = names(inputFiles),
  minTSS = 1, 
  addTileMat = TRUE,TileMatParams=list(tileSize = 5000,blacklist = NULL),
  addGeneScoreMat = TRUE ,GeneScoreMatParams = list(blacklist = NULL)
)


dir.create("out_dir")
proj1 <- ArchRProject(
  ArrowFiles = paste0("LOOP_",prefix,".arrow"),
  outputDirectory = "out_dir",
  copyArrows = TRUE
)

tmp_name <- paste0("LOOP_",prefix,"#")
tmp <- gsub(tmp_name,'',proj1$cellNames)
celltype <- sapply(strsplit(tmp,'_'),'[',1)
proj1$celltype <- celltype


p1 <- plotGroups(
    ArchRProj = proj1,
    groupBy = "celltype",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "ridges",baseSize = 10
   )

p1_name <- paste0(prefix,"_6celllines.Rep",rep,"_readsnum",cell_num,"_TSSEnrichment.pdf")
plotPDF(p1, name = p1_name, ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

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

color <- c("#E41A1C","#0072B5FF",
          "#4DAF4A",
          "#B45DC2",
          "#FF7F00",
          "#FF3DA8")
names(color) <- proj1@cellColData$celltype %>% unique %>%
  as.character %>% gtools::mixedsort(.)



p1 <-  plotEmbedding(ArchRProj = proj1,
                    colorBy = "cellColData",
                    name = "celltype",
                    embedding = "UMAP",
                    size = 2,pal=color)

p2<- plotEmbedding(ArchRProj = proj1,
                    colorBy = "cellColData",
                    name = "Clusters",
                    embedding = "UMAP",
                    size = 2)


p3_name <- paste0(prefix,"_6celllines.Rep",rep,"_readsnum",cell_num,"_umap.pdf")
plotPDF(p1,p2, name = p3_name, ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)


saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-Proj1", load = FALSE)


cluster_data <- as.data.frame(table(proj1$Clusters,proj1$celltype))
data <- spread(cluster_data, Var2,Freq)
xls_name <- paste0(prefix,"_6celllines.Rep",rep,"_readsnum",cell_num,"_clustering_result.xls")
write.table(data,xls_name,quote = F,sep = '\t',row.names = F)

