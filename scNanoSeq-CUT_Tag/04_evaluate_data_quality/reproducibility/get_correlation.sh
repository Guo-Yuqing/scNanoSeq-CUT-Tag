#!bin/bash
cellline=$1

rootdir=~
multiBigwigSummary=$rootdir/software/miniconda3/envs/DNAloop/bin/multiBigwigSummary
plotCorrelation=$rootdir/software/miniconda3/envs/DNAloop/bin/plotCorrelation
 
$multiBigwigSummary bins  -b   \
	*.bw  \
	--binSize 5000  \
    --numberOfProcessors 10 \
    --outRawCounts ${cellline}_5000.txt \
    -o ${cellline}_5000.npz \


$plotCorrelation \
    -in ${cellline}_5000.npz   --removeOutliers \
    --corMethod pearson \
    --skipZeros  --zMax 1  \
    --plotTitle "Pearson Correlation of Read Counts" \
    --whatToPlot heatmap \
    --colorMap RdBu_r \
    -o ${cellline}_5000_heatmap_pearsonCorr.pdf 




