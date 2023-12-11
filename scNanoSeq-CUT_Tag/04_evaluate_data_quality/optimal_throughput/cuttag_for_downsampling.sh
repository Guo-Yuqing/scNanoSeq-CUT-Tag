#!/bin/bash
prefix=$1
rep=$2
throughput=$3

rootdir=~
workdir=$rootdir/H3K4me3_downsampling

cd $workdir
reads_num=${throughput}000
outdir=$workdir/downsampling/readsnum${throughput}K/Rep${rep}_readsnum${throughput}K
mkdir -p $outdir

cellname=$workdir/downsample_cell_${prefix}

for sample in K562 GM12878 H9 HFF1 293T HG002
do 

    less $cellname|grep ${sample}|while read cell
    do    
        less $workdir/${sample}_${prefix}_merge_all_cell.rmdup_fragments.site.bed| \
		grep $cell| \
		shuf -n $reads_num  - >>$outdir/${prefix}_6cellline.rmdup_fragments.site.Rep${rep}_readsnum${throughput}K.bed    
    done
done

cd $outdir
less ${prefix}_6cellline.rmdup_fragments.site.Rep${rep}_readsnum${throughput}K.bed | \
awk '{if($2-50>0)print $1"\t"$2-50"\t"$2+50"\t"$4"\t"$5;else print $0}' >${prefix}_6cellline.rmdup_fragments.site.Rep${rep}_readsnum${throughput}K.tmp1

less ${prefix}_6cellline.rmdup_fragments.site.Rep${rep}_readsnum${throughput}K.bed | \
awk '{if($3-50>0)print $1"\t"$3-50"\t"$3+50"\t"$4"\t"$5;else print $0}' >${prefix}_6cellline.rmdup_fragments.site.Rep${rep}_readsnum${throughput}K.tmp2

cat *.tmp* |sort -k 1,1V -k 2,2n -k 3,3n -k 4,4V >${prefix}_6celllines.Rep${rep}_readsnum${throughput}K.bed

bgzip=$rootdir/software/miniconda3/bin/bgzip
tabix=$rootdir/software/miniconda3/bin/tabix

$bgzip ${prefix}_6celllines.Rep${rep}_readsnum${throughput}K.bed
$tabix -p bed ${prefix}_6celllines.Rep${rep}_readsnum${throughput}K.bed.gz

rm ${prefix}_6cellline.rmdup_fragments.site.Rep${rep}_readsnum${throughput}K.tmp1 ${prefix}_6cellline.rmdup_fragments.site.Rep${rep}_readsnum${throughput}K.tmp2

archr_scripts=$workdir/archr_arrow_cluster.R
Rscript=$rootdir/software/miniconda3/envs/r403/bin/Rscript 

$Rscript $archr_scripts $prefix $rep $throughput

