#!/bin/bash
workdir=`pwd`
mkdir -p $workdir/L1HS_longreads
cd $workdir/raw_file

bin_file=$workdir/L1HS/blast_result/scNanoCuttag/QUERY_index_Bin_300bp

ls *.bed|while read file
do

sample=$(echo $file|sed 's/_merge_all_cell.rmdup_fragments.site.filter.bed//g')

/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/guoyuqing/software/liftover/liftOver ${sample}_merge_all_cell.rmdup_fragments.site.filter.bed  \
/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/guoyuqing/reference/T2T/LIHS_liftover/grch38-chm13v2.chain  ${sample}_merge_all_cell.rmdup_fragments.site.hg38toT2T.bed   unmap

less ${sample}_merge_all_cell.rmdup_fragments.site.hg38toT2T.bed|  \
bedtools intersect -a  $bin_file -b - -wo | \
awk '{print $0"\t"$8-$7}'|awk '{if ($11>=150)print}'| \
awk -vOFS='\t'  '{print $1,$2,$3,$4,$5,1}'|bedtools groupby -g 1,2,3,4,5 -c 6 -o sum|\
bedtools intersect -a $bin_file -b - -loj|\
awk -vOFS='\t'  '{if ($7>0)print $1,$2,$3,$4,$5,$11;else print $1,$2,$3,$4,$5,0}'| \
sort -k 4,4V -k 5,5V |uniq >$workdir/L1HS_longreads/${sample}_L1HS_longreads.bed

done

