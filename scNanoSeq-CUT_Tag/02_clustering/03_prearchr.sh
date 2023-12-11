#!/bin/bash
antibody=$1
flank_size=100
for cell in GM12878 HG002 K562 293T H9 HFF1 
do 
  prefix=${cell}_${antibody}
  cd ${cell}
  
  cat ${cell}_*/${cell}_${antibody}_*.rmdup_fragments|grep -v '\.' | \
  awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5}'| \
  sort -k 1,1V -k 2,2n -k 3,3n >${prefix}_merge_all_cell.rmdup_fragments.site.bed
  
  cd ../
done

mkdir arrow_file  && cd arrow_file
ln -s ../*/*_${antibody}_merge_all_cell.rmdup_fragments.site.bed  .

prefix=${antibody}
for i in K562 GM12878 293T HFF1 H9 HG002
do

  less ${i}_${prefix}_merge_all_cell.rmdup_fragments.site.bed| \
  awk '{if($2-50>0)print $1"\t"$2-50"\t"$2+50"\t"$4"\t"$5;else print $0}' >${i}_${prefix}_merge_all_cell.rmdup_fragments.site.tmp1
  less ${i}_${prefix}_merge_all_cell.rmdup_fragments.site.bed| \
  awk '{if($3-50>0)print $1"\t"$3-50"\t"$3+50"\t"$4"\t"$5;else print $0}' >${i}_${prefix}_merge_all_cell.rmdup_fragments.site.tmp2

done

cat *.tmp* |sort -k 1,1V -k 2,2n -k 3,3n -k 4,4V >${prefix}_6celllines.bed
bgzip ${prefix}_6celllines.bed
tabix -p bed ${prefix}_6celllines.bed.gz


