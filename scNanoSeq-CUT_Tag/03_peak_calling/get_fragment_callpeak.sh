#!/bin/bash
cell=$1
antibody=$2
prefix=${cell}_${antibody}

rootdir=~

flank_size=100
chr_size=$rootdir/reference/minimap2/genome/Homo_sapiens.GRCh38.dna.chr.size.txt
cellname=celllines_filter_TSS1_cellname.txt

sampledir=$rootdir/project/$antibody/filter_fragments
mkdir -p $sampledir/$cell
cd $sampledir/$cell

do_flank_site(){
less $cellname|grep $cell|  \
while read id ;do (less $sampledir/../raw_fragments/$cell/${cell}_${antibody}_merge_all_cell.rmdup_fragments.site.bed| \
grep $id >>${cell}_${antibody}_merge_all_cell.rmdup_fragments.site.tmp.bed);done

less ${cell}_${antibody}_merge_all_cell.rmdup_fragments.site.tmp.bed|sort -k 1,1V -k 2,2n -k 3,3n >${cell}_${antibody}_merge_all_cell.rmdup_fragments.site.filter.bed

less ${cell}_${antibody}_merge_all_cell.rmdup_fragments.site.filter.bed| \
        grep -v '\.'| \
        awk -vOFS='\t' \
            '{
            chr=$1
            st=$2
            ed=$3
            st=st+ext/2
            ed=ed-ext/2
            if (st<ed){
                print chr,st,ed
            }
            }' ext=$flank_size| \
        bedtools flank -i - \
        -b $flank_size  \
        -g  $chr_size |sort -k 1,1V -k 2,2n -k 3,3n >${cell}_${antibody}_merge_all_cell.short.flank${flank_size}.bed

}

workdir=$sampledir/$cell
cd $workdir

do_threshold(){
flank_size=100
threshold_size=0.05
proportion=0.015

peak_dir=$workdir/SEACR/flank${flank_size}_threshold${threshold_size}
mkdir -p $peak_dir

}


do_call_SEACR_peak(){

peak_dir=$workdir/SEACR/flank${flank_size}_threshold${threshold_size}
cd $peak_dir
mkdir log

call_peak_file=$workdir/${prefix}_merge_all_cell.short.flank${flank_size}.bed
SEACR_scripts=$rootdir/SEACR-master/SEACR_1.3.sh
chr_size=$rootdir/reference/minimap2/genome/Homo_sapiens.GRCh38.dna.chr.size.txt

bedtools genomecov -bg -i $call_peak_file -g $chr_size > ${prefix}_merge_all_cell.short.flank${flank_size}.bedgraph
bdg_file=${prefix}_merge_all_cell.short.flank${flank_size}.bedgraph
bash  $SEACR_scripts  $bdg_file  ${threshold_size} non stringent ${prefix}_merge_all_cell.sc_seacr_top${threshold_size}.peaks

}

do_cellname(){

      allcell=$(less $cellname|grep $prefix|wc -l)
      cell_pert=$(printf "%.0f" `echo "scale=4;${proportion}*${allcell}"|bc`)
      add=`awk -v num1=$cell_pert -v num2=${allcell} 'BEGIN{print(num1<num2)?"1":"0"}'`

      cellnumber_tmp=$((${cell_pert} + ${add}))

      if [ $cellnumber_tmp -lt 5 ]
      then
        cell_num=5
      elif  [ $cellnumber_tmp -gt 15 ]
      then
        cell_num=$cellnumber_tmp
      else
        cell_num=$cellnumber_tmp
      fi
     echo $cell_num


}

do_filter_peak(){

peak_dir=$workdir/SEACR/flank${flank_size}_threshold${threshold_size}
cd $peak_dir

reads_site_file=$workdir/${prefix}_merge_all_cell.rmdup_fragments.site.filter.bed
peak_file=${prefix}_merge_all_cell.sc_seacr_top${threshold_size}.peaks.stringent.bed

less $reads_site_file|awk -vOFS='\t' '{print $1,$3,$3,$4"_end",$5}' >${prefix}_all_cells_site.end.tmp
less $reads_site_file|awk -vOFS='\t' '{print $1,$2,$2,$4"_start",$5}' >${prefix}_all_cells_site.start.tmp
cat ${prefix}_all_cells_site.start.tmp ${prefix}_all_cells_site.end.tmp|sort -k 1,1V -k 2,2n -k 3,3n >${prefix}_all_cells_site.start.end.sort

bedtools intersect -a ${prefix}_all_cells_site.start.end.sort \
-b $peak_file -wa -wb|  \
awk -vOFS='\t'  '{print $6,$7,$8,$9,$10,$11,$4,$5}'|\
sort -k 1,1V -k 2,2n -k 3,3n -k 7,7V >${prefix}_Peak_readsSite_sort
less ${prefix}_Peak_readsSite_sort|uniq| \
    bedtools groupby -g 1,2,3,4,5,6 -c 8 -o sum| \
    awk '{print $0"\t""'$cell_num'"}'|awk -vOFS='\t' '{if ($7>=$8)print $1,$2,$3,$4,$5,$6,$7 }'| \
    uniq >${prefix}_sc_seacr_top${threshold_size}.peaks.stringent.${cell_num}cell_support.tmp.bed

less ${prefix}_sc_seacr_top${threshold_size}.peaks.stringent.${cell_num}cell_support.tmp.bed | \
sort -k 1,1V -k 2,2n -k 3,3n|\
awk '$0=NR"\t"$0' - | \
awk -vOFS='\t' '{print $2,$3,$4,"SEACRpeak_"$1,$5,$6}'>${prefix}_sc_seacr_top${threshold_size}.peaks.stringent.${cell_num}cell_support.bed

chmod +x ${prefix}_sc_seacr_top${threshold_size}.peaks.stringent.${cell_num}cell_support.bed

rm ${prefix}_all_cells_site.end.tmp ${prefix}_all_cells_site.start.tmp ${prefix}_sc_seacr_top${threshold_size}.peaks.stringent.${cell_num}cell_support.tmp.bed


}


run(){
    do_flank_site
    do_threshold
    do_call_SEACR_peak
    do_cellname
    do_filter_peak
}
run


