#!/bin/bash
cell=$1
antibody=$2
chr_=$3
prefix=${cell}_${antibody}

root_dir=~
workdir=$root_dir/project/co_occupancy
mkdir -p $workdir/$cell
cd $workdir/$cell

do_get_flank(){

mkdir flank_frag
frag_file=$root_dir/project/filter_fragments/${cell}/${prefix}_merge_all_cell.rmdup_fragments.site.filter.bed

  less $frag_file |grep -w $chr_ >${prefix}_${chr_}
  awk '$0=NR"\t"$0'  ${prefix}_${chr_}|awk -vOFS='\t' '{print $2,$3-1,$3,"Reads_"$1,60,"-"}' >${prefix}_${chr_}_tmp1
  awk '$0=NR"\t"$0'  ${prefix}_${chr_}|awk -vOFS='\t' '{print $2,$4-1,$4,"Reads_"$1,60,"+"}' >>${prefix}_${chr_}_tmp1
  sort -k 1,1V -k 2,2n -k 3,3n -k 4,4V ${prefix}_${chr_}_tmp1|bgzip >flank_frag/${prefix}_${chr_}_sort.gz
  rm ${prefix}_${chr_}_tmp1  ${prefix}_${chr_}

}

do_get_peakpair(){
mkdir peakpair
cd peakpair

frag_file=$root_dir/project/filter_fragments/${cell}/${prefix}_merge_all_cell.rmdup_fragments.site.filter.bed
peak_file=$root_dir/project/filter_fragments/$cell/SEACR/flank100_threshold0.05/${prefix}_sc_seacr_top0.05.peaks.stringent.*cell_support.bed

less $frag_file|grep -w $chr_ |awk '$0=NR"\t"$0'|awk -vOFS='\t' '{print $2,$3,$4,"Reads_"$1,$5}' >${prefix}_rmdup_fragments.site.${chr_}
less $peak_file|grep -w $chr_ >${prefix}_peak.${chr_}

bedtools intersect -a ${prefix}_peak.${chr_} -b ${prefix}_rmdup_fragments.site.${chr_} -wo| \
awk '{print $0"\t"$3-$2}'|awk '{if ($12<$13)print}'|sort -k 10,10V -k 4,4V|cut -f 1,2,3,4,5,6,7,8,9,10| \
awk '{print $0"\t"1}'|bedtools groupby -g 10 -c 11,4 -o sum,collapse|awk '{if ($2>1)print}'| \
cut -f3|sed 's/,/\t/g'|sort -k 1,1V -k 2,2V|uniq >${chr_}_peakpairs


less ${chr_}_peakpairs|while read line 
do 
    echo "$line">./${chr_}_${date}_temp.txt
    peak1=`awk '{print $1}' ./${chr_}_${date}_temp.txt`
    peak2=`awk '{print $2}' ./${chr_}_${date}_temp.txt`

    grep -w $peak1 ${prefix}_peak.${chr_}>${chr_}_tmp1
    grep -w $peak2 ${prefix}_peak.${chr_}>${chr_}_tmp2

    paste -d '\t' ${chr_}_tmp1 ${chr_}_tmp2 >>${chr_}_peakpairs.bed

done 

rm ${chr_}_tmp1 ${chr_}_tmp2 ./${chr_}_${date}_temp.txt
rm ${prefix}_rmdup_fragments.site.${chr_}  ${prefix}_peak.${chr_}

}


run(){
	do_get_flank
	do_get_peakpair
}

run



