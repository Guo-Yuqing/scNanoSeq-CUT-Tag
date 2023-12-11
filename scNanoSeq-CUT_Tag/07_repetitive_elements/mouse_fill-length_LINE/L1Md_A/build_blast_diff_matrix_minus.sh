#!/bin/bash
workdir=`pwd`
total_blast=$workdir/trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast
motif_scripts=$workdir/motif_ref_query_diff.new.R
Rscript=/datg/guoyuqing/software/miniconda3/envs/r403/bin/Rscript

file1=$workdir/trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast_index_anno
file2=$workdir/QUERY_index

less $file1|grep 'Minus' |grep -v '#'|while read line 
do 

    echo "$line">./${date}_temp.txt
    start_index=`awk '{print $1}' ./${date}_temp.txt`
    end_index=`awk '{print $2}' ./${date}_temp.txt`
    region=`awk '{print $3}' ./${date}_temp.txt|sed 's/Query=//g'`
    rm ./${date}_temp.txt
    region_range=$(($end_index-$start_index))
    query=$(less $file2|grep $region|cut -f7)

    mkdir -p detail/$query
    cd detail/$query
	
    less $total_blast|tail -n  +${start_index} -|head -n $region_range|grep Query -A3|grep -v = |grep -v Score|grep -v '^-' >${query}
    awk '$0=NR"\t"$0'  ${query}|grep Sbjct|awk -vOFS='\t' '{print $1,$2}' >${query}_index

    less ${query}_index|cut -f1|while read index
    do

    end=$index
    str=$(($index-2))

    less ${query}|head -${end}|tail -n +${str}|head -1|awk '{print $3}'|sed 's/./&\n/g'|awk '{print toupper($0)}'|grep -v '^$' >tmp_query
    less ${query}|head -${end}|tail -n +${str}|tail -1|awk '{print $3}'|sed 's/./&\n/g'|awk '{print toupper($0)}'|grep -v '^$' >tmp_ref

    paste -d '\t' tmp_ref tmp_query|awk '{if ($1==$2)print $0"\t"".";else if ($2=="-") print $0"\t""GAP";else if ($1=="-") print $0"\t""INS";else print $0"\t"$2}'  >>ref_query_diff

    done

    start=$(less ${query}|grep -v '^$'|head -3|tail -1|awk '{print $2}')
    diff=$((4961-$start))

    end=$(less ${query}|grep -v '^$'|tail -1|awk '{print $4}')
    diff2=$(($end-1))

    if [ $diff -gt 0 ]
    then

        for (( c=1; c<=diff; c++)) ; do echo -e "X""\t""-""\t""GAP" >>miss_line; done
        cat miss_line ref_query_diff >${query}_ref_query_diff1
    else
        cp ref_query_diff ${query}_ref_query_diff1
    fi


    if [ $diff2 -gt 0 ]
    then

        for (( c=1; c<=diff2; c++)) ; do echo -e "X""\t""-""\t""GAP" >>miss_line2; done
        cat ${query}_ref_query_diff1 miss_line2 >${query}_ref_query_diff_final_tmp
    else
        cp ${query}_ref_query_diff1 ${query}_ref_query_diff_final_tmp
    fi

    tac ${query}_ref_query_diff_final_tmp|sed 's/T/M/g'|sed 's/C/N/g'|sed 's/A/T/g'|sed 's/G/C/g'|sed 's/M/A/g'|sed 's/N/G/g'|sed 's/CTP/GAP/g'|sed 's/IGS/INS/g'>${query}_ref_query_diff_final

    rm miss_line2 ${query}_ref_query_diff1
    rm ${query}_ref_query_diff_final_tmp
    rm  miss_line  ref_query_diff
    rm tmp_query  tmp_ref miss_line

    $Rscript $motif_scripts ${query}

    cd ../../

done 
