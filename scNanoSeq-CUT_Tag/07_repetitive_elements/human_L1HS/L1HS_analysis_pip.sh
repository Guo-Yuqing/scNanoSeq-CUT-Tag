#!/bin/bash
root_dir=~
T2T_ReferenceL1HS_T2T=T2T_fulllengthL1HS.trans.chrm13v2.0.txt
T2T_fa=$root_dir/reference/T2T/chm13v2.0.fa
T2T_repeat_ref=$root_dir/reference/T2T/T2T_repeat.bed

workdir=$root_dir/L1HS
fmt_downdir=$workdir/format/fmt_down_analysis

active_chrX_bed=$workdir/T2T_fulllengthL1HS.trans.chrm13v2.0.activate_chrX.bed
all_L1_bed=$workdir/T2T_fulllengthL1HS.trans.chrm13v2.0.txt
active_chrX_fa=$workdir/T2T_fulllengthL1HS.trans.chrm13v2.0.activate_chrX.fasta
all_L1_fa=$workdir/T2T_fulllengthL1HS.trans.chrm13v2.0.fasta


do_get_fasta(){

	cd $workdir
	
	less $T2T_repeat_ref|grep LINE| \
	sort -k 1,1V -k 2,2n -k 3,3n |awk '{print $0"\t"$3-$2}' >T2T_repeat.LINE.sort.bed
	
	less T2T_repeat.LINE.sort.bed|grep L1HS>T2T_repeat.LINE.L1HS.sort.bed
	
	less T2T_fulllengthL1HS.chrm13v1.0.txt|grep -v CHROM|grep chr|cut -f 1,2,3|sort -k 1,1V -k 2,2n | \
	awk '{print $1"\t"$2-1000"\t"$3+1000"\t"$0}'| \
	bedtools intersect -a - -b T2T_repeat.LINE.L1HS.sort.bed -wo| \
	grep L1HS|awk '{if ($18>5500)print}'|cut -f 7,8,9,10,11,12,13,14,15,16,17 >${all_L1_bed}
	
	less ${all_L1_bed}|grep chrX|head -1 >${active_chrX_bed}
	
	bedtools getfasta -fi  $T2T_fa -bed ${active_chrX_bed} -fo ${active_chrX_fa}
	bedtools getfasta -fi  $T2T_fa -bed ${all_L1_bed} -fo ${all_L1_fa}


}
##############################################################################################################################################################################################################################################################
do_blast(){

    cd $workdir

    makeblastdb=$root_dir/software/miniconda3/bin/makeblastdb
    blastn=$root_dir/software/miniconda3/bin/blastn
    blast_database=$workdir/T2T_T2T_ReferenceL1HS.fulllength_active_chrX.fasta_database

    $makeblastdb -in ${active_chrX_fa} -dbtype nucl  -parse_seqids  -out ${blast_database}

    diff_format=$workdir/format
    mkdir -p $diff_format


  $blastn  -query ${all_L1_fa}  \
    -out $diff_format/T2T_fulllength_L1HS.trans.chrm13v2.0.blast_info   -db ${blast_database}  -outfmt "7 std stitle"  -max_hsps 1 -xdrop_gap 1000  

   $blastn  -query ${all_L1_fa}  \
    -out $diff_format/T2T_fulllength_L1HS.trans.chrm13v2.0.blast   -db ${blast_database}   -max_hsps 1  -xdrop_gap 1000 
}
##############################################################################################################################################################################################################################################################

do_build_index_new(){
    fmt_downdir=$workdir/format/fmt_down_analysis
    mkdir -p $fmt_downdir
    cd $fmt_downdir

    blast_file=$workdir/format/T2T_fulllength_L1HS.trans.chrm13v2.0.blast
    ln ${blast_file} .

    less  -S  ${active_chrX_fa}|  \
    grep -v ">"|sed 's/./&\n/g'|awk '{print toupper($0)}'|  \
    awk '$0=NR"\t"$0'|awk '{print "Base"$0}'|  \
    sed '1i#Base\tref_active_chrX' >$fmt_downdir/T2T_ReferenceL1HS.fulllength_active_chrX.matrix_file

    less ${all_L1_bed}|awk '$0=NR"\t"$0' -|awk -vOFS='\t'  '{print $2":"$3"-"$4,$5,$6,$7,$8,$9,"Query_"$1}' >$fmt_downdir/QUERY_index

     ln -s $workdir/format/T2T_fulllength_L1HS.trans.chrm13v2.0.blast_info .
    less T2T_fulllength_L1HS.trans.chrm13v2.0.blast_info|grep -v '#'|awk '{if ($4>5000)print}'|cut -f1|sort -k 1,1V|uniq >true_blast_query

    end=$(less T2T_fulllength_L1HS.trans.chrm13v2.0.blast|wc -l)
    awk '$0=NR"\t"$0'  ${blast_file} | \
    grep Query=  >T2T_fulllength_L1HS.trans.chrm13v2.0.blast_index_tmp

    less T2T_fulllength_L1HS.trans.chrm13v2.0.blast_index_tmp|cut -f1| \
    tail -n +2|sed  '$a'$end'' >T2T_fulllength_L1HS.trans.chrm13v2.0.blast_index_tmp_end

    paste -d '\t' T2T_fulllength_L1HS.trans.chrm13v2.0.blast_index_tmp T2T_fulllength_L1HS.trans.chrm13v2.0.blast_index_tmp_end| \
    awk -vOFS='\t' '{print $1,$4,$2,$3}' |sed 's/Query=//g'|  \
    awk 'NR==FNR{a[$1]=1}NR!=FNR{if($3 in a)print $0}' true_blast_query - |\
    awk '{print $1"\t"$2"\t""Query=""\t"$3}'|awk '{if ($2-$1>10)print}' >T2T_fulllength_L1HS.trans.chrm13v2.0.blast_index
    
    less T2T_fulllength_L1HS.trans.chrm13v2.0.blast_index|while read line
    do
    echo "$line">./${date}_temp.txt
    start_index=`awk '{print $1}' ./${date}_temp.txt`
    end_index=`awk '{print $2}' ./${date}_temp.txt`
    query=`awk '{print $3}' ./${date}_temp.txt`
    resion=`awk '{print $4}' ./${date}_temp.txt`
    rm ./${date}_temp.txt

    region_range=$(($end_index-$start_index))
    strand=$(less ${blast_file}|tail -n  +${start_index} -|head -n ${region_range}|grep Strand|sed 's/"Strand=Minus\/"//g')

    echo -e  $start_index"\t"$end_index"\t"$resion"\t"$strand >>T2T_fulllength_L1HS.trans.chrm13v2.0.blast_index_anno
    done   

     rm T2T_fulllength_L1HS.trans.chrm13v2.0.blast_index_tmp  T2T_fulllength_L1HS.trans.chrm13v2.0.blast_index_tmp_end

}


##############################################################################################################################################################################################################################################################
do_get_blast_diff_matrix(){
    cd $fmt_downdir

    plus_scripts=$workdir/build_blast_diff_matrix.plus.sh
    minus_scripts=$workdir/build_blast_diff_matrix_minus.sh
    r_scripts=$workdir/motif_ref_query_diff.new.R

    cp $plus_scripts .
    cp $minus_scripts .
    cp $r_scripts .

    bash $plus_scripts
    bash $minus_scripts

}

##############################################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################################
do_merge(){
    merge_dir=$workdir/format/fmt_down_analysis/detail/merge
    mkdir -p $merge_dir
    cd $merge_dir
    wc -l ../*/*_ref_query_diff_final_modify |grep _modify|awk '{if ($1==6032)print $2}'|while read id ;do (ln -s $id .);done
    ln -s $fmt_downdir/T2T_ReferenceL1HS.fulllength_active_chrX.matrix_file .

    cp T2T_ReferenceL1HS.fulllength_active_chrX.matrix_file T2T_ReferenceL1HS.fulllength_all_blast

    ls *_modify|while read id
    do

    less $id|cut -f3|paste -d '\t' T2T_ReferenceL1HS.fulllength_all_blast - >tmp

    mv tmp T2T_ReferenceL1HS.fulllength_all_blast

    done


}


do_get_Bin_info(){

    bin_region=$workdir/Repeat_bin/Repeat_region.300bin.bed
    outdir=$workdir/format/fmt_down_analysis/detail/merge/repeat_bin/300bp/
    cd $workdir/format/fmt_down_analysis/detail/merge
    mkdir -p $outdir
    ls *_modify|while read id
    do
        name=$(echo $id|sed 's/_ref_query_diff_final_modify//g')
        echo $name

        for bin in "Bin_"{1..20}
        do
            mut=$(less $id|grep -v 'Ref'|awk '$0=NR"\t"$0'| \
            awk -vOFS='\t' '{print "chrRepeat""\t"$1-1"\t"$1"\t"$4}'| \
            bedtools intersect -a - -b  $bin_region  -wo| \
            grep -w ${bin}|awk -vOFS='\t' '{print $1,$2,$3,$4,"Base"$3"_"$4,$6,$7}'|grep -v '\.'|wc -l)

            if  [ $mut = 0 ]
            then

                less $id |grep -v 'Ref'|awk '$0=NR"\t"$0'|  \
                awk -vOFS='\t' '{print "chrRepeat""\t"$1-1"\t"$1"\t"$4}'|  \
                bedtools intersect -a - -b  $bin_region  -wo|grep -w ${bin}| \
                awk -vOFS='\t' '{print $1,$2,$3,$4,"Base"$3"_"$4,$6,$7}'|head -1|awk  -vOFS='\t'  '{print $1,$6,$7,"'$name'","NA"}' >>$outdir/${bin}
            else
                less $id |grep -v 'Ref'|awk '$0=NR"\t"$0'|  \
                awk -vOFS='\t' '{print "chrRepeat""\t"$1-1"\t"$1"\t"$4}'|  \
                bedtools intersect -a - -b  $bin_region  -wo|grep -w ${bin}| \
                awk -vOFS='\t' '{print $1,$2,$3,$4,"Base"$3"_"$4,$6,$7}'|grep -v '\.'| \
                awk '{print $0"\t""'$name'"}'|bedtools groupby -g 1,6,7,8 -c 5 -o collapse >>$outdir/${bin}

            fi
        done
    done

    cd $outdir/
    for bin in "Bin_"{1..20}
    do
    not=$(less ${bin}|sed 's/,/_/g'|sort -k 5,5V|uniq|awk '{print $0"\t"1}'|bedtools groupby -g 1,2,3,5 -c 6 -o sum|sort -k 5,5rn|awk '{if ($5>1)print $5}'|awk '{sum+=$1}END{print sum}')
    yes=$(less ${bin}|sed 's/,/_/g'|sort -k 5,5V|uniq|awk '{print $0"\t"1}'|bedtools groupby -g 1,2,3,5 -c 6 -o sum|sort -k 5,5rn|awk '{if ($5==1||$5==NA)print $5}'|awk '{sum+=$1}END{print sum}')

    echo -e ${bin}"\t"${not}"\t"${yes} >>final_file
    done


}

do_judgement(){
cd $workdir/format/fmt_down_analysis/detail/merge/repeat_bin/300bp/
mkdir judgement

for bin in "Bin_"{1..20}
do

less ${bin}|sed 's/,/_/g'|sort -k 5,5V|  \
uniq|awk '{print $0"\t"1}'|  \
bedtools groupby -g 1,2,3,5 -c 6,4 -o sum,collapse|  \
sort -k 5,5rn|awk '{if ($5==1||$5==NA)print $6"\t""YES"}' >judgement/${bin}_judgement_tmp

less ${bin}|sed 's/,/_/g'|sort -k 5,5V| \
uniq|awk '{print $0"\t"1}'| \
bedtools groupby -g 1,2,3,5 -c 6,4 -o sum,collapse| \
sort -k 5,5rn|awk '{if ($5>1)print $6}'|sed 's/,/\n/g'|awk '{print $1"\t""NO"}' >>judgement/${bin}_judgement_tmp

less judgement/${bin}_judgement_tmp|sort -k 1,1V|sed '1iQuery\t'$bin'_judgement' >judgement/${bin}_judgement

rm judgement/${bin}_judgement_tmp


done

cd judgement

paste -d '\t' Bin_1_judgement Bin_2_judgement Bin_3_judgement Bin_4_judgement Bin_5_judgement Bin_6_judgement Bin_7_judgement Bin_8_judgement Bin_9_judgement Bin_10_judgement Bin_11_judgement Bin_12_judgement Bin_13_judgement Bin_14_judgement Bin_15_judgement Bin_16_judgement Bin_17_judgement Bin_18_judgement Bin_19_judgement Bin_20_judgement|cut -f 1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40|grep -v judgement|sort -k 2,2rV  -k 3,3rV -k 4,4rV -k 5,5rV -k 6,6rV -k 7,7rV -k 8,8rV -k 9,9rV -k 10,10rV -k 11,11rV -k 12,12rV -k 13,13rV -k 14,14rV -k 15,15rV -k 16,16rV -k 17,17rV -k 18,18rV -k 19,19rV -k 20,20rV|sed 's/YES/1/g'|sed 's/NO/0/g'  >merge_judgement.bed


}


run(){
	do_get_fasta
	do_blast
	do_build_index_new
	do_get_blast_diff_matrix
	do_merge
	do_get_Bin_info
	do_judgement

}

run
