#!/bin/bash
workdir=~

do_prepare(){
    zless -S  mouse_rmsk.txt.gz| \
    awk -vOFS='\t' '{print $6,$7,$8,$10,$11,$12,$13}'|  \
    grep LINE|bedtools intersect -a - -b l1base_mouse.sort.bed -wo|  \
    awk '{if ($17>3500)print}'|bedtools groupby -g 1,2,3,4,5,6,7 -c 17 -o max|  \
    bedtools intersect -a - -b l1base_mouse.sort.bed -wo|  \
    awk '{if ($8==$18)print}'|grep -E "A|T|Gf"| \
    cut -f 1,2,3,4,5,6,7,8|sort -k 1,1V -k 2,2n -k 3,3n|uniq >trans_rmsk_l1base_mouse.sort.bed

    mkdir Repeat_A
    less trans_rmsk_l1base_mouse.sort.bed|grep A|sed 's/chr//g' >Repeat_A/trans_rmsk_l1base_mouse.sort.L1Md_A.bed
    cd Repeat_A
    makeblastdb -in AY053456.fasta -dbtype nucl  -parse_seqids  -out  AY053456_database
    less trans_rmsk_l1base_mouse.sort.L1Md_A.bed|awk -vOFS='\t' '{print $1,$2-2000,$3+2000,$4,$5,$6,$7}'  >trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.bed


    bedtools getfasta -fi  ../mm10.fa -bed trans_rmsk_l1base_mouse.sort.L1Md_A.bed -fo trans_rmsk_l1base_mouse.sort.L1Md_A.bed.fasta
    bedtools getfasta -fi  ../mm10.fa -bed  trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.bed -fo trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.bed.fasta

    blastn  -query trans_rmsk_l1base_mouse.sort.L1Md_A.bed.fasta -out trans_rmsk_l1base_mouse.sort.L1Md_A.blast_result -db AY053456_database -xdrop_gap 4000  -max_hsps 1 -outfmt "7 std stitle"  &
    blastn -query  trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.bed.fasta  -out  trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast_result  -db AY053456_database -xdrop_gap 4000  -max_hsps 1 -outfmt "7 std stitle"  &
}


do_blast(){
    cd $workdir/Repeat_A
    mkdir flank2K
    blastn -query trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.bed.fasta -out flank2K/trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast -db  AY053456_database -xdrop_gap 4000  -max_hsps 1   
    mkdir flank0
    blastn -query  trans_rmsk_l1base_mouse.sort.L1Md_A.bed.fasta -out  flank0/trans_rmsk_l1base_mouse.sort.L1Md_A.blast -db  AY053456_database -xdrop_gap 4000  -max_hsps 1   


}


active_fa=$workdir/Repeat_A/AY053456.fasta
all_L1_bed=$workdir/Repeat_A/trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.bed


do_build_index_new(){
    fmt_downdir=$workdir/Repeat_A/flank2K/fmt_down_analysis
    mkdir -p $fmt_downdir
    cd $fmt_downdir

    blast_file=$workdir/Repeat_A/flank2K/trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast
    ln ${blast_file} .

    less  -S  ${active_fa}|  \
    grep -v ">"|sed 's/./&\n/g'|awk '{print toupper($0)}'|grep -v '^$'|  \
    awk '$0=NR"\t"$0'|awk '{print "Base"$0}'|  \
    sed '1i#Base\tref_active' >$fmt_downdir/mm10_L1Md_A.fulllength_active.matrix_file

    less ${all_L1_bed}|awk '$0=NR"\t"$0' -|awk -vOFS='\t'  '{print $2":"$3"-"$4,$5,$6,$7,$8,$9,"Query_"$1}' >$fmt_downdir/QUERY_index

    ln -s $workdir/Repeat_A/trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast_result .
    less trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast_result|grep -v '#'|awk '{if ($4>4800)print}'|cut -f1|sort -k 1,1V|uniq|awk '{print $0}' >true_blast_query


    end=$(less trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast|wc -l)
    awk '$0=NR"\t"$0'  ${blast_file} | \
    grep Query=  >trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast_index_tmp

    less trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast_index_tmp|cut -f1| \
    tail -n +2|sed  '$a'$end'' >trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast_index_tmp_end

    paste -d '\t' trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast_index_tmp trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast_index_tmp_end| \
    awk -vOFS='\t' '{print $1,$4,$2,$3}' |sed 's/Query=//g'|  \
    awk 'NR==FNR{a[$1]=1}NR!=FNR{if($3 in a)print $0}' true_blast_query - |\
    awk '{print $1"\t"$2"\t""Query=""\t"$3}'|awk '{if ($2-$1>10)print}' >trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast_index


    less trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast_index|while read line
    do
    echo "$line">./${date}_temp.txt
    start_index=`awk '{print $1}' ./${date}_temp.txt`
    end_index=`awk '{print $2}' ./${date}_temp.txt`
    query=`awk '{print $3}' ./${date}_temp.txt`
    resion=`awk '{print $4}' ./${date}_temp.txt`
    rm ./${date}_temp.txt

    region_range=$(($end_index-$start_index))
    strand=$(less ${blast_file}|tail -n  +${start_index} -|head -n ${region_range}|grep Strand|sed 's/"Strand=Minus\/"//g')

    echo -e  $start_index"\t"$end_index"\t"$resion"\t"$strand >>trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast_index_anno
    done
    rm trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast_index_tmp  trans_rmsk_l1base_mouse.sort.L1Md_A.flank2K.blast_index_tmp_end

}

##############################################################################################################################################################################################################################################################
do_get_blast_diff_matrix(){
   fmt_downdir=$workdir/Repeat_A/flank2K/fmt_down_analysis 
   cd $fmt_downdir

    plus_scripts=build_blast_diff_matrix.plus.sh
    minus_scripts=build_blast_diff_matrix_minus.sh
    r_scripts=motif_ref_query_diff.new.R


    bash $plus_scripts
    bash $minus_scripts

}


do_merge(){
    fmt_downdir=$workdir/Repeat_A/flank2K/fmt_down_analysis
    merge_dir=$workdir/Repeat_A/flank2K/fmt_down_analysis/detail/merge
    mkdir -p $merge_dir
    cd $merge_dir
    wc -l ../*/*_ref_query_diff_final_modify |grep _modify|awk '{if ($1==4962)print $2}'|while read id ;do (ln -s $id .);done
    ln -s $fmt_downdir/mm10_L1Md_A.fulllength_active.matrix_file .

    cp mm10_L1Md_A.fulllength_active.matrix_file mm10_L1Md_A.fulllength_all_blast

    ls *_modify|while read id
    do

    less $id|cut -f3|paste -d '\t' mm10_L1Md_A.fulllength_all_blast - >tmp

    mv tmp mm10_L1Md_A.fulllength_all_blast

    done

}

do_get_Bin_info(){

    bin_region=$workdir/Repeat_A/Repeat_bin/Repeat_region.new2.bin.bed 
    outdir=$workdir/Repeat_A/flank2K/fmt_down_analysis/detail/merge/repeat_bin/
    cd $workdir/Repeat_A/flank2K/fmt_down_analysis/detail/merge
    mkdir -p $outdir
    ls *_modify|while read id
    do
        name=$(echo $id|sed 's/_ref_query_diff_final_modify//g')
        echo $name

        for bin in "Bin_"{1..17}
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
    for bin in "Bin_"{1..17}
    do
    not=$(less ${bin}|sed 's/,/_/g'|sort -k 5,5V|uniq|awk '{print $0"\t"1}'|bedtools groupby -g 1,2,3,5 -c 6 -o sum|sort -k 5,5rn|awk '{if ($5>1)print $5}'|awk '{sum+=$1}END{print sum}')
    yes=$(less ${bin}|sed 's/,/_/g'|sort -k 5,5V|uniq|awk '{print $0"\t"1}'|bedtools groupby -g 1,2,3,5 -c 6 -o sum|sort -k 5,5rn|awk '{if ($5==1||$5==NA)print $5}'|awk '{sum+=$1}END{print sum}')

    echo -e ${bin}"\t"${not}"\t"${yes} >>final_file
    done

}




run(){
	do_prepare
	do_blast
	do_build_index_new
	do_get_blast_diff_matrix
	do_merge
	do_get_Bin_info
}


run

