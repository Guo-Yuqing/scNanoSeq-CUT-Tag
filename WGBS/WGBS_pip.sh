#!/bin/bash
sample=$1
rootdir=~
workdir=$rootdir/project/WGBS/
cd $workdir

fastp=$rootdir/software/fastp #version 0.23.0
bismark=$rootdir/software/miniconda3/bin/bismark
samtools=$rootdir/software/miniconda3/bin/samtools
bedtools=$rootdir/software/miniconda3/bin/bedtools
bismark_methylation_extractor=$rootdir/software/miniconda3/bin/bismark_methylation_extractor
coverage2cytosine=$rootdir/software/miniconda3/bin/coverage2cytosine

raw_dir=$workdir/${sample}
clean_dir=$workdir/${sample}/01_clean_fq
map_dir=$workdir/${sample}/02_Mapping
map_lamada_dir=$workdir/${sample}/02_Mapping_lamada
meth_dir=$workdir/${sample}/03_methylation
meth_lamada_dir=$workdir/${sample}/03_methylation_lamada
SS_count_dir=$workdir/SS_count
mkdir -p $clean_dir
mkdir -p $map_dir $map_lamada_dir
mkdir -p $meth_dir  $meth_lamada_dir
mkdir -p $SS_count_dir

raw_fq1=$workdir/${sample}/${sample}_merge_1.fq.gz
raw_fq2=$workdir/${sample}/${sample}_merge_2.fq.gz
clean_fq1=$clean_dir/${sample}_1.Clean.fq.gz
clean_fq2=$clean_dir/${sample}_2.Clean.fq.gz
map_R1=$map_dir/${sample}_1.Clean_bismark_bt2.bam
map_R2=$map_dir/${sample}_2.Clean_bismark_bt2.bam
map_sort_R1=$map_dir/${sample}_1.Clean_bismark_bt2.sort.bam
map_sort_R2=$map_dir/${sample}_2.Clean_bismark_bt2.sort.bam
map_sort_rmdup_R1=$map_dir/${sample}_1.Clean_bismark_bt2.sort.rmdup.bam
map_sort_rmdup_R2=$map_dir/${sample}_2.Clean_bismark_bt2.sort.rmdup.bam
merge_bam=$map_dir/${sample}_merge.sort.rmdup.bam
merge_bed=$map_dir/${sample}_merge.sort.rmdup.bed
unmap_R1=$map_dir/${sample}_1.Clean.fq.gz_unmapped_reads.fq.gz
unmap_R2=$map_dir/${sample}_2.Clean.fq.gz_unmapped_reads.fq.gz
lamada_map_R1=$map_lamada_dir/${sample}_1.Clean.fq.gz_unmapped_reads_bismark_bt2.bam
lamada_map_R2=$map_lamada_dir/${sample}_2.Clean.fq.gz_unmapped_reads_bismark_bt2.bam
lamada_map_R1_sort=$map_lamada_dir/${sample}_1.Clean.fq.gz_unmapped_reads_bismark_bt2.sort.bam
lamada_map_R2_sort=$map_lamada_dir/${sample}_2.Clean.fq.gz_unmapped_reads_bismark_bt2.sort.bam
lamada_map_R1_sort_rmdup=$map_lamada_dir/${sample}_1.Clean.fq.gz_unmapped_reads_bismark_bt2.sort.rmdup.bam
lamada_map_R2_sort_rmdup=$map_lamada_dir/${sample}_2.Clean.fq.gz_unmapped_reads_bismark_bt2.sort.rmdup.bam
lamada_merge_bam=$map_lamada_dir/${sample}_lamada_merge.sort.rmdup.bam
lamada_merge_bed=$map_lamada_dir/${sample}_lamada_merge.sort.rmdup.bed

do_clean(){

    $fastp  -i $raw_fq1                                       \
            -I $raw_fq2                                     \
            -o $clean_fq1                                           \
            -O $clean_fq2

}
do_mapping(){

    $bismark --genome_folder $rootdir/reference/genome/hg38_no_alt  \
            --bowtie2 --fastq --non_directional --unmapped  \
            --path_to_bowtie2 $rootdir/software/bowtie2-2.3.5.1-linux-x86_64  \
            --samtools_path $rootdir/software/miniconda3/bin/samtools  \
            --output_dir $map_dir  \
            $clean_fq1
    $samtools sort $map_R1 -o $map_sort_R1
    $samtools rmdup -s  $map_sort_R1 $map_sort_rmdup_R1

    $bismark --genome_folder $rootdir/reference/genome/hg38_no_alt  \
            --bowtie2 --fastq --non_directional --unmapped  \
            --path_to_bowtie2 $rootdir/software/bowtie2-2.3.5.1-linux-x86_64  \
            --samtools_path $rootdir/software/miniconda3/bin/samtools  \
            --output_dir $map_dir  \
            $clean_fq2
    $samtools sort $map_R2 -o $map_sort_R2
    $samtools rmdup -s  $map_sort_R2 $map_sort_rmdup_R2

}

do_merge_bam(){

    $samtools merge $merge_bam  $map_sort_rmdup_R1  $map_sort_rmdup_R2
    $samtools index $merge_bam
    $bedtools bamtobed -i  $merge_bam >${merge_bed}

}

do_mapping_lamada(){

    $bismark --genome_folder $rootdir/reference/genome/hg38_lamada  \
            --bowtie2 --fastq --non_directional --unmapped  \
            --path_to_bowtie2 $rootdir/software/bowtie2-2.3.5.1-linux-x86_64  \
            --samtools_path $rootdir/software/miniconda3/bin/samtools  \
            --output_dir $map_lamada_dir  \
            $unmap_R1
    $samtools sort $lamada_map_R1 -o $lamada_map_R1_sort
    $samtools rmdup -s  $lamada_map_R1_sort $lamada_map_R1_sort_rmdup

    $bismark --genome_folder $rootdir/reference/genome/hg38_lamada  \
            --bowtie2 --fastq --non_directional --unmapped  \
            --path_to_bowtie2 $rootdir/software/bowtie2-2.3.5.1-linux-x86_64  \
            --samtools_path $rootdir/software/miniconda3/bin/samtools  \
            --output_dir $map_lamada_dir  \
            $unmap_R2

    $samtools sort $lamada_map_R2 -o $lamada_map_R2_sort
    $samtools rmdup -s  $lamada_map_R2_sort $lamada_map_R2_sort_rmdup


    $samtools merge $lamada_merge_bam  $lamada_map_R1_sort_rmdup  $lamada_map_R2_sort_rmdup
    $samtools index $lamada_merge_bam
    $bedtools bamtobed -i  $lamada_merge_bam >${lamada_merge_bed}

}


do_extract_meth(){
    $bismark_methylation_extractor   -s --ignore 6 --bedGraph  \
                                    --CX_context --zero_based  \
                                    --samtools_path $rootdir/software/miniconda3/bin/samtools  \
                                    --buffer_size 10G \
                                    --output $meth_dir \
                                    --genome_folder $rootdir/reference/genome/hg38_no_alt  \
                                    $merge_bam

    $coverage2cytosine  --dir $meth_dir   \
                        --gc \
                         --nome-seq \
                        --genome_folder $rootdir/reference/genome/hg38_no_alt  \
                        -o ${sample}   \
                        $meth_dir/${sample}_merge.sort.rmdup.bismark.cov.gz

}


do_extract_meth_lamada(){

    $bismark_methylation_extractor   -s --ignore 6 --bedGraph  \
                                    --CX_context --zero_based  \
                                    --samtools_path $rootdir/software/miniconda3/bin/samtools  \
                                    --buffer_size 10G \
                                    --output $meth_lamada_dir \
                                    --genome_folder $rootdir/reference/genome/hg38_lamada  \
                                    $lamada_merge_bam

    $coverage2cytosine  --dir $meth_lamada_dir   \
                        --gc \
                         --nome-seq \
                        --genome_folder $rootdir/reference/genome/hg38_lamada  \
                        -o ${sample}_lamada   \
                        $meth_lamada_dir/${sample}_lamada_merge.sort.rmdup.bismark.cov.gz

}


do_count_state(){

    Raw1_Reads_num=$(zcat $raw_fq1 | grep '@'| wc -l)
    Raw2_Reads_num=$(zcat $raw_fq2 | grep '@'| wc -l)
    Total_raw_Reads_num=$((${Raw1_Reads_num} + ${Raw2_Reads_num}))
    Clean1_Reads_num=$(zcat $clean_fq1 | grep '@'| wc -l)
    Clean2_Reads_num=$(zcat $clean_fq2 | grep '@'| wc -l)
    Total_clean_Reads_num=$((${Clean1_Reads_num} + ${Clean2_Reads_num}))
    Reads_clean_percent=$(printf "%.2f" `echo "scale=2;100*${Total_clean_Reads_num}/${Total_raw_Reads_num}"|bc`)
    map_R1_num=$(less $map_dir/${sample}_1.Clean_bismark_bt2_SE_report.txt|sed -n '8,8p'|awk '{print $13}')
    map_R2_num=$(less $map_dir/${sample}_2.Clean_bismark_bt2_SE_report.txt|sed -n '8,8p'|awk '{print $13}')
    Total_mapping_reads_num=$(($map_R1_num + $map_R2_num))
    Reads_mapping_ratio=$(printf "%.2f" `echo "scale=2;100*${Total_mapping_reads_num}/${Total_clean_Reads_num}"|bc`)
    rmdup_Reads_mapping_number=$($samtools view -c $merge_bam)
    Duplication=$(printf "%.2f" `echo "scale=2;100*${rmdup_Reads_mapping_number}/${Total_mapping_reads_num}"|bc`)
    cover_size=$(less ${merge_bed}|mergeBed -i - |awk '{print $3-$2}'|awk '{sum+=$1}END{print sum}')
    genome_size=$(less $rootdir/reference/genome/hg38_no_alt/hg38_no_alt_length.txt|awk '{sum+=$2}END{print sum}')
    Coverage=$(printf "%.2f" `echo "scale=2;100*${cover_size}/${genome_size}"|bc`)
    tmp1=$(less $meth_dir/${sample}.NOMe.CpG_report.txt|grep -v lambda|awk '{sum+=$4}END{print sum}')
    tmp2=$(less $meth_dir/${sample}.NOMe.CpG_report.txt|grep -v lambda|awk '{sum+=$5}END{print sum}')
    WCG_meth_MeanRatio=$(printf "%.2f" `echo "scale=2;100*${tmp1}/(${tmp1}+${tmp2})"|bc`)
    WCG_total_sites_num=$(less $meth_dir/${sample}.NOMe.CpG_report.txt|grep chr| wc -l)
    tmp3=$(less $meth_lamada_dir/${sample}*NOMe.CpG_report.txt|grep lambda|awk '{sum+=$4}END{print sum}')
    tmp4=$(less $meth_lamada_dir/${sample}*NOMe.CpG_report.txt|grep lambda|awk '{sum+=$5}END{print sum}')
    WCG_meth_MeanRatio_lamada=$(printf "%.2f" `echo "scale=2;100*${tmp4}/(${tmp3}+${tmp4})"|bc`)
    WCG_total_sites_num_lamada=$(less $meth_lamada_dir/${sample}*NOMe.CpG_report.txt|grep lambda| wc -l)

    echo -e ${sample}'\t'${Raw1_Reads_num}'\t'${Raw2_Reads_num}'\t'${Total_raw_Reads_num}'\t'${Clean1_Reads_num}'\t'${Clean2_Reads_num}'\t'${Total_clean_Reads_num}'\t'${Reads_clean_percent}'\t'${Total_mapping_reads_num}'\t'${Reads_mapping_ratio}'\t'${rmdup_Reads_mapping_number}'\t'${Duplication}'\t'${Coverage}'\t'${WCG_meth_MeanRatio}'\t'${WCG_total_sites_num}'\t'${WCG_meth_MeanRatio_lamada}'\t'${WCG_total_sites_num_lamada} > $SS_count_dir/${sample}_ReadCount.txt

}

run(){
    do_clean
    do_mapping
    do_merge_bam
    do_mapping_lamada
    do_extract_meth
    do_extract_meth_lamada
    do_count_state
}

run

