#!/bin/bash
batch=$1
sample1=$2
sample2=$3

rootdir=~
workdir=$rootdir/project/GIAB_GM12878_WGS
cd $workdir

SS_count_dir=$workdir/SS_count_dir

name1=$(echo ${sample1} |awk -F'\/' '{print $12}'|sed 's/.fastq.gz//g')
name2=$(echo ${sample2} |awk -F'\/' '{print $12}'|sed 's/.fastq.gz//g')

sample=$(echo ${name1}|sed 's/_R[0-9]//g')

fastp=$rootdir/software/fastp #version 0.23.0
bwa=$rootdir/software/bwa-0.7.17/bwa #0.7.17
sambamba=$rootdir/software/miniconda3/bin/sambamba
bedtools=$rootdir/software/miniconda3/bin/bedtools
picard=$rootdir/software/miniconda3/bin/picard
samtools=$rootdir/software/miniconda3/bin/samtools

ref_fa=$rootdir/reference/index/bwa/T2T/chm13v2.0.fa

raw_fq1=$workdir/raw_fq/${batch}/${name1}.fastq.gz
raw_fq2=$workdir/raw_fq/${batch}/${name2}.fastq.gz
clean_fq1=$workdir/raw_fq/${batch}/${name1}.Clean.fastq.gz
clean_fq2=$workdir/raw_fq/${batch}/${name2}.Clean.fastq.gz
map_Q30=$workdir/Mapping_dir/${batch}/${sample}/${sample}.Q30.bam
map_Q30_unique=$workdir/Mapping_dir/${batch}/${sample}/${sample}.Q30.unique.bam
map_Q30_unique_sort=$workdir/Mapping_dir/${batch}/${sample}/${sample}.Q30.unique.sort.bam
map_Q30_unique_rmdup=$workdir/Mapping_dir/${batch}/${sample}/${sample}.Q30.unique.sort.rmdup.bam



do_clean(){
	$fastp -i $raw_fq1                                       \
	-I $raw_fq2                                     \
    -o $clean_fq1                                           \
    -O $clean_fq2

}


do_mapping(){
#+++++++++++++++++++++
# BWA Mapping
#+++++++++++++++++++++
	mkdir -p $workdir/Mapping_dir/${batch}/${sample}
	cd $workdir/Mapping_dir/${batch}/${sample}
    echo -e "["$(date)"]\tStart $sample aligning.."
    $bwa mem -M -t 6                                                           \
            $ref_fa                                                           \
            $clean_fq1                                            \
            $clean_fq2   | \
         $samtools view -bS -F2308 -q 30 - > $map_Q30

	$sambamba view -t 8 -h  \
	-f bam -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)"  \
	$map_Q30   -o $map_Q30_unique

	$samtools sort -m 10G -o $map_Q30_unique_sort  $map_Q30_unique
	$bedtools bamtobed -i $map_Q30_unique_sort  -tag NM  >$workdir/Mapping_dir/${batch}/${sample}/${sample}.Q30.unique.bed

	$picard MarkDuplicates                                          \
      I=$map_Q30_unique_sort                            \
      O=$map_Q30_unique_rmdup           \
      M=${sample}_dup.metrics                                     \
      CREATE_INDEX=true                                                    \
      ASSUME_SORTED=true                                                   \
      VALIDATION_STRINGENCY=SILENT                                         \
      REMOVE_DUPLICATES=true

	$bedtools bamtobed -i $map_Q30_unique_rmdup  -tag NM  >$workdir/Mapping_dir/${batch}/${sample}/${sample}.Q30.unique.sort.rmdup.bed

    echo -e "["$(date)"]\t$sample BWA done!"

}

do_count(){
    Raw1_Reads_num=$(zcat $raw_fq1 |  grep HWI| wc -l)
    Raw2_Reads_num=$(zcat $raw_fq2 |  grep HWI| wc -l)
    Total_raw_Reads_num=$((${Raw1_Reads_num} + ${Raw2_Reads_num}))
    Clean1_Reads_num=$(zcat $clean_fq1 |  grep HWI| wc -l)
    Clean2_Reads_num=$(zcat $clean_fq1 |  grep HWI| wc -l)
    Total_clean_Reads_num=$((${Clean1_Reads_num} + ${Clean2_Reads_num}))
    Reads_clean_percent=$(printf "%.2f" `echo "scale=2;100*${Total_clean_Reads_num}/${Total_raw_Reads_num}"|bc`)
    Mapped_Q30_Reads=$($samtools view -F2308 -c $map_Q30)
    Mapped_percent=$(printf "%.2f" `echo "scale=2;100*${Mapped_Q30_Reads}/${Total_clean_Reads_num}"|bc`)
    Mapped_Q30_unique_Reads=$($samtools view -F2308 -c $map_Q30_unique_sort)
    Mapped_Q30_rmdup_Reads=$($samtools view -F2308 -c $map_Q30_unique_rmdup)
    Mapped_Q30_rmdup_percent=$(printf "%.2f" `echo "scale=2;100*${Mapped_Q30_rmdup_Reads}/${Mapped_Q30_Reads}"|bc`)
	genome_cov_dedup=$(less $workdir/Mapping_dir/${batch}/${sample}/${sample}.Q30.unique.sort.rmdup.bed|mergeBed -i -|awk '{print $3-$2}' - | awk '{sum+=$1}END{print sum}')
    echo -e ${batch}_${sample}'\t'${Raw1_Reads_num}'\t'${Raw2_Reads_num}'\t'${Total_raw_Reads_num}'\t'${Clean1_Reads_num}'\t'${Clean2_Reads_num}'\t'${Total_clean_Reads_num}'\t'${Reads_clean_percent}'\t'${Mapped_Q30_Reads}'\t'${Mapped_percent}'\t'${Mapped_Q30_rmdup_Reads}'\t'${Mapped_Q30_rmdup_percent}'\t'${genome_cov_dedup}> $SS_count_dir/${batch}_${sample}_ReadCount.txt

}


do_rm(){

rm $clean_fq1 $clean_fq2 $workdir/Mapping_dir/${batch}/${sample}/${sample}.Q30.bam $workdir/Mapping_dir/${batch}/${sample}/${sample}.Q30.unique.bam
rm $map_Q30_unique_sort  
rm $workdir/Mapping_dir/${batch}/${sample}/${sample}.Q30.unique.bam 
}



do_L1HS_merge_noerror(){

cd $workdir/Mapping_dir/${batch}/${sample}


bin_file=/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/guoyuqing/DNAloop/project/000_Merge_analysis/02_L1HS_fulllength_bin_new/L1HS/blast_result/QUERY_index_Bin_300bp

less -S ${sample}.Q30.unique.sort.rmdup.bed|  \
awk '{if ($5==0)print}'| \
bedtools intersect -a $bin_file -b - -wo| \
awk '{print $0"\t"$8-$7}'|awk '{if ($12==$13)print}'| \
#cut -f 1,2,3,4,5|uniq| \
awk -vOFS='\t'  '{print $1,$2,$3,$4,$5,1}'|bedtools groupby -g 1,2,3,4,5 -c 6 -o sum|\
bedtools intersect -a $bin_file -b - -loj|\
awk -vOFS='\t'  '{if ($7>0)print $1,$2,$3,$4,$5,$11;else print $1,$2,$3,$4,$5,0}'| \
sort -k 4,4V -k 5,5V |uniq>${sample}.error0_L1HS_300bp_bin.bed

}



run(){
	do_clean
	do_mapping
	do_count
	do_rm
	do_L1HS_merge_noerror
}

run

