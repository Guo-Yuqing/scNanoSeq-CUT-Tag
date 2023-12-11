#!/bin/sh
rootdir=~
workdir=$rootdir/project
cd $workdir

cell=$1
batch=$2
prefix=$3

demultiplex_dir=$workdir/01_demultiplex/$batch
######################################################################################################################################################
QC_dir=$workdir/02_trim/$batch
align_dir=$workdir/03_alignment/$batch
mkdir -p $QC_dir $align_dir

ref=$rootdir/reference/minimap2/genome/Homo_sapiens.GRCh38.dna.primary_assembly.mmi
ref_genome_fa=$rootdir/reference/minimap2/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
NanoFilt=$rootdir/software/miniconda3/bin/NanoFilt
NanoStat=$rootdir/software/miniconda3/bin/NanoStat
cutadapt=$rootdir/software/miniconda3/bin/cutadapt
minimap2=$rootdir/software/minimap2-2.26_x64-linux/minimap2
samtools=$rootdir/software/miniconda3/bin/samtools
cuteSV=$rootdir/software/miniconda3/envs/DNAloop/bin/cuteSV
scripts=$workdir/run_scripts/02_trim_mapping.sh


do_QC(){
# QC and Trim
cd $QC_dir/
mkdir -p ${cell}
cd ${cell}/
cp $scripts .

if [ -f ${demultiplex_dir}/${cell}.fastq.gz ];then
    gunzip ${demultiplex_dir}/${cell}.fastq.gz
else
    echo "done"
fi

echo -e "["$(date)"]\tTrimming adaptors..."
$NanoFilt -q 7 -l 100 ${demultiplex_dir}/${cell}.fastq > ${prefix}_${cell}_QC.fastq
$NanoStat --fastq ${prefix}_${cell}_QC.fastq -t 8 -o ./ -n ${prefix}_${cell}_stat.txt

$cutadapt \
    -e 0.25 \
    -j 8 \
    -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG...CTGTCTCTTATACACATCTGACGCTGCCGACGA \
    -o ${prefix}_${cell}_trimmed.fastq \
    ${prefix}_${cell}_QC.fastq

rm -f ${prefix}_${cell}_QC.fastq
cd ../

}

do_mapping(){
# minimap2 mapping
mkdir -p $align_dir/${prefix}_${cell}
echo -e "["$(date)"]\tMapping to reference genome..."

cd $align_dir/${prefix}_${cell}/

cd $QC_dir/${cell}
cp $scripts .
####remove supplementary mapping
$minimap2 -t 8 \
        --MD \
		-ax map-ont \
         ${ref} \
         $QC_dir/${cell}/${prefix}_${cell}_trimmed.fastq \
       |  $samtools view  -ShuF 2308 -q 30 -|\
         $samtools sort - -o $align_dir/${prefix}_${cell}/${prefix}_${cell}_q30.bam

gzip $QC_dir/${cell}/${prefix}_${cell}_trimmed.fastq

$samtools rmdup -s $align_dir/${prefix}_${cell}/${prefix}_${cell}_q30.bam $align_dir/${prefix}_${cell}/${prefix}_${cell}_q30.rmdup.temp.bam
$samtools view -F 2308 -c $align_dir/${prefix}_${cell}/${prefix}_${cell}_q30.bam > $align_dir/${prefix}_${cell}/${prefix}_${cell}_MappedReads

#rename bam header
$rootdir/software/miniconda3/bin/picard AddOrReplaceReadGroups \
       I=$align_dir/${prefix}_${cell}/${prefix}_${cell}_q30.rmdup.temp.bam \
       O=$align_dir/${prefix}_${cell}/${prefix}_${cell}_q30.rmdup.bam \
       RGID=${prefix}_${cell} \
       RGLB=$batch \
       RGPL=ONT \
       RGPU=unit1 \
       SORT_ORDER=coordinate \
       RGSM=${prefix}_${cell}

rm $align_dir/${prefix}_${cell}/${prefix}_${cell}_q30.rmdup.temp.bam
# bam to bed(fragment)
echo -e "["$(date)"]\tCreating fragments file..."
bedtools bamtobed -i $align_dir/${prefix}_${cell}/${prefix}_${cell}_q30.rmdup.bam | awk '{print $1"\t"$2"\t"$3"\t""'"$prefix"'""_""'"$cell"'""\t"1}' > $align_dir/${prefix}_${cell}/${prefix}_${cell}.rmdup_fragments
echo -e "["$(date)"]\tFinished!"

if [[ "${prefix}" =~ "GM12878" ]]
then
    echo "GM12878 cellline"
else
    rm -f $align_dir/${prefix}_${cell}/${prefix}_${cell}_q30.bam
fi


cd $workdir
}


do_countStat(){
demultiplex_dir=$workdir/01_demultiplex/$batch
QC_dir=$workdir/02_trim/$batch
align_dir=$workdir/03_alignment/$batch
count_dir=$workdir/SS_count/$batch
mkdir -p $count_dir

samtools=$rootdir/software/miniconda3/bin/samtools

Raw_Reads=$(zcat ${demultiplex_dir}/${cell}.fastq.gz | grep runid| wc -l)
Q7_Reads=$(cat ${QC_dir}/${cell}/${prefix}_${cell}_stat.txt | grep Q7 | awk 'BEGIN{IFS='\t'} {print $2}')
Q7_Percent=$(printf "%.2f" `echo "scale=2;100*${Q7_Reads}/${Raw_Reads}"|bc`)
Q10_Reads=$(cat ${QC_dir}/${cell}/${prefix}_${cell}_stat.txt | grep Q10 | awk 'BEGIN{IFS='\t'} {print $2}')
Q10_Percent=$(printf "%.2f" `echo "scale=2;100*${Q10_Reads}/${Raw_Reads}"|bc`)
Q15_Reads=$(cat ${QC_dir}/${cell}/${prefix}_${cell}_stat.txt | grep Q15 | awk 'BEGIN{IFS='\t'} {print $2}')
Q15_Percent=$(printf "%.2f" `echo "scale=2;100*${Q15_Reads}/${Raw_Reads}"|bc`)

Reads_median_length=$(cat ${QC_dir}/${cell}/${prefix}_${cell}_stat.txt |grep 'Median read length' | awk 'BEGIN{IFS='\t'} {print $4}'|sed 's/,//g')
total_base=$(cat ${QC_dir}/${cell}/${prefix}_${cell}_stat.txt |grep 'Total bases' | awk 'BEGIN{IFS='\t'} {print $3}'|sed 's/,//g')
Seq_data=$(printf "%.4f" `echo "scale=4;${total_base}/1000000000"|bc`)

Mapped_Reads=$(cat ${align_dir}/${prefix}_${cell}/${prefix}_${cell}_MappedReads)
Mapped_rmdup_Reads=$($samtools view -F2308 -c ${align_dir}/${prefix}_${cell}/${prefix}_${cell}_q30.rmdup.bam)


Mapped_percent=$(printf "%.2f" `echo "scale=2;100*${Mapped_Reads}/${Q7_Reads}"|bc`)
Mapped_rmdup_percent=$(printf "%.2f" `echo "scale=2;100*${Mapped_rmdup_Reads}/${Mapped_Reads}"|bc`)

mt_rmdup_num=$(less ${align_dir}/${prefix}_${cell}/${prefix}_${cell}.rmdup_fragments|grep -w 'MT'|wc -l)
mt_rmdup_pert=$(printf "%.4f" `echo "scale=4;100*${mt_rmdup_num}/${Mapped_rmdup_Reads}"|bc`)


Coverage=$(less ${align_dir}/${prefix}_${cell}/${prefix}_${cell}.rmdup_fragments|awk '{print "chr"$1"\t"$2"\t"$3}'|mergeBed -i - |awk '{print $3-$2+1}'|awk '{sum +=$1}END{print sum}')
Depth=$(printf "%.4f" `echo "scale=4;(${Mapped_Reads}*${Reads_median_length})/${Coverage}"|bc`)

echo -e ${prefix}"_"${cell}'\t'${Raw_Reads}'\t'${Q7_Reads}'\t'${Q7_Percent}'\t'${Q10_Reads}'\t'${Q10_Percent}'\t'${Q15_Reads}'\t'${Q15_Percent}'\t'${Reads_median_length}'\t'${total_base}'\t'${Seq_data}'\t'${Mapped_Reads}'\t'${Mapped_rmdup_Reads}'\t'${Mapped_percent}'\t'${Mapped_rmdup_percent}'\t'${mt_rmdup_num}'\t'${mt_rmdup_pert}'\t'${Coverage}'\t'${Depth} > ${count_dir}/${cell}_ReadCount.txt


}


run(){
    do_QC
    do_mapping
	do_countStat
}
run

