#!/bin/bash
celltype=$1
sample=$2

rootdir=~
workdir=$rootdir/project/ASP
raw_map=$workdir/Map_dir/$celltype
SNP_dir=$workdir/SNP_dir/$celltype
SNP_filter_dir=$workdir/SNP_filter_dir/$celltype/${sample}

mkdir -p $raw_map $SNP_dir $SNP_filter_dir

ref=$rootdir/reference/minimap2/genome/Homo_sapiens.GRCh38.dna.primary_assembly.mmi
ref_genome_fa=$rootdir/reference/minimap2/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
minimap2=$rootdir/software/miniconda3/bin/minimap2
samtools=$rootdir/software/miniconda3/bin/samtools
picard=$rootdir/software/miniconda3/bin/picard
bedtools=$rootdir/software/miniconda3/bin/bedtools
whatshap=$rootdir/software/miniconda3/envs/DNAloop/bin/whatshap
HG001_SNP=$rootdir/project/data_base/GM12878/HG001_GRCh38_GIAB_highconf_nochr_noPATMAT.vcf.gz
hetSNP_dir=$rootdir/project/data_base/GM12878/split_hetSNP_bed
python_tools=$rootdir/software/miniconda3/bin/python
subset_scripts=$workdir/pull_reads_new.py


do_subset_hetSNP(){
cd $raw_map/
	######################
    #split pat mat reads
    ######################
    sample_dir=$SNP_dir/$sample
    mkdir -p $sample_dir
    cd $sample_dir
    ln -s $raw_map/${sample}_q30.bam .
     $samtools index ${sample}_q30.bam

    for chr_ in  `seq 1 22` X
    do
        $whatshap haplotag \
            -o ${sample}.${chr_}.bam \
            --reference $ref_genome_fa \
            --ignore-read-groups \
            --regions $chr_ \
            ${HG001_SNP} \
            ${sample}_q30.bam

        $bedtools bamtobed -i ${sample}.${chr_}.bam >${sample}.${chr_}.tmp.bed
		
		$samtools view -h ${sample}.${chr_}.bam | \
		awk -vOFS='\t' 'match($0,/PS:i:([0-9])/,PS);match($0,/HP:i:([0-9])/,HP);{if (PS[1] != "" && HP[1] != "")print $1,$2,$3,$4,PS[1],HP[1]}'| \
		awk '{if ($5==1)print "chr"$3"\t"$4-1"\t"$4"\t"$1"\t""'$sample'""\t"$5"\t"$6}'>${sample}.${chr_}.PATMAT

        less ${sample}.${chr_}.tmp.bed|awk '{print "chr"$0}'|$bedtools intersect -a - -b ${sample}.${chr_}.PATMAT -wa -wb| \
        awk -vOFS='\t' '{if ($2==$8 && $4==$10)print $1,$2,$3,$4,$11,$12,$13}'|sed 's/_q30.rmdup//g' >${sample}.${chr_}.PATMAT.bed

        rm ${sample}.${chr_}.tmp.bed  ${sample}.${chr_}.PATMAT
         $samtools index ${sample}.${chr_}.bam

    ######################
    #subset reads SNPs
    ######################
        $python_tools $subset_scripts  \
        --bam ${sample}.${chr_}.bam  \
        --vcf $hetSNP_dir/chr${chr_}_hetSNP.forVCF.bed  \
        --out ${sample}.${chr_}
    done

cd $raw_map/
}


do_filter_HetSNPs(){
cd $SNP_dir/${sample}

for chr_ in  `seq 1 22` X
do
######################
    #filter reads SNPs
######################
less ${sample}.${chr_}.tsv|cut -f 1,2,6,7| \
sort -k 1,1V -k 3,3V -k 4,4V|grep -v  Error| \
awk '{print $0"\t"1}'|  \
awk '{if ($4=="PAT")print $0"\t"1;else  print $0"\t"0}'|  \
bedtools groupby -g 1,3 -c 5,6 -o sum,sum|  \
awk '{if ($3>1)print}'|awk '{if ($4/$3>0.9 ||($3-$4)/$3>0.9)print}'| \
awk 'NR==FNR{a[$2]=1}NR!=FNR{if($4 in a)print $0}' - ${sample}.${chr_}.PATMAT.bed  >tmp1


less ${sample}.${chr_}.tsv|cut -f 1,2,6,7| \
sort -k 1,1V -k 3,3V -k 4,4V|grep -v  Error|awk '{print $0"\t"1}'| \
awk '{if ($4=="PAT")print $0"\t"1;else  print $0"\t"0}'|bedtools groupby -g 1,3 -c 5,6 -o sum,sum| \
awk '{if ($3==1)print}' | \
awk 'NR==FNR{a[$2]=1}NR!=FNR{if($6 in a)print $0}' - ${sample}.${chr_}.tsv | \
awk '{print $0"\t"$2"_"$7}'>tmp2_1snp_reads


awk 'NR==FNR{a[$2]=1}NR!=FNR{if($2 in a)print $0}' tmp2_1snp_reads ${sample}.${chr_}.tsv | \
cut -f 1,2,7 |grep -v Error|awk '{if ($3=="PAT")print $0"\t"1;else print $0"\t"0}'| \
bedtools groupby -g 1,2 -c 3,4 -o count,sum| \
awk '{if ($4/$3<=0.1)print $0"\t""MAT";else if($4/$3>=0.9)print $0"\t""PAT"}'| \
awk '{print $2"_"$5}'|awk 'NR==FNR{a[$1]=1}NR!=FNR{if($8 in a)print $0}' - tmp2_1snp_reads| \
awk 'NR==FNR{a[$6]=1}NR!=FNR{if($4 in a)print $0}' - ${sample}.${chr_}.PATMAT.bed >>tmp1

less tmp1|cut -f 1,2,3,5,6,7|sort -k 1,1V -k 2,2n -k 3,3n -k 6,6n|uniq| \
awk '$0=NR"\t"$0' - |awk -vOFS='\t' '{print $2,$3-50,$3+50,"Reads_"$1,$5,$6,$7}' >tmp_unique
less tmp1|cut -f 1,2,3,5,6,7|sort -k 1,1V -k 2,2n -k 3,3n -k 6,6n|uniq| \
awk '$0=NR"\t"$0' - |awk -vOFS='\t' '{print $2,$4-50,$4+50,"Reads_"$1,$5,$6,$7}' >>tmp_unique

cp tmp1 $SNP_filter_dir/${sample}.${chr_}.PATMAT.filter.reads
cp tmp_unique $SNP_filter_dir/${sample}.${chr_}.PATMAT.filter.bed
rm tmp1 tmp2_1snp_reads tmp_unique
rm  ${sample}.${chr_}.bam  ${sample}.${chr_}.bam.bai 

done


}

#after all single cells be done

do_merge(){
cd $SNP_filter_dir
cat */*/*.bed |sort -k 1,1V -k 2,2n -k 3,3n|bgzip  >GM12878_H3K4me3_PATMAT_flank100.gz

}




run(){
	do_subset_hetSNP
	do_filter_HetSNPs
#	do_merge
}

run







