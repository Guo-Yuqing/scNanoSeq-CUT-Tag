#!/bin/bash
samdir=sample_list_1_batch1
#samdir=sample_list_2_batch2
#samdir=sample_list_3_batch3
#samdir=sample_list_4_batch4
workdir=`pwd`
cd $workdir
scripts=$workdir/GIAB_12878_WGS_Unique_mapping_count.sh

cat $workdir/samplelist/$samdir|grep -v "^#"|grep gz|while read line
do
 echo "$line">./${date}_temp.txt
 batch1=`awk '{print $1}' ./${date}_temp.txt`
 sample1=`awk '{print $2}' ./${date}_temp.txt`
 sample2=`awk '{print $3}' ./${date}_temp.txt`
 rm ./${date}_temp.txt

run_log=$workdir/run_scripts/run_log/
mkdir -p $run_log

name1=$(echo ${sample1} |awk -F'\/' '{print $12}'|sed 's/.fastq.gz//g')
name2=$(echo ${sample2} |awk -F'\/' '{print $12}'|sed 's/.fastq.gz//g')


echo "#!/bin/bash
sh $scripts ${batch1}  $sample1 $sample2
" > $run_log/GIAB_12878_WGS_Unique_mapping_${batch1}_${name1}_${name2}.sh


chmod +x $run_log/GIAB_12878_WGS_Unique_mapping_${batch1}_${name1}_${name2}.sh

pkubatch \
    -J GYQ_WGS \
    -p cn-long \
    -A tangfuchou_g1 \
    -c 10  \
    --qos=tangfuchoucnl \
    -o $run_log/GIAB_12878_WGS_Unique_mapping_${batch1}_${name1}_${name2}.o -e $run_log/GIAB_12878_WGS_Unique_mapping_${batch1}_${name1}_${name2}.e \
    $run_log/GIAB_12878_WGS_Unique_mapping_${batch1}_${name1}_${name2}.sh


done
