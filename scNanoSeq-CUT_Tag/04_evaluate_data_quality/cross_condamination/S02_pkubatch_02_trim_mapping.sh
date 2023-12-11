#!/bin/bash
samdir=20230331_GM12878_3T3_H3K4me3_sample.list
workdir=`pwd`
cd $workdir
scripts=$workdir/run_scripts/cross_condamination/02_trim_mapping_hg38_mm10.sh

cat $workdir/samplelist/$samdir|grep -v "^#"|while read line
do
  echo "$line">./${date}_temp.txt
  cell=`awk '{print $1}' ./${date}_temp.txt`
  batch=`awk '{print $2}' ./${date}_temp.txt`
  prefix=`awk '{print $3}' ./${date}_temp.txt`
  
  rm ./${date}_temp.txt
 
 run_log=$workdir/run_scripts/cross_condamination/run_log/$batch
 mkdir -p $run_log
 
 echo "#!/bin/bash
 sh $scripts  $cell $batch $prefix
 " > $run_log/${cell}_02_trim_mapping.sh
 
 chmod +x $run_log/${cell}_02_trim_mapping.sh
 
 pkubatch \
    -J GYQ_cross_${cell} \
    -p cn-long \
    -A tangfuchou_g1 \
    -c 10  \
    --qos=tangfuchoucnl \
    -o $run_log/${cell}_02_trim_mapping.o -e $run_log/${cell}_02_trim_mapping.e \
    $run_log/${cell}_02_trim_mapping.sh

done
