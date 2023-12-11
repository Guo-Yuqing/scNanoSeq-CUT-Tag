#!/bin/bash
workdir=`pwd`
cd $workdir
scripts=$workdir/pre_ASP_filter_reads.sh

cat $workdir/sample.list|grep -v "^#"|while read line
do

 echo "$line">./${date}_temp.txt
 batch=`awk '{print $1}' ./${date}_temp.txt`
 cell=`awk '{print $2}' ./${date}_temp.txt`

 rm ./${date}_temp.txt

run_log=$workdir/$batch
mkdir -p $run_log

echo "#!/bin/bash
sh $scripts  $batch $cell
" > $run_log/${cell}_SNP.sh

chmod +x $run_log/${cell}_SNP.sh

pkubatch \
   -J GYQ_${cell} \
   -p cn-long \
   -A tangfuchou_g1 \
   -c 10  \
   --qos=tangfuchoucnl \
   -o $run_log/${cell}_SNP.o -e $run_log/${cell}_SNP.e \
   $run_log/${cell}_SNP.sh


done






