#!/bin/bash
samdir=tmp_nanoplex
workdir=`pwd`
cd $workdir
scripts=$workdir/run_scripts/01_demultiplex.sh

cat $workdir/samplelist/$samdir|grep -v "^#"|while read line
do
 echo "$line">./${date}_temp.txt
 sample=`awk '{print $1}' ./${date}_temp.txt`
 batch=`awk '{print $2}' ./${date}_temp.txt`
 rm ./${date}_temp.txt

run_log=$workdir/run_scripts/run_log/$batch
mkdir -p $run_log

echo "#!/bin/bash
sh $scripts  $sample $batch
" > $run_log/${batch}_01_demultiplex.sh

chmod +x $run_log/${batch}_01_demultiplex.sh

pkubatch \
    -J GYQ_demultiplex \
    -p cn-long \
    -A tangfuchou_g1 \
    -c 20  \
    --qos=tangfuchoucnl \
    -o $run_log/${batch}_01_demultiplex.o -e $run_log/${batch}_01_demultiplex.e \
    $run_log/${batch}_01_demultiplex.sh

done
