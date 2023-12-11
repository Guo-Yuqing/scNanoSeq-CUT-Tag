#!/bin/bash
workdir=`pwd`
cd $workdir
scripts=$workdir/WGBS_pip.sh

cat $workdir/sample.list|grep -v "^#"|while read sample
do

    run_log=$workdir/run_log/
    mkdir -p $run_log

    echo "#!/bin/bash
    sh $scripts  $sample 
    " > $run_log/${sample}_WGBS_pip.sh

    chmod +x $run_log/${sample}_WGBS_pip.sh

    pkubatch \
        -J GYQ_${sample}_WGBS \
        -p cn-long \
        -A tangfuchou_g1 \
        -c 20  \
        --qos=tangfuchoucnl \
        -o $run_log/${sample}_WGBS_pip.o -e $run_log/${sample}_WGBS_pip.e \
        $run_log/${sample}_WGBS_pip.sh

done


