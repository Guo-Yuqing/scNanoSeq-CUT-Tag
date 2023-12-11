#!/bin/bash
workdir=`pwd`
cd $workdir
scripts=$workdir/get_fragment_callpeak.sh

antibody='H3K4me3'   #H3K27ac H3K36me3 H3K27me3 H3K9me3 CTCF RAD21

for cellline in K562  GM12878 HG002 293T H9 HFF1
do 
    run_log=$workdir/run_log/$cellline
    mkdir -p $run_log

    echo "#!/bin/bash
    sh $scripts  $cellline $antibody
    " > $run_log/${cellline}_get_fragment_callpeak.sh

    chmod +x $run_log/${cellline}_get_fragment_callpeak.sh

    pkubatch \
        -J GYQ_${cellline}_${antibody}_callpeak \
        -p cn-long \
        -A tangfuchou_g1 \
        -c 15  \
        --qos=tangfuchoucnl \
        -o $run_log/${cellline}_get_fragment_callpeak.o -e $run_log/${cellline}_get_fragment_callpeak.e \
        $run_log/${cellline}_get_fragment_callpeak.sh

done

