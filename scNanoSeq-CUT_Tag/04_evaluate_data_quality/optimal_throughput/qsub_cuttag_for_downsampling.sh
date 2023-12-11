#!/bin/bash
workdir=`pwd`
cd $workdir
scripts=$workdir/cuttag_for_downsampling.sh
mkdir run_log


for antibody in H3K4me3
do 
    for rep in   `seq 1 15`
    do 
        for throughput in  `seq 1 10`
        do 
            log=./run_log/${antibody}_rep${rep}_throughput${throughput}K
            mkdir -p $log
            echo "#!/bin/bash
            sh $scripts $antibody $rep $throughput
            " >$log/${antibody}_rep${rep}_throughput${throughput}K.sh
            cd $log
           qsub -o ${antibody}_rep${rep}_throughput${throughput}K.o  -e ${antibody}_rep${rep}_throughput${throughput}K.e -cwd -l vf=10G -V ${antibody}_rep${rep}_throughput${throughput}K.sh
            cd $workdir
        done
    done 
done

