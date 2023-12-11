#!/bin/bash
workdir=/date/guoyuqing/DNAloop/project/000_Merge_analysis/0420_6_celllines_H3K4me3_batch1_2/co_accessibility
cd $workdir
scripts=$workdir/co_occupancy.step2.sh

antibody='H3K4me3'

for cellline in K562 GM12878 HG002 293T H9 HFF1
do

    run_log=$workdir/run_log/$cellline
    mkdir -p $run_log

    echo "#!/bin/bash
    sh $scripts  $cellline $antibody 
    " > $run_log/${cellline}_co_accessibility.step2.sh

    chmod +x $run_log/${cellline}_co_accessibility.step2.sh


    qsub -o $run_log/${cellline}_co_accessibility.step2.o  \
    -e $run_log/${cellline}_co_accessibility.step2.e -cwd -l vf=20G -V $run_log/${cellline}_co_accessibility.step2.sh


done


