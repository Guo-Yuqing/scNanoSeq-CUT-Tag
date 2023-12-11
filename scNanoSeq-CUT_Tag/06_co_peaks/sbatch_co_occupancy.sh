#!/bin/bash
workdir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/guoyuqing/DNAloop/project/000_Merge_analysis/0420_6_celllines_H3K4me3_batch1_2/co_accessibility
cd $workdir
scripts=$workdir/co_occupancy.step1.sh

antibody='H3K4me3'

for cellline in K562 GM12878 HG002 293T H9 HFF1
do 
	for chr_ in chr{1..22} chrX
	do 

    run_log=$workdir/run_log/$cellline
    mkdir -p $run_log

    echo "#!/bin/bash
    sh $scripts  $cellline $antibody $chr_
    " > $run_log/${cellline}_${chr_}_co_accessibility.sh

    chmod +x $run_log/${cellline}_${chr_}_co_accessibility.sh

    pkubatch \
        -J GYQ_${cellline}_${antibody}_coa \
        -p cn-long \
        -A tangfuchou_g1 \
        -c 5  \
        --qos=tangfuchoucnl \
        -o $run_log/${cellline}_${chr_}_co_accessibility.o -e $run_log/${cellline}_${chr_}_co_accessibility.e \
        $run_log/${cellline}_${chr_}_co_accessibility.sh
	done
done

