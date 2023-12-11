#!/bin/bash
cell=$1
antibody=$2

prefix=${cell}_${antibody}
rootdir=~
workdir=$rootdir/project/co_occupancy
ks_test_scripts=$workdir/co_occupancy_ks_test.R
summary_bedpe=$workdir/summary_coa.R
R_run=$rootdir/software/miniconda3/envs/r403/bin/Rscript

cd $workdir/${cell}/

mkdir -p $workdir/${cell}/coa_ks_test
mkdir -p $workdir/${cell}/coa_bedpe

do_merge(){

    less $rootdir/project/03_all_peak/scNanoSeqCUTTag/${prefix}_sc_seacr_top0.05.peaks.stringent.*cell_support.bed| \
    cut -f 1,2,3 >$workdir/${cell}/coa_ks_test/${prefix}_peak.bed

}

do_KS_test(){
    cd $workdir/${cell}/coa_ks_test

    for chr_ in chr{1..22} chrX
    do 
        $R_run $ks_test_scripts $cell $antibody $chr_
    done 
}

do_summary(){
    cd $workdir/${cell}/coa_bedpe
    $R_run $summary_bedpe $cell $antibody
}


run(){
    do_merge
    do_KS_test
    do_summary
}

run 

