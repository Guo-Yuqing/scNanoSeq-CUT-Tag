#!/bin/sh
rootdir=~
workdir=$rootdir/project
cd $workdir

fastq=$1
batch=$2

barcode_dir=$workdir/barcode/$batch

barcode_inner=${barcode_dir}/barcode_inner.fa
barcode_outer=${barcode_dir}/barcode_outer.fa

nanoplexer=$rootdir/software/miniconda3/bin/nanoplexer
seqkit=$rootdir/software/miniconda3/bin/seqkit
demultiplex_dir=$workdir/01_demultiplex/$batch


mkdir -p $demultiplex_dir

# outer barcode
$nanoplexer -b ${barcode_outer} \
         -p $demultiplex_dir \
         -t 8 \
         ${fastq}

cd $demultiplex_dir
ls $demultiplex_dir | grep Barcode > file.list

# inner barcode
while read fastq
do
cell=${fastq%.*}
Bc=${cell#*e}
mkdir $cell
$nanoplexer -b ${barcode_inner} -p $cell -t 8 ${cell}.fastq
cd $cell
rename Barcode Bc${Bc}_Bc Barcode*
cd ..
done < file.list

mv */Bc*fastq ./
cat */unclassified.fastq >> unclassified.fastq

$seqkit stat *.fastq >stat

rm -rf Barcode*

