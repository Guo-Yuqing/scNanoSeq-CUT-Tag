# scNanoSeq-CUT-Tag
scNanoSeq-CUT&Tag is a streamlined method by adapting a modified CUT&Tag protocol to Oxford Nanopore sequencing platform for efficient chromatin modification profiling at single-cell resolution. scNanoSeq-CUT&Tag can accurately distinguish different cell types in vitro and in vivo. Moreover, scNanoSeq-CUT&Tag enables to effectively map the allele-specific epigenomic modifications in the human genome and allows to analyze co-occupancy of histone modifications. Taking advantage of long-read sequencing, scNanoSeq-CUT&Tag can sensitively detect epigenomic state of repetitive elements. In addition, by applying scNanoSeq-CUT&Tag to testicular cells of adult mouse B6D2F1, we demonstrated that scNanoSeq-CUT&Tag maps dynamic epigenetic state changes during mouse spermatogenesis.

 ![image](https://github.com/Guo-Yuqing/scNanoSeq-CUT-Tag/assets/153424744/b3949a8d-72d8-48fa-acc9-995fce69463a)![image](https://github.com/Guo-Yuqing/scNanoSeq-CUT-Tag/assets/153424744/e6262d93-02e2-4554-a010-f76476fd7cae)
```
├── scNanoSeq-CUT_Tag
│   ├── 01_preprocess
│   │   ├── 01_demultiplex.sh
│   │   ├── 02_trim_mapping.sh
│   │   ├── barcode_inner.fa
│   │   ├── barcode_outer.fa
│   │   ├── S01_pkubatch_01_demultiplex.sh
│   │   └── S02_pkubatch_02_trim_mapping.sh
│   ├── 02_clustering
│   │   ├── 03_prearchr.sh
│   │   └── 04_archr.R
│   ├── 03_peak_calling
│   │   ├── get_fragment_callpeak.sh
│   │   └── sbatch_get_fragment_callpeak.sh
│   ├── 04_evaluate_data_quality
│   │   ├── cross_condamination
│   │   │   ├── 02_trim_mapping_hg38_mm10.sh
│   │   │   └── S02_pkubatch_02_trim_mapping.sh
│   │   ├── optimal_throughput
│   │   │   ├── archr_arrow_cluster.R
│   │   │   ├── cuttag_for_downsampling.sh
│   │   │   └── qsub_cuttag_for_downsampling.sh
│   │   └── reproducibility
│   │       └── get_correlation.sh
│   ├── 05_ASPs
│   │   ├── Identification
│   │   │   ├── ASP_identification.R
│   │   │   ├── pre_ASP_filter_reads.sh
│   │   │   ├── pull_reads_new.py
│   │   │   └── sbatch_pre_ASP_filter_reads.sh
│   │   └── Validation
│   │       └── ENCODE_validation.R
│   ├── 06_co_peaks
│   │   ├── co_occupancy_ks_test.R
│   │   ├── co_occupancy.step1.sh
│   │   ├── co_occupancy.step2.sh
│   │   ├── qsub_co_occupancy.step2.sh
│   │   ├── sbatch_co_occupancy.sh
│   │   └── summary_coa.R
│   ├── 07_repetitive_elements
│   │   ├── human_L1HS
│   │   │   ├── build_blast_diff_matrix_minus.sh
│   │   │   ├── build_blast_diff_matrix.plus.sh
│   │   │   ├── L1HS_analysis_pip.sh
│   │   │   ├── motif_ref_query_diff.new.R
│   │   │   └── T2T_fulllengthL1HS.trans.chrm13v2.0.txt
│   │   ├── human_L1HS_judgeable
│   │   │   ├── NGS_WGS
│   │   │   │   ├── get_bin_bulk.sh
│   │   │   │   ├── GIAB_12878_WGS_Unique_mapping_count.sh
│   │   │   │   └── S01_pkubatch_GIAB_12878_WGS_Unique_mapping.sh
│   │   │   └── TGS_scNanoSeq-CUT_Tag
│   │   │       └── get_L1HS_300bp_bin.sh
│   │   ├── mouse_fill-length_LINE
│   │   │   ├── L1Md_A
│   │   │   │   ├── AY053456.fasta
│   │   │   │   ├── build_blast_diff_matrix_minus.sh
│   │   │   │   ├── build_blast_diff_matrix.plus.sh
│   │   │   │   ├── motif_ref_query_diff.new.R
│   │   │   │   └── mouse_mm10_fulllength_LINE1_analysis_pip_L1Md_A.sh
│   │   │   └── L1Md_T
│   │   │       ├── AF016099.ORF1_ORF2.fasta
│   │   │       ├── build_blast_diff_matrix_minus.sh
│   │   │       ├── build_blast_diff_matrix.plus.sh
│   │   │       ├── motif_ref_query_diff.new.R
│   │   │       └── mouse_mm10_fulllength_LINE1_analysis_pip_L1Md_T.sh
│   │   └── mouse_fill-length_LINE_judgeable
│   │       ├── NGS_WGS
│   │       │   └── get_MGP_data.sh
│   │       └── TGS_scNanoSeq-CUT_Tag
│   │           └── get_cuttag_data.sh
│   └── 08_mouse_spermatogenesis
│       └── mouse_spermatogenesis_analysis.r
└── WGBS
    ├── qsub_WGBS_pip.sh
    └── WGBS_pip.sh
```
# Requirements
Please have the following softwares installed first:
```
###scNanoSeq-CUT&Tag
nanoplexer (v0.1) <https://github.com/hanyue36/nanoplexer>
NanoFilt (v2.8.0) 
cutadapt (v3.5)
minimap (v2.26-r1175)
samtools (v1.3)
bedtools (v2.30.0)
SEACR (v1.3) 
deepTools (v3.5.1)
whatshap (v1.6)
blast (v2.5.0)
GenMap (v1.3.0)
ArchR (v1.0.1)

(and all the packages declared at the head of R scripts...)

###NGS WGS
fastp (v0.23.0)
bwa mem (v0.7.17)
sambamba (v1.0.0)

###NGS WGBS
bismark (v0.23.1)

```

# scNanoSeq-CUT&Tag pipeline
## 01_preprocess
These scripts in the `01_preprocess` folder were used to demultiplex raw sequence data into single cell files according to the single-cell barcodes sequence and the single cell reads were next mapped to reference genome of human (hg38) and mouse (mm10).
### Step1: demultiplex
```
sh S01_pkubatch_01_demultiplex.sh
```
```
#samdir example
#absolute_directory_of_pass.fastq.gz    subfolder_containing_fastq_file
/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/guoyuqing/scnanoseqcuttag/20221223_LP_221208_C6/pass.fastq.gz 20221223_LP_221208_C6
```
pass.fastq.gz is the raw sequencing data from Nanopore platform
### Step2: trim, align, filter
```
sh S02_pkubatch_02_trim_mapping.sh
```
```
#samdir example
head GM12878_H3K4me3_T10_sample.list
#barcode_name directory new_name
Bc81_Bc01       20221223_LP_221208_C6   GM12878_H3K4me3_T10
Bc81_Bc02       20221223_LP_221208_C6   GM12878_H3K4me3_T10
Bc81_Bc03       20221223_LP_221208_C6   GM12878_H3K4me3_T10
Bc81_Bc04       20221223_LP_221208_C6   GM12878_H3K4me3_T10
Bc81_Bc05       20221223_LP_221208_C6   GM12878_H3K4me3_T10
```
## 02_clustering
These scripts in the `02_clustering` folder were used to analyse chromatin modification signal of scNanoSeq-CUT&Tag, including extraction of chromatin modification signal; clustering of scNanoSeq-CUT&Tag profiles; generation of genomic coverage track; identification of cell type mark genes.
### Step1: prepare files for archr
Extraction of chromatin modification signal:
```
sh 03_prearchr.sh $antibody
```
### Step2: Clustering of scNanoSeq-CUT&Tag profiles
For active chromatin marks (H3K4me3, H3K27ac, H3K36me3, CTCF, RAD21), the parameters were setting as 'minTSS = 1'. For repressive chromatin marks (H3K27me3, H3K9me3), the parameters were setting as 'minTSS = 0.1'. 
```
Rscript 04_archr.R $antibody
```
## 03_peak_calling
These scripts in the `03_peak_calling` folder were used to peak calling. For scNanoSeq-CUT&Tag has a high signal-to-noise ratio, peak calling was performed using SEACR (v1.3) with parameters 'non stringent 0.05'. To further improve the precision of peak calling, we retained the peaks supported by at least (1.5% * No. of total cells) cells, setting a minimum threshold of 5 cells for support. .
```
sh sbatch_get_fragment_callpeak.sh
```
## 04_evaluate_data_quality
### cross_condamination
These scripts in the `04_evaluate_data_quality/cross_condamination` folder were used to evaluate human-mouse cross-contamination libraries for scNanoSeq-CUT&Tag.
```
sh S02_pkubatch_02_trim_mapping.sh
```
```
head 20230331_GM12878_3T3_H3K4me3_sample.list
#barcode_name directory new_name
Bc77_Bc01       20230316_bulk_newRAD21_HumanMouse_cross_C1      GM12878_3T3_H3K4me3
Bc77_Bc02       20230316_bulk_newRAD21_HumanMouse_cross_C1      GM12878_3T3_H3K4me3
Bc77_Bc03       20230316_bulk_newRAD21_HumanMouse_cross_C1      GM12878_3T3_H3K4me3
Bc77_Bc04       20230316_bulk_newRAD21_HumanMouse_cross_C1      GM12878_3T3_H3K4me3
Bc77_Bc05       20230316_bulk_newRAD21_HumanMouse_cross_C1      GM12878_3T3_H3K4me3
```
### optimal_throughput
These scripts in the `04_evaluate_data_quality/optimal_throughput` folder were used to evaluate the optimal throughput of scNanoSeq-CUT&Tag by downsampling analysis.
```
sh qsub_cuttag_for_downsampling.sh
```
### reproducibility
These scripts in the `04_evaluate_data_quality/reproducibility` folder were used to evaluate the correlation among different cell lines with the same antibody and the correlation in one cell line with different antibodies.
```
sh get_correlation.sh $cellline
```
The bw files were generated using `getGroupBW` in Archr during step 02_clustering.
## 05_ASPs
### Identification of ASPs
These scripts in the `05_ASPs/Identification` folder were used to identify allele-specific chromatin modifications peak from scNanoSeq-CUT&Tag data.
* step1: Split pat mat reads
* step2: Subset each reads HetSNPs
* step3: Filter reads according to HetSNPs and duplicate reads
* step4: Merge all single cell filter reads
* step5: Identify allele-specific chromatin modifications peak using `binom.test`
```
sh sbatch_pre_ASP_filter_reads.sh
```
The input bam files were raw mapping reads without remove duplicate reads, since high error rate of third generation sequencing, we corrected the phasing result of each reads by HetSNPs and duplicate reads.
```
Rscript ASP_identification.R
```
### Validation of ASPs
These scripts in the `05_ASPs/Validation` folder were used to validate allele-specific chromatin modifications peak from scNanoSeq-CUT&Tag data using ChIP-seq data from the ENCODE database.
```
Rscript ENCODE_validation.R
```
## 06_co_peaks
These scripts in the `06_co_peaks` were used to detect chromatin modification co-occupacy peaks events from scNanoSeq-CUT&Tag data.
* step1: Find reads directly support peak pair as candidate co-occupancy peak (Such as peak1-peak2; peak1-peak3 )
* step2: Calculate reads strand and length distribution in each peak (Only save reads length >1 kb for calculation)
* step3: Filter candidate peak pairs:  
         1. left peak pvalue adj <0.01  
         2. Right peak pvalue adj <0.01   
         3. Filter peak length > 1 kb (if a peak is long ,the resolution is poor, it may has many regulatory elements)
```
sh sbatch_co_occupancy.sh
sh qsub_co_occupancy.step2.sh
```
## 07_repetitive_elements
These scripts in the `07_repetitive_elements` were used to compare and detect full-length LINE in human and mouse genome.  

*for human*  
Based on the third-generation sequencing technology, T2T human reference genome has been assembled and 320 full-length L1Hs have been identified. So for human, we used the T2T reference genome for the analysis of this section. 
```
sh L1HS_analysis_pip.sh
```
*for mouse*  
The full-length LINEs of mouse were collected from L1Base.
```
sh mouse_mm10_fulllength_LINE1_analysis_pip_L1Md_A.sh
sh mouse_mm10_fulllength_LINE1_analysis_pip_L1Md_T.sh
```
## 08_mouse_spermatogenesis
These scripts in the `08_mouse_spermatogenesis` were used to analyse the dynamic epigenetic state changes during mouse spermatogenesis.  
```
Rscript mouse_spermatogenesis_analysis.r
```

# WGBS pipeline
These scripts in the `WGBS` were used to analysis whole genome bisulfite sequencing.
```
sh qsub_WGBS_pip.sh
```
# Contact
guoyuqing@stu.pku.edu.cn



