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
│   └── 07_repetitive_elements
│       ├── human_L1HS
│       │   ├── build_blast_diff_matrix_minus.sh
│       │   ├── build_blast_diff_matrix.plus.sh
│       │   ├── L1HS_analysis_pip.sh
│       │   ├── motif_ref_query_diff.new.R
│       │   └── T2T_fulllengthL1HS.trans.chrm13v2.0.txt
│       ├── human_L1HS_judgeable
│       │   ├── NGS_WGS
│       │   │   ├── get_bin_bulk.sh
│       │   │   ├── GIAB_12878_WGS_Unique_mapping_count.sh
│       │   │   └── S01_pkubatch_GIAB_12878_WGS_Unique_mapping.sh
│       │   └── TGS_scNanoSeq-CUT_Tag
│       │       └── get_L1HS_300bp_bin.sh
│       ├── mouse_fill-length_LINE
│       │   ├── L1Md_A
│       │   │   ├── AY053456.fasta
│       │   │   ├── build_blast_diff_matrix_minus.sh
│       │   │   ├── build_blast_diff_matrix.plus.sh
│       │   │   ├── motif_ref_query_diff.new.R
│       │   │   └── mouse_mm10_fulllength_LINE1_analysis_pip_L1Md_A.sh
│       │   └── L1Md_T
│       │       ├── AF016099.ORF1_ORF2.fasta
│       │       ├── build_blast_diff_matrix_minus.sh
│       │       ├── build_blast_diff_matrix.plus.sh
│       │       ├── motif_ref_query_diff.new.R
│       │       └── mouse_mm10_fulllength_LINE1_analysis_pip_L1Md_T.sh
│       └── mouse_fill-length_LINE_judgeable
│           ├── NGS_WGS
│           │   └── get_MGP_data.sh
│           └── TGS_scNanoSeq-CUT_Tag
│               └── get_cuttag_data.sh
└── WGBS
    ├── qsub_WGBS_pip.sh
    └── WGBS_pip.sh
```
# scNanoSeq-CUT&Tag pipeline
## 01_preprocess
These scripts in the 01_preprocess folder were used to demultiplex raw sequence data into single cell files according to the single-cell barcodes sequence and the single cell reads were next mapped to reference genome of human (hg38) and mouse (mm10).
### Step1: demultiplex

```
sh S01_pkubatch_01_demultiplex.sh
```
```
#samdir example
absolute_directiry_of_pass.fastq.gz    pass.fastq.gz
```
pass.fastq.gz is the raw sequencing data from Nanopore platform



