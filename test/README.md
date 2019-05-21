# Test data
This folder contains tiny test data for rapid testing of the pipeline. The folder structure is as follows:

**./test/contamination_Hsapiens_Mmusculus/**
These are files of 100 RNA-seq reads each, 80 of which come from human and 20 of which come from mouse.
They are paired end reads run on the Illumina machines. They come from different experiments, taken from
reads that map only to MT, so they have been adapter trimmed and quality controlled already. They are
from a RNA-seq dataset but are useful for testing contamination detection in any NGS data.
# Mitocondrial reads from WGS

**./test/wgs/**
The reads in `mt_1.fq.gz` and `mt_2.fq.gz` comes from the normal sample of the synthetic dataset 3 of the ICGC-TCGA DREAM challenge. Reads aligning to the mitochontrion were extracted from the aligned BAM file, and converted to fastq format. The dataset is 2x101 paired illumina reads, In total 4597 read pairs. Average coverage after removing duplicate is around 110x, but varies wildly. This provides suffcient depth to test variant calling downstream. 

**./test/hg19/**
This directory contains a subset of the human reference sequence hg19 (GRCh37), and the index files for bowtie2 and bwa as well as the dict files required to run GATK.

File Tree:
test
├── README.md
├── contamination_Hsapiens_Mmusculus
│   ├── Hsapiens_Mmusculus_1.fq.gz
│   ├── Hsapiens_Mmusculus_2.fq.gz
│   └── README.md
├── hg19
│   ├── bowtie2
│   │   ├── hg19.1.bt2
│   │   ├── hg19.2.bt2
│   │   ├── hg19.3.bt2
│   │   ├── hg19.4.bt2
│   │   ├── hg19.rev.1.bt2
│   │   └── hg19.rev.2.bt2
│   ├── bwa
│   │   ├── hg19.fa.amb
│   │   ├── hg19.fa.ann
│   │   ├── hg19.fa.bwt
│   │   ├── hg19.fa.pac
│   │   └── hg19.fa.sa
│   ├── seq
│   │   ├── hg19-resources.yaml
│   │   ├── hg19.dict
│   │   ├── hg19.fa
│   │   └── hg19.fa.fai
│   └── ucsc
│       └── hg19.2bit
└── wgs
    ├── README.md
    ├── mt_1.fq.gz
    └── mt_2.fq.gz
