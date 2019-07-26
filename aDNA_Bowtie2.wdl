## # WDL Workflow for aDNA processing

## Created by: Yassine Souilmi <yassine.souilmi@adelaide.edu.au>
## Last edited: Tuesday, 21 May 2019

## # Runs:
##  - FASTP quality control and filtering on input FastQ files
##  - Bowtie2 to map paired and collapsed reads
##  - Mark duplicates using samblaster
##  - Sort with samtools
##  - Generate mapping stats using samtools idxstats

##########################################################################
## Workflow

workflow AncientDNA_bowtie2 {
  # Indexes
  String ref_fasta
  String ref_fasta_basename
  String ref_dict

  # Data
  File samplesInfoTSV
  Array[Array[String]] inputSamples = read_tsv(samplesInfoTSV)

  scatter(sampleRow in inputSamples) {

    call FASTP {
      input:
        experimentName = sampleRow[0],
        sampleName = sampleRow[1],
        libraryName = sampleRow[2],
        runName = sampleRow[5],
        fastqRead1 = sampleRow[6],
        fastqRead2 = sampleRow[7]
    }

    call Bowtie2Collapsed {
        input:
        experimentName = sampleRow[0],
        sampleName = sampleRow[1],
        libraryName = sampleRow[2],
        platformName = sampleRow[3],
        unitName = sampleRow[4],
        runName = sampleRow[5],
        ref_fasta_basename = ref_fasta_basename,
        fastqFilteredCollapsed = FASTP.fastqFilteredCollapsed
    }

    call SamtoolsIdxstats {
        input:
        collapsed_mapped_markdup_bam =  Bowtie2Collapsed.collapsed_mapped_markdup_bam,
        experimentName = sampleRow[0],
        ref_fasta_basename = ref_fasta_basename,
        sampleName = sampleRow[1],
        libraryName = sampleRow[2],
        runName = sampleRow[5]
    }
  } #scatter  
} #workflow

##########################################################################
## Tasks

## Fastp statistics on fastq file and filtering
task FASTP {
  String sampleName
  String experimentName
  String runName
  String libraryName
  File fastqRead1
  File fastqRead2
  Int cores=4

  command {
    fastp \
    --thread ${cores} \
    -g -x -y -p -V -c \
    --in1 ${fastqRead1} \
    --in2 ${fastqRead2} \
    -R '${sampleName} – ${experimentName}' \
    -h ${sampleName}_${experimentName}_${libraryName}_${runName}_fastp.html \
    -j ${sampleName}_${experimentName}_${libraryName}_${runName}_fastp.json \
    --merge \
    --merged_out ${sampleName}_${experimentName}_${libraryName}_${runName}_collapsed_fastp.fastq.gz \
    --out1 ${sampleName}_${experimentName}_${libraryName}_${runName}_fastp_R1.fastq.gz \
    --out2 ${sampleName}_${experimentName}_${libraryName}_${runName}_fastp_R2.fastq.gz
  }

  output {
    File reportHTML = "${sampleName}_${experimentName}_${libraryName}_${runName}_fastp.html"
    File reportJSON = "${sampleName}_${experimentName}_${libraryName}_${runName}_fastp.json"
    File fastqFilteredCollapsed = "${sampleName}_${experimentName}_${libraryName}_${runName}_collapsed_fastp.fastq.gz"
    File fastqFilteredRead1 = "${sampleName}_${experimentName}_${libraryName}_${runName}_fastp_R1.fastq.gz"
    File fastqFilteredRead2 = "${sampleName}_${experimentName}_${libraryName}_${runName}_fastp_R2.fastq.gz"
  }

  runtime {
    cores: cores
  }
}

## BWA Mapping Collapsed reads
task Bowtie2Collapsed {
  File fastqFilteredCollapsed
  String ref_fasta_basename
  String sampleName
  String experimentName
  String runName
  String libraryName
  String platformName
  String unitName
  Int cores=4

  command {
    set -o pipefail
    set -e
    
    Bowtie2 \
        -x ${ref_fasta_basename} \
        --very-sensitive \
        --threads ${cores} \
        -U ${fastqFilteredCollapsed} \
        --rg-id ${sampleName}_collapsed_${experimentName} \
        --rg SM:${sampleName} \
        --rg LB:${libraryName} \
        --rg PL:${platformName} \
        --rg PU:${unitName} \
        | samblaster \
        | samtools sort \
            -@ ${cores} \
            -O BAM  \
            -o ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted.bam \
            - 
    samtools index ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted.bam
  }

  output {
    File collapsed_mapped_markdup_bam = "${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup.bam"
  }

  runtime {
    cores: cores
  }
}

## IdxStats on bam file
task SamtoolsIdxstats {
  File collapsed_mapped_markdup_bam
  String ref_fasta_basename
  String sampleName
  String experimentName
  String runName
  String libraryName

  command {
    samtools idxstats ${collapsed_mapped_markdup_bam} \
          > ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted.statis 
  }

  output {
    File collapsed_mapped_markdup_stats="${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted.statis"
  }
}