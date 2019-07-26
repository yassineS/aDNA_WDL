## # WDL Workflow for aDNA processing

## Created by: Yassine Souilmi <yassine.souilmi@adelaide.edu.au>
## Last edited: Tuesday, 21 May 2019

## # Runs:
##  - FASTP quality control and filtering on input FastQ files
##  - Bowtie2 to map paired and collapsed reads

##########################################################################
## Workflow

workflow AncientDNA_bowtie2 {
  # Indexes
  String ref_fasta
  String ref_fasta_basename
  String ref_fasta_index
  String ref_dict
  String ref_bowtie2_1
  String ref_bowtie2_2
  String ref_bowtie2_3
  String ref_bowtie2_4
  String ref_bowtie2_rev_1
  String ref_bowtie2_rev_2

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
        ref_fasta = ref_fasta,
        ref_fasta_basename = ref_fasta_basename,
        ref_bowtie2_1 = ref_bowtie2_1,
        ref_bowtie2_2 = ref_bowtie2_2,
        ref_bowtie2_3 = ref_bowtie2_3,
        ref_bowtie2_4 = ref_bowtie2_4,
        ref_bowtie2_rev_1 = ref_bowtie2_rev_1,
        ref_bowtie2_rev_2 = ref_bowtie2_rev_2,
        fastqFilteredCollapsed = FASTP.fastqFilteredCollapsed
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
  File ref_fasta_basename
  File ref_fasta
  File ref_bowtie2_1
  File ref_bowtie2_2
  File ref_bowtie2_3
  File ref_bowtie2_4
  File ref_bowtie2_rev_1
  File ref_bowtie2_rev_2
  String sampleName
  String experimentName
  String runName
  String libraryName
  String platformName
  String unitName
  Int cores = 4

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
  }

  output {
    File collapsed_mapped_bam = "${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted.bam"
  }

  runtime {
    cores: cores
  }
}


## Merge the two BAM files (Collapsed & PE)
task MergeBamAlignment {
  Array[File]+ inputCollapsedBams
  Array[File]+ inputPeBams

  command {
    samtools merge \
      -@ ${cores} \
      -s 100000 \
      -b \
      -p \
      -O BAM \
      -o ${sampleName}_${experimentName}_MEM_Merged.bam \
      ${sep=" " inputCollapsedBams} ${sep=" " inputPeBams}
  }
  runtime {
    cores: cores
  }
  output {
    File output_bam = "${sampleName}_${experimentName}_MEM_Merged.bam"
  }
}