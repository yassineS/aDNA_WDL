## # WDL Workflow for aDNA processing

## Created by: Yassine Souilmi <yassine.souilmi@adelaide.edu.au>
## Last edited: Tuesday, 21 May 2019

## # Runs:
##  - FASTP quality control and filtering on input FastQ files
##  - BWA MEM to map paired and collapsed reads

##########################################################################
## Workflow

workflow AncientDNA {
  # Indexes
  String ref_fasta
  String ref_fasta_index
  String ref_dict
  String ref_bwt
  String ref_amb
  String ref_ann
  String ref_pac
  String ref_sa

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

    call BwaMemCollapsed {
        input:
        experimentName = sampleRow[0],
        sampleName = sampleRow[1],
        libraryName = sampleRow[2],
        platformName = sampleRow[3],
        unitName = sampleRow[4],
        runName = sampleRow[5],
        ref_fasta = ref_fasta,
        ref_bwt = ref_bwt,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        fastqFilteredCollapsed = FASTP.fastqFilteredCollapsed
    }
    
    call BwaMemPE {
        input:
        experimentName = sampleRow[0],
        sampleName = sampleRow[1],
        libraryName = sampleRow[2],
        platformName = sampleRow[3],
        unitName = sampleRow[4],
        runName = sampleRow[5],
        ref_fasta = ref_fasta,
        ref_bwt = ref_bwt,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        fastqFilteredRead1 = FASTP.fastqFilteredRead1,
        fastqFilteredRead2 = FASTP.fastqFilteredRead2
    }
  } #scatter on input files
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
task BwaMemCollapsed {
  File fastqFilteredCollapsed
  File ref_fasta
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  String sampleName
  String experimentName
  String runName
  String libraryName
  String platformName
  String unitName
  Int cores = 4

  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit), 
  # listing the reference contigs that are "alternative". Leave blank in JSON for legacy 
  # references such as b37 and hg19.
  ###File? ref_alt


  command {
    set -o pipefail
    set -e
    
    bwa mem \
        -K 100000000 \
        -v 3 \
        -M \
        -t ${cores} \
        -R "@RG\tID:${sampleName}_collapsed_${experimentName}\tLB:${libraryName}\tPL:${platformName}\tPU:${unitName}\tSM:${sampleName}" \
        ${ref_fasta} \
        ${fastqFilteredCollapsed} | \
        samtools view -1 - > ${sampleName}_${experimentName}_${libraryName}_${runName}_collapsed_MEM.bam
  }

  output {
    File collapsed_mapped_bam = "${sampleName}_${experimentName}_${libraryName}_${runName}_collapsed_MEM.bam"
  }

  runtime {
    cores: cores
  }
}

## BWA Mapping PE reads (non collapsed)
task BwaMemPE {
  File fastqFilteredRead1
  File fastqFilteredRead2
  File ref_fasta
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  String sampleName
  String experimentName
  String runName
  String libraryName
  String platformName
  String unitName
  Int cores = 4

  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit), 
  # listing the reference contigs that are "alternative". Leave blank in JSON for legacy 
  # references such as b37 and hg19.
  ###File? ref_alt


  command {
    set -o pipefail
    #set -e
    
    bwa mem \
        -K 100000000 \
        -v 3 \
        -M \
        -t ${cores} \
        -R "@RG\tID:${sampleName}_PE_${experimentName}\tLB:${libraryName}\tPL:${platformName}\tPU:${unitName}\tSM:${sampleName}" \
        ${ref_fasta} \
        ${fastqFilteredRead1} ${fastqFilteredRead2} | \
        samtools view -1 - > ${sampleName}_${experimentName}_${libraryName}_${runName}_PE_MEM.bam
  }

  output {
    File collapsed_mapped_bam = "${sampleName}_${experimentName}_${libraryName}_${runName}_PE_MEM.bam"
  }

  runtime {
    cores: cores
  }
}

## Get version of BWA
task GetBwaVersion {
  String docker
  command {
    # Not setting "set -o pipefail" here because /bwa has a rc=1 and we don't want to allow rc=1 to succeed 
    # because the sed may also fail with that error and that is something we actually want to fail on.
    /usr/gitc/bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
  runtime {
    memory: "1 GB"
  }
  output {
    String version = read_string(stdout())
  }
}

## Merge the two BAM files (Collapsed & PE)
task MergeBamAlignment {
  File collapsed_mapped_bam
  File pe_mapped_bam
  String bwa_version
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  command {
    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    bwa_commandline="bwa mem -K 100000000 -v 3 -M"

    java -Xmx2500m -jar /usr/gitc/picard.jar \
      MergeBamAlignment \
      VALIDATION_STRINGENCY=SILENT \
      EXPECTED_ORIENTATIONS=FR \
      ATTRIBUTES_TO_RETAIN=X0 \
      ALIGNED_BAM=${aligned_bam} \
      UNMAPPED_BAM=${unmapped_bam} \
      OUTPUT=${output_bam_basename}.bam \
      REFERENCE_SEQUENCE=${ref_fasta} \
      PAIRED_RUN=true \
      SORT_ORDER="unsorted" \
      IS_BISULFITE_SEQUENCE=false \
      ALIGNED_READS_ONLY=false \
      CLIP_ADAPTERS=false \
      MAX_RECORDS_IN_RAM=2000000 \
      ADD_MATE_CIGAR=true \
      MAX_INSERTIONS_OR_DELETIONS=-1 \
      PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
      PROGRAM_RECORD_ID="bwamem" \
      PROGRAM_GROUP_VERSION="${bwa_version}" \
      PROGRAM_GROUP_COMMAND_LINE="${bwa_commandline}" \
      PROGRAM_GROUP_NAME="bwamem" \
      UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
      ALIGNER_PROPER_PAIR_FLAGS=true \
      UNMAP_CONTAMINANT_READS=true
  }
  runtime {
    memory: "3500 MB"
    cpu: "1"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
  }
}