## # WDL Workflow for aDNA processing

## Created by: Yassine Souilmi <yassine.souilmi@adelaide.edu.au>
## Last edited: Tuesday, 21 May 2019

## # Runs:
##  - FASTP quality control and filtering on input FastQ files


workflow AncientDNA {
  # Annotation and indexes
  ##String ReferenceFile
  ##String BWAIndexDir

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
  }
}

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

## BWA Mapping
task BwaMemCollapsed {
  File fastqFilteredCollapsed
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit), 
  # listing the reference contigs that are "alternative". Leave blank in JSON for legacy 
  # references such as b37 and hg19.
  ###File? ref_alt

  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa

  command {
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}

  	java -Xmx3000m -jar /usr/gitc/picard.jar \
    	SamToFastq \
    	INPUT=${input_bam} \
    	FASTQ=/dev/stdout \
    	INTERLEAVE=true \
    	NON_PF=true | \
  	/usr/gitc/${bwa_commandline} /dev/stdin -  2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) | \
  	samtools view -1 - > ${output_bam_basename}.bam

  }
  runtime {
    cores: cores
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
  }
}
