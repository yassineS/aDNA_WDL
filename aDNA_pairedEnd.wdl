## # WDL Workflow for aDNA processing

## Created by: Yassine Souilmi <yassine.souilmi@adelaide.edu.au>
## Last edited: Tuesday, 21 May 2019

## # Runs:
##  - FASTP quality control and filtering on input FastQ files


workflow StandardRNAseqPaired {
  # Annotation and indexes
  String annotationFile
  String STARIndexDir
  Int fastQCKmer

  # Data
  File samplesInfoTSV
  Array[Array[String]] inputSamples = read_tsv(samplesInfoTSV)

  scatter(sampleRow in inputSamples) {

    call FastQC as FastQCRaw {
      input:
        fastqRead1 = sampleRow[6],
        fastqRead2 = sampleRow[8],
        fastQCKmer = fastQCKmer
    }
  }
}

## Tasks

## Fastp statistics on fastq file and filtering
task FastQC {
  String sampleName
  String experimentName
  File fastqRead1
  File fastqRead2
  Int fastQCKmer
  Int cores=4

  command {
    fastp \
    --thread ${cores} \
    -g -x -y -p -V \
    --in1 ${fastqRead1} \
    --in2 ${fastqRead2} \
    -R '${s} – AHP High-Coverage' \
    -h ${s}_paired_fastp_report.html \
    -j ${s}_paired_fastp_report.json \
    --out1 ${s}_paired_fastp_R1.fq.gz \
    --out2 ${s}_paired_fastp_R2.fq.gz
    
    fastqc -t ${cores} \
      --outdir "." \
      -k ${fastQCKmer} \
      --noextract \
      -f fastq ${fastqRead1} ${fastqRead2}
  }

  output {
    File reportR1 = sub(basename(fastqRead1), "\\..+", "_fastqc.html")
    File reportZipR1 = sub(basename(fastqRead1), "\\..+", "_fastqc.zip")
    File reportR2 = sub(basename(fastqRead2), "\\..+", "_fastqc.html")
    File reportZipR2 = sub(basename(fastqRead2), "\\..+", "_fastqc.zip")
  }

  runtime {
    cores: cores
  }
}
