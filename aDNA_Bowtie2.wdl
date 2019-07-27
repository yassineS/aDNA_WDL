## # WDL Workflow for aDNA processing

## Created by: Yassine Souilmi <yassine.souilmi@adelaide.edu.au>
## Last edited: Tuesday, 21 May 2019

## # Runs:
##  - FASTP quality control and filtering on input FastQ files
##  - Bowtie2 to map paired and collapsed reads
##  - Mark duplicates using samblaster
##  - Sort with samtools
##  - Generate mapping stats using samtools idxstats
##  - Estimate complexity, coverage and sequencing yield extrapolation
##  - Indel realignement with GATK
##  - Postmortem damage and base requalibration with mapDamage
##  - Softclipping 3bp from both sides of each mapped read
##  - Estimating coverage and QC mapped bams (postprocessed bams)

##########################################################################
## Workflow

workflow AncientDNA_bowtie2 {
  # Indexes
  File ref_fasta
  String ref_fasta_basename
  File gatk3_jar

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

    call preseq {
        input:
          collapsed_mapped_markdup_bam =  Bowtie2Collapsed.collapsed_mapped_markdup_bam,
          experimentName = sampleRow[0],
          ref_fasta_basename = ref_fasta_basename,
          sampleName = sampleRow[1],
          libraryName = sampleRow[2],
          runName = sampleRow[5]
    }

    call IndelRealignment {
      input:
          collapsed_mapped_markdup_bam =  Bowtie2Collapsed.collapsed_mapped_markdup_bam,
          experimentName = sampleRow[0],
          ref_fasta_basename = ref_fasta_basename,
          sampleName = sampleRow[1],
          libraryName = sampleRow[2],
          runName = sampleRow[5],
          gatk3_jar = gatk3_jar,
          ref_fasta = ref_fasta

    }

    call mapDamage {
      input:
          collapsed_mapped_markdup_IndelReal_bam = IndelRealignment.collapsed_mapped_markdup_IndelReal_bam,
          experimentName = sampleRow[0],
          ref_fasta_basename = ref_fasta_basename,
          sampleName = sampleRow[1],
          libraryName = sampleRow[2],
          runName = sampleRow[5],
          ref_fasta = ref_fasta
    } 

    call trimBam {
      input:
          rescaled_bam = mapDamage.rescaled_bam,
          experimentName = sampleRow[0],
          ref_fasta_basename = ref_fasta_basename,
          sampleName = sampleRow[1],
          libraryName = sampleRow[2],
          runName = sampleRow[5]
    }

    call Qualimap {
      input:
          trimmed_bam = trimBam.trimmed_bam,
          experimentName = sampleRow[0],
          ref_fasta_basename = ref_fasta_basename,
          sampleName = sampleRow[1],
          libraryName = sampleRow[2],
          runName = sampleRow[5]
    }

    call FreeBayes {
      input
    }
    ## SequenceTools
    ## Shmutzi
    ## MultiQC
    ## Stats
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
    cpu: cores
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
    cpu: cores
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
          > ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted.stats 
  }

  output {
    File collapsed_mapped_markdup_stats="${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted.stats"
  }
}

## Preseq
task preseq {
  File collapsed_mapped_markdup_bam
  String ref_fasta_basename
  String sampleName
  String experimentName
  String runName
  String libraryName

  command {
    preseq c_curve \
        -seed 1234 \
        -bam \
        -output ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted.ComplexityCurve.txt
        ${collapsed_mapped_markdup_bam}

    preseq lc_extrap \
        -seed 1234 \
        -bam \
        -output ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted.YieldCurve.txt
        ${collapsed_mapped_markdup_bam}

    preseq gc_extrap \
        -seed 1234 \
        -output ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted.CoverageCurve.txt \
        ${collapsed_mapped_markdup_bam}
  }
  output {
    File complexityCurve = "${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted.ComplexityCurve.txt"
    File yieldCurve = "${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted.YieldCurve.txt"
    File coverageCurve = "${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted.CoverageCurve.txt"
  }
}

## IndelRealignment
task IndelRealignment  {
  File collapsed_mapped_markdup_bam
  File ref_fasta
  File gatk3_jar
  String sampleName
  String experimentName
  String runName
  String libraryName
  String ref_fasta_basename
  Int cores=4
  String mem='12G'

  command {
    java -jar ${gatk3_jar} \
        -Xmx${mem} \
        -T RealignerTargetCreator \
        -R ${ref_fasta} \
        --num_threads ${cores} \
        --mismatchFraction 0.30 \
        --maxIntervalSize 650 \
        --allow_potentially_misencoded_quality_scores \
        -I ${collapsed_mapped_markdup_bam} \
        -o ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted_IndelReal.intervals 


    java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
        -Xmx${mem} \
        -T IndelRealigner \
        -R ${ref_fasta} \
        -model USE_READS \
        -compress 0 \
        --filter_bases_not_stored \
        --allow_potentially_misencoded_quality_scores \
        -I ${collapsed_mapped_markdup_bam} \
        -targetIntervals ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted_IndelReal.intervals \
        -o ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted_IndelReal.bam
  }

  output {
    File collapsed_mapped_markdup_IndelReal_bam = "${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted_IndelReal.bam"
  }

  runtime {
    cpu: cores
    memory: mem
  }
}

## mapDamage
task mapDamage {
  File collapsed_mapped_markdup_IndelReal_bam
  File ref_fasta
  String sampleName
  String experimentName
  String runName
  String libraryName
  String ref_fasta_basename

  command {
    mapDamage \
        -i ${collapsed_mapped_markdup_IndelReal_bam}  \
        -r ${ref_fasta} \
        -m 25 \
        -d ./ \
        --title='${sampleName}_${experimentName}_${libraryName}_${runName}' \
        --rescale \
        --rescale-out=${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted_IndelReal.mapDamage.bam

    samtools index ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted_IndelReal.mapDamage.bam
  }

  output {
    File rescaled_bam = "${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted_IndelReal.mapDamage.bam"
    File smiley_plot = "Fragmisincorporation_plot.pdf"
    File length_plot = "Length_plot.pdf"
  }
}

## BamUtil trimbam
task trimBam {
  File rescaled_bam
  String sampleName
  String experimentName
  String runName
  String libraryName
  String ref_fasta_basename
  Int cores=4

  command {
    bam trimBam \
        --ignoreStrand \
        --clip \
        ${rescaled_bam} \
        ${sampleName}_${experimentName}_${libraryName}_${runName}_tmp.bam \
        3

    samtools sort -@ ${cores} \
        -O BAM \
        -l 7 \
        -o ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted_IndelReal.mapDamage.trim3_2ends.bam \
        ${sampleName}_${experimentName}_${libraryName}_${runName}_tmp.bam

	  samtools index ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted_IndelReal.mapDamage.trim3_2ends.bam
  }
  output {
    File trimmed_bam = "${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted_IndelReal.mapDamage.trim3_2ends.bam"
  }

  runtime {
    cpu: cores
  }
}

## Qualimap
task Qualimap {
  File rescaled_bam
  String sampleName
  String experimentName
  String runName
  String libraryName
  String ref_fasta_basename
  Int cores=4
  Int mem=8

  command {
    qualimap bamqc \
        -c \
        -nt ${cores} \
        --skip-duplicated --skip-dup-mode 0 \
        --java-mem-size='${mem}G' \
        -outdir ./ \    
        -bam ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bowtie2_markdup_sorted_IndelReal.mapDamage.trim3_2ends.bam
  }
  output {
    File qualimap_html = "qualimapReport.html"
    File qualimap_report_txt = "genome_results.txt"
  }

  runtime {
    cpu: cores
    memory: mem + "GB"
  }
}
## FreeBayes
## SequenceTools
## Shmutzi
## MultiQC
## Stats