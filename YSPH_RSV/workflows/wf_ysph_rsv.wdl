version 1.0

import "../tasks/task_bwa.wdl" as bwa_task
import "../tasks/task_bcftools.wdl" as bcftools_task
import "../tasks/task_ivar.wdl" as ivar_task
import "../tasks/task_mash.wdl" as mash_task
import "../tasks/task_versioning.wdl" as versioning

## rsveillance

workflow ysph_rsv {
  input {
    Array[File] read1
    Array[File] read2
    Array[String] samplename

    # reference files -- will he know in advance if rsva or rsvb?
    File reference_fasta
    File reference_mash_index
    File reference_bwa_index
    File reference_gff
    File primer_bed
    File amplicon_bed
    
    # questions for seth: 
    #  does he still want it to run on the set or is individual okay for terra?
    #  does he want commented out tasks implemented?
    #  do you know if it will be rsva or rsvb in advance?
    #    would you prefer to have a boolean option to indicate either rsva or rsvb?
    #  do you want to save the indexes for the references instead of regenerating them each time? (mash/bwa)
    #  mash questions:
    #   is there a reason why you're only mashing the top 10000 reads instead of all of them?
    #   what is the goal of the four functions at the end of mash.smk (line 88 - 115) -- do you have example output and is this still something you want?
    #  bwa questions:
    #   why are you subsampling the bam files? is this mainly a size requirement?
    #   this file is used in the get_depth_distribution.py script -- what is the reasoning?
    #  variant calling quesitons:
    #   do you want to continue calling variants on the untrimmed bam files as well?
    #  stats questions:
    #   do you want all of these outputs or are there a few that are most important?
    #   what is the output format of these files?
    # other questions:
    #  do you have example input/output files for the workflow?
    #  anything else?
    

    ### references: RSVA.fasta & RSVB.fasta & assoc. gffs
    ### amplicons/primers: RSVA & RSVB
    ### metadata: r script to pull from gisaid and google sheets -- individual files exist in repo

  }
  call versioning.version_capture {
    input:
  }
  call mash_task.mash {
    input:
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      reference_mash_index = reference_mash_index
  }
  call bwa_task.bwa {
    input:
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      reference_fasta = reference_fasta,
      reference_bwa_index = reference_bwa_index
  }
  call ivar_task.ivar_trim {
    input:
      aligned_bam = bwa.sorted_bam,
      samplename = samplename,
      reference_fasta = reference_fasta,
      primer_bed = primer_bed
  }
  call ivar_task.ivar_consensus {
    input:
      mpileup = ivar_trim.mpileup,
      samplename = samplename
  }
  call ivar_task.ivar_variants {
    input:
      mpileup = ivar_trim.mpileup,
      samplename = samplename,
      reference_fasta = reference_fasta
  }
  call bcftools_task.bcftools as trimmed_bcftools {
    input:
      input_bam = ivar_trim.trimmed_bam,
      input_bai = ivar_trim.input_bai,
      samplename = samplename,
      reference_fasta = reference_fasta
  }
  call bcftools_task.bcftools as untrimmed_bcftools {
    input:
      input_bam = bwa.sorted_bam,
      input_bai = bwa.sorted_bai,
      samplename = samplename,
      reference_fasta = reference_fasta
  }



  output {
    String ysph_rsv_analysis_date = version_capture.date
    String ysph_rsv_wf_version = version_capture.version


# alignstats.txt ->
# consensus.fasta -> from get_ivar_report -> input: ivar (ivariants.tsv & consensus.fa)
# genotype-summary.txt -> from get_genotyping_report -> input: bcftools (all.vcf.gz & untrim_all.vcf.gz) & primer.bed
# final_calls.txt ->
# final_alignstats.tzt -> from call_rsvab -> input: alignstats.txt
# ampdepths -> 

  }
}