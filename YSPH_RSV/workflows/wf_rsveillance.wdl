version 1.0

import "../tasks/task_bwa.wdl" as bwa_task
import "../tasks/task_bcftools.wdl" as bcftools_task
import "../tasks/task_ivar_consensus.wdl" as ivar_consensus_task
import "../tasks/task_ivar_trim.wdl" as ivar_trim_task
import "../tasks/task_ivar_variants.wdl" as ivar_variants_task
import "../tasks/task_mash.wdl" as mash_task
import "../tasks/task_mash_index.wdl" as mash_index_task
import "../tasks/task_versioning.wdl" as versioning

workflow rsveillance {
  input {
    File read1
    File read2
    String samplename

    String reference_prefixes # **space**-delimited list of the prefixes of the potential references
    String reference_files_gcuri

    File primer_bed
    File amplicon_bed
  }
  call versioning.version_capture {
    input:
  }
  call mash_index_task.mash_index {
    input:
      reference_prefixes = reference_prefixes,
      reference_files_gcuri = reference_files_gcuri
  }
  call mash_task.mash {
    input:
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      reference_prefixes = reference_prefixes,
      reference_files_gcuri = reference_files_gcuri,
      reference_mash_index = mash_index.mash_all_index
  }
  call bwa_task.bwa {
    input:
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      reference_fasta = mash.reference_fasta
  }
  call ivar_trim_task.ivar_trim {
    input:
      aligned_bam = bwa.sorted_bam,
      samplename = samplename,
      primer_bed = primer_bed
  }
  call ivar_consensus_task.ivar_consensus {
    input:
      trimmed_bam = ivar_trim.trimmed_bam,
      trimmed_bai = ivar_trim.trimmed_bai,
      samplename = samplename,
      reference_fasta = mash.reference_fasta,
  }
  call ivar_variants_task.ivar_variants {
    input:
      mpileup = ivar_trim.mpileup,
      samplename = samplename,
      reference_fasta = mash.reference_fasta
  }
  # to-do: add stats rules here



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