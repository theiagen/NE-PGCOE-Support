version 1.0

import "../tasks/task_bcftools.wdl" as bcftools_task
import "../tasks/task_bwa.wdl" as bwa_task
import "../tasks/task_ivar_consensus.wdl" as ivar_consensus_task
import "../tasks/task_ivar_trim.wdl" as ivar_trim_task
import "../tasks/task_ivar_variants.wdl" as ivar_variants_task
import "../tasks/task_mash.wdl" as mash_task
import "../tasks/task_mash_index.wdl" as mash_index_task
import "../tasks/task_stats.wdl" as stats
import "../tasks/task_utility.wdl" as utility
import "../tasks/task_versioning.wdl" as versioning

workflow rsveillance {
  input {
    Array[File] read1
    Array[File] read2
    Array[String] samplenames

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
  scatter (triplet in zip(zip(read1, read2), samplenames)) {
    File current_read1 = triplet.left.left
    File current_read2 = triplet.left.right
    String current_sample = triplet.right

    call mash_task.mash {
      input:
        read1 = current_read1,
        read2 = current_read2,
        samplename = current_sample,
        reference_prefixes = reference_prefixes,
        reference_files_gcuri = reference_files_gcuri,
        reference_mash_index = mash_index.mash_all_index
    }
    scatter (index in range(length(mash.reference_fastas))) {
      File first_current_reference_fasta = mash.reference_fastas[index]
      File first_current_reference_gff = mash.reference_gffs[index]
      String first_current_reference_name = mash.reference_names[index]

      call bwa_task.bwa {
        input:
          read1 = current_read1,
          read2 = current_read2,
          output_prefix = current_sample + "_" + first_current_reference_name,
          reference_fasta = first_current_reference_fasta,
      }
      call ivar_trim_task.ivar_trim {
        input:
          aligned_bam = bwa.sorted_bam,
          aligned_bai = bwa.sorted_bai,
          output_prefix = current_sample + "_" + first_current_reference_name,
          primer_bed = primer_bed
      }
      call ivar_consensus_task.ivar_consensus {
        input:
          trimmed_bam = ivar_trim.trimmed_bam,
          trimmed_bai = ivar_trim.trimmed_bai,
          output_prefix = current_sample + "_" + first_current_reference_name,
          reference_fasta = first_current_reference_fasta,
      }
      call ivar_variants_task.ivar_variants {
        input:
          mpileup = ivar_consensus.mpileup,
          output_prefix = current_sample + "_" + first_current_reference_name,
          reference_fasta = first_current_reference_fasta
      }
      call stats.get_depth_histogram {
        input:
          depth_file = bwa.depth_file,
          reference_name = first_current_reference_fasta,
          associated_sample = current_sample
      }
      call stats.get_alignment_stats {
        input:
          alignment_flagstat = bwa.flagstat,
          depth_histogram = get_depth_histogram.depth_histogram,
          reference_name = first_current_reference_fasta,
          associated_sample = current_sample
      }
      call stats.get_depths {
        input:
          depth_file = bwa.depth_file,
          amplicon_bed = amplicon_bed,
          reference_gff = first_current_reference_gff,
          reference_name = first_current_reference_name,
          associated_sample = current_sample
      }
      call utility.create_sample_output_group {
        input:
          trimmed_bam = ivar_trim.trimmed_bam,
          trimmed_bai = ivar_trim.trimmed_bai,
          untrimmed_bam = bwa.sorted_bam,
          untrimmed_bai = bwa.sorted_bai,
          reference_fasta = first_current_reference_fasta,
          reference_name = first_current_reference_name,
          depth_windows = get_depth_histogram.depth_windows,
          depth_histograms = get_depth_histogram.depth_histogram,
          alignment_stats = get_alignment_stats.alignment_stats,
          amplicon_depths = get_depths.amplicon_depths,
          gene_depths = get_depths.gene_depths,
          sample_name = current_sample
      }
    }
  }
  call utility.create_reference_groups {
    input:
      reference_info_files = flatten(create_sample_output_group.sample_output_group)
  }

  scatter (i in range(length(create_reference_groups.unique_reference_names))) {
    String current_reference_name = create_reference_groups.unique_reference_names[i]
    File current_reference_fasta = create_reference_groups.unique_reference_fastas[i]
    Array[File] current_trimmed_bams = read_lines(create_reference_groups.all_trimmed_bams[i])
    Array[File] current_trimmed_bais = read_lines(create_reference_groups.all_trimmed_bais[i])
    Array[File] current_untrimmed_bams = read_lines(create_reference_groups.all_untrimmed_bams[i])
    Array[File] current_untrimmed_bais = read_lines(create_reference_groups.all_untrimmed_bais[i])
    Array[String] current_sample_names = read_lines(create_reference_groups.all_sample_names[i])

    call bcftools_task.bcftools as bcftools_trimmed {
      input:
        input_bam = current_trimmed_bams, 
        input_bai = current_trimmed_bais,
        reference_fasta = current_reference_fasta,
        output_prefix = current_reference_name + "_trimmed"
    }
    call bcftools_task.bcftools as bcftools_untrimmed {
      input:
        input_bam = current_untrimmed_bams,
        input_bai = current_untrimmed_bais,
        reference_fasta = current_reference_fasta,
        output_prefix = current_reference_name + "_untrimmed"
    } 
    call stats.get_diversity {
      input:
        bcftools_vcf = bcftools_trimmed.filtered_vcf,
        amplicon_bed = reference_files_gcuri + "/" + current_reference_name + "_amplicon.bed",
        reference_gff = reference_files_gcuri + "/" + current_reference_name + ".gff3",
        reference_name = current_reference_name
    } 
    call utility.concatenate_stats_by_reference {
      input:
        depth_windows = read_lines(create_reference_groups.all_depth_windows[i]),
        depth_histograms = read_lines(create_reference_groups.all_depth_histograms[i]),
        alignment_stats =  read_lines(create_reference_groups.all_alignment_stats[i]),
        amplicon_depths =  read_lines(create_reference_groups.all_amplicon_depths[i]),
        gene_depths = read_lines(create_reference_groups.all_gene_depths[i]),
        reference_name = current_reference_name
    }
    # call stats.get_genotyping_report {
    #   input:
    #     trimmed_vcf = bcftools_trimmed.filtered_vcf,
    #     untrimmed_vcf = bcftools_untrimmed.filtered_vcf,
    #     primer_bed = primer_bed,
    #     reference_name = current_reference_name
    # }
    # call stats.get_diversity_metrics {
    #   input:
    #     bcftools_vcf = bcftools_trimmed.filtered_vcf,
    #     depth_file = concatenate_stats_by_reference.depth
    # }
    call stats.call_rsvab {
      input:
        concatenated_alignment_stats = concatenate_stats_by_reference.concatenated_alignment_stats,
    }
  }
  output {
    String rsveillance_version = version_capture.version
    String rsveillance_analysis_date = version_capture.date
    Array[File] final_calls = call_rsvab.final_calls
    Array[File] final_alignment_stats = call_rsvab.final_alignment_stats
    Array[File] alignment_stats = concatenate_stats_by_reference.concatenated_alignment_stats
    Array[File] assembly_fasta = flatten(ivar_consensus.assembly_fasta)
    Array[File] amplicon_depths = concatenate_stats_by_reference.concatenated_amplicon_depths
  }
}