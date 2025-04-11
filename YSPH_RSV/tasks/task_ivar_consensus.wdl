version 1.0

task ivar_consensus {
  input {
    File trimmed_bam
    File trimmed_bai
    String output_prefix
    File reference_fasta

    Int mpileup_min_base_quality = 0
    Int mpileup_max_depth = 10000

    Float min_frequency_threshold = 0.75
    Int min_consensus_depth = 20

    Int cpu = 2
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ivar:1.4.2"
    Int memory = 8
  }
  command <<<
    set -euo pipefail

    ivar version | head -n1 | tee IVAR_VERSION
    samtools --version | head -n1 | tee SAMTOOLS_VERSION

    echo "DEBUG: indexing reference genome"
    bwa index ~{reference_fasta}

    echo "DEBUG: generating mpileup"
    samtools mpileup -aa -A \
      -Q ~{mpileup_min_base_quality} \
      -d ~{mpileup_max_depth} \
      -f ~{reference_fasta} \
      -o ~{output_prefix}.mpileup \
      ~{trimmed_bam}

    echo "DEBUG: generating consensus sequence"
    cat ~{output_prefix}.mpileup | ivar consensus \
      -t ~{min_frequency_threshold} \
      -m ~{min_consensus_depth} \
      -p ~{output_prefix}.consensus \
      -i ~{output_prefix}

  >>>
  output {
    String ivar_version = read_string("IVAR_VERSION")
    String samtools_version = read_string("SAMTOOLS_VERSION")
    File mpileup = "~{output_prefix}.mpileup"
    File assembly_fasta = "~{output_prefix}.consensus.fa"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}