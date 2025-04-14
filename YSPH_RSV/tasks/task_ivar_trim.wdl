version 1.0

task ivar_trim {
  input {
    File aligned_bam
    File aligned_bai
    String output_prefix
    File primer_bed

    Int cpu = 2
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ivar:1.4.2"
    Int memory = 8
  }
  command <<<
    set -euo pipefail

    ivar version | head -n1 | tee IVAR_VERSION
    samtools --version | head -n1 | tee SAMTOOLS_VERSION

    echo "DEBUG: trimming primers"
    ivar trim -i ~{aligned_bam} -b ~{primer_bed} -p ~{output_prefix}_trimmed -e 

    echo "DEBUG: sorting and indexing trimmed bam"
    samtools sort -@ 4 -o ~{output_prefix}_trimmed_sorted.bam ~{output_prefix}_trimmed.bam
    samtools index ~{output_prefix}_trimmed_sorted.bam
  >>>
  output {
    String ivar_version = read_string("IVAR_VERSION")
    String samtools_version = read_string("SAMTOOLS_VERSION")
    File trimmed_bam = "~{output_prefix}_trimmed_sorted.bam"
    File trimmed_bai = "~{output_prefix}_trimmed_sorted.bam.bai"
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