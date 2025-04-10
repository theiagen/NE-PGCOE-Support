version 1.0

task bwa {
  input {
    File read1
    File read2
    String output_prefix
    File reference_fasta
    
    Int cpu = 6
    Int memory = 16
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ivar:1.4.2"
  }
  command <<<
    set -euo pipefail

    echo "BWA $(bwa 2>&1 | grep Version )" | tee BWA_VERSION
    samtools --version | head -n1 | tee SAMTOOLS_VERSION

    echo "DEBUG: indexing reference genome ~{reference_fasta}"
    bwa index ~{reference_fasta}

    echo "DEBUG: aligning reads to reference genome"
    bwa mem \
      -t ~{cpu} \
      ~{reference_fasta} \
      ~{read1} ~{read2} \
      -o ~{output_prefix}.bam

    echo "DEBUG: generating flagstat"
    samtools flagstat \
      -O tsv ~{output_prefix}.bam \
      > ~{output_prefix}.flagstat.txt

    echo "DEBUG: removing unaligned reads from bam"
    samtools view -b \
      -F 4 \
      -F 2048 \
      -o ~{output_prefix}_aligned.bam \
      ~{output_prefix}.bam

    echo "DEBUG: sorting and indexing aligned bamfile"
    samtools sort \
      -@ ~{cpu} \
      -o ~{output_prefix}_sorted.bam \
      ~{output_prefix}_aligned.bam
    samtools index ~{output_prefix}_sorted.bam

    echo "DEBUG: calculating depth"
    samtools depth -a -H ~{output_prefix}_sorted.bam > ~{output_prefix}_depth.txt
  >>>
  output {
    String bwa_version = read_string("BWA_VERSION")
    String samtools_version = read_string("SAMTOOLS_VERSION")
    File flagstat = "~{output_prefix}.flagstat.txt"
    File sorted_bam = "~{output_prefix}_sorted.bam"
    File sorted_bai = "~{output_prefix}_sorted.bam.bai"
    File depth_file = "~{output_prefix}_depth.txt"
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