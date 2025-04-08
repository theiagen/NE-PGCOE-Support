version 1.0

task ivar_variants {
  input {
    File mpileup
    String samplename
    String reference_fasta

    Int min_quality_threshold= 2
    Int min_frequency_threshold = 0.2
    Int min_read_depth = 20

    Int cpu = 2
    Int disk-size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ivar:1.4.2"
    Int memory = 8
  }
  command <<<
    set -euo pipefail

    ivar version | head -n1 | tee IVAR_VERSION

    cat ~{mpileup} | \
      ivar variants \
      -q ~{min_quality_threshold} \
      -t ~{min_frequency_threshold} \
      -m ~{min_read_depth} \
      -r ~{reference_fasta} \
      -p ~{samplename}.ivariants
  >>>
  output {
    String ivar_version = read_string("IVAR_VERSION")
    File ivar_variants = "~{samplename}.ivariants.tsv"
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