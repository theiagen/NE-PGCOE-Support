version 1.0

task mash_index {
  input {
    String reference_prefixes
    String reference_files_gcuri

    String genome_size = "11k"

    Int cpu = 4
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/mash-gcloud:2.3"
    Int memory = 4
  }
  command <<<
    set -euo pipefail

    mash --version | tee MASH_VERSION

    echo "DEBUG: Creating mash indices for the indicated reference prefixes (~{reference_prefixes}) in ~{reference_files_gcuri}"
    FASTAS=(~{reference_prefixes})

    for index in "${!FASTAS[@]}"; do
      echo "DEBUG: Localizing $i reference file: ${FASTAS[$index]}"
      gcloud storage cp ~{reference_files_gcuri}/${FASTAS[$index]}.fasta .

      echo "DEBUG: Creating mash index for ${FASTAS[$index]}"
      mash sketch -g ~{genome_size} ${FASTAS[$index]}.fasta

      echo "DEBUG: updating reference array to point to the sketch"
      FASTAS[$index]=${FASTAS[$index]/.fasta/.fasta.msh}
    done

    echo "DEBUG: merging mash indices"
    mash paste -o index_all.msh ${FASTAS[*]}
  >>>
  output {
    String mash_version = read_string("MASH_VERSION")
    File mash_all_index = "index_all.msh"
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
