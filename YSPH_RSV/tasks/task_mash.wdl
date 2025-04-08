version 1.0

task mash {
  input {
    File read1
    File read2
    String samplename

    File reference_mash_index
    String reference_prefixes
    String reference_files_gcuri

    Int number_of_reads = 10000
    Int bloom_filter = 10
    String genome_size = "11k"
    Float max_mash_prob = 1e-50
    Float max_mash_dist = 0.25

    Int cpu = 2
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/mash-gcloud:2.3"
    Int memory = 4
  }
  command <<<
    set -euo pipefail

    mash --version | tee MASH_VERSION

    echo "DEBUG: Extracting top ~{number_of_reads} reads from input files"
    let "LINES = ~{number_of_reads} / 2 * 4"
    for file in ~{read1} ~{read2}; do
      BASENAME=$(basename ${file/fastq.gz/head.fastq})
      gunzip -dc $file | head -n $LINES >> ~{samplename}_~{number_of_reads}.fastq
    done

    echo "DEBUG: mashing reads against the reference hashes"
    mash dist \
      -m ~{bloom_filter} \
      -r -g ~{genome_size} \
      ~{reference_mash_index} \
      ~{samplename}_~{number_of_reads}.fastq \
      > ~{samplename}_mash.txt

    echo "DEBUG: Pulling matches below ~{max_mash_dist}/~{max_mash_prob} from mash output"
    awk -v dist=~{max_mash_dist} -v prob=~{max_mash_prob} -v sample=~{samplename} '($3+0 < dist+0 && $4+0 < prob+0) {sub(".fasta","",$1); print sample, $1, $3}' ~{samplename}_mash.txt | sort -k3,3n > ~{samplename}_calls.txt

    echo "DEBUG: Identifying the most likely reference"
    FILENAME=$(head -n1 ~{samplename}_calls.txt | cut -f2)
    REF_BASENAME=$(basename ${FILENAME/.fasta/})

    echo "DEBUG: extracting path of the reference files for downstream usage"
    REFERENCE_FASTA="~{reference_files_gcuri}/${REF_BASENAME}.fasta"
    REFERENCE_GFF="~{reference_files_gcuri}/${REF_BASENAME}.gff"
  >>>
  output {
    String mash_version = read_string("MASH_VERSION")
    File mash_calls = "~{samplename}_calls.txt"
    File mash_output = "~{samplename}_mash.txt"
    File reference_fasta = read_string("REFERENCE_FASTA")
    File reference_gff = read_string("REFERENCE_GFF")
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