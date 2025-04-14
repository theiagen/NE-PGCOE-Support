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
    String max_mash_prob = "1e-50"
    Float max_mash_dist = 0.25

    Int cpu = 2
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/mash-gcloud:2.3"
    Int memory = 4
  }
  command <<<
    mash --version | tee MASH_VERSION

    echo "DEBUG: Extracting top ~{number_of_reads} reads from input files"
    let "LINES = ~{number_of_reads} / 2 * 4"
    for file in ~{read1} ~{read2}; do
      gunzip -dc $file | head -n $LINES >> ~{samplename}_~{number_of_reads}.fastq
    done

    set -euo pipefail

    echo "DEBUG: mashing reads against the reference hashes"
    mash dist \
      -m ~{bloom_filter} \
      -r -g ~{genome_size} \
      ~{reference_mash_index} \
      ~{samplename}_~{number_of_reads}.fastq \
      > ~{samplename}_mash.txt

    echo "DEBUG: Pulling matches below ~{max_mash_dist}/~{max_mash_prob} from mash output"
    awk -v dist=~{max_mash_dist} -v prob=~{max_mash_prob} -v sample=~{samplename} '($3+0 < dist+0 && $4+0 < prob+0) {sub(".fasta","",$1); print sample, $1}' ~{samplename}_mash.txt > ~{samplename}_calls.txt

    while read -r line; do
      FILENAME=$(echo "$line" | cut -d' ' -f2)

      echo "DEBUG: extracting path of the reference files for downstream usage"
      REFERENCE_FASTA="~{reference_files_gcuri}/${FILENAME}.fasta"
      REFERENCE_GFF="~{reference_files_gcuri}/${FILENAME}.gff"
      echo "$REFERENCE_FASTA" >> ~{samplename}_fastas.txt
      echo "$REFERENCE_GFF" >> ~{samplename}_gffs.txt
      echo "${FILENAME}" >> ~{samplename}_references.txt
      echo "~{samplename}" >> ~{samplename}_samples.txt
    done < ~{samplename}_calls.txt
  >>>
  output {
    String mash_version = read_string("MASH_VERSION")
    File mash_calls = "~{samplename}_calls.txt"
    File mash_output = "~{samplename}_mash.txt"
    Array[File] reference_fastas = read_lines("~{samplename}_fastas.txt")
    Array[File] reference_gffs = read_lines("~{samplename}_gffs.txt")
    Array[String] reference_names = read_lines("~{samplename}_references.txt")
    Array[String] associated_samples = read_lines("~{samplename}_samples.txt")
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