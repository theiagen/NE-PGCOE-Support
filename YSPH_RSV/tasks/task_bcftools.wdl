version 1.0

task bcftools {
  input {
    Array[File] input_bam
    Array[File] input_bai
    File reference_fasta
    String output_prefix

    Int ploidy = 1

    Int cpu = 4
    Int threads = 3
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/bcftools:1.21"
    Int memory = 24
  }
  command <<<
    echo ~{sep="\n" input_bam} > bam_list.txt

    echo "DEBUG: creating mpileup file"
    bcftools mpileup \
      -Ov \
      --threads ~{threads} \
      -o ~{output_prefix}_variants.vcf \
      -f ~{reference_fasta} \
      -b bam_list.txt

    echo "DEBUG: calling variants"
    bcftools call \
      --threads ~{threads} \
      --ploidy ~{ploidy} \
      -A -vcO z \
      -o ~{output_prefix}_all_unfiltered.vcf.gz \
      ~{output_prefix}_variants.vcf

    echo "DEBUG: indexing variants"
    tabix -p vcf ~{output_prefix}_all_unfiltered.vcf.gz

    echo "DEBUG: filtering variants"
    bcftools filter \
      -O z \
      -o ~{output_prefix}_all.vcf.gz \
      -i 'QUAL>10 && DP>10' \
      ~{output_prefix}_all_unfiltered.vcf.gz

    echo "DEBUG: indexing filtered variants"
    tabix -p vcf ~{output_prefix}_all.vcf.gz
  >>>
  output {
    File filtered_vcf = "~{output_prefix}_all.vcf.gz"
    File filtered_vcf_index = "~{output_prefix}_all.vcf.gz.tbi"
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

