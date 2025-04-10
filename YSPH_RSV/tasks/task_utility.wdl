version 1.0

task create_reference_groups {
  input {
    Array[File] reference_info_files

    Int cpu = 1
    Int memory = 1
    Int disk_size = 10
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0"
  }
  command <<<
    for file in ~{sep=" " reference_info_files}; do
      reference_name=$(cut -f1 $file)
      reference_fasta=$(cut -f2 $file)
      trimmed_bam=$(cut -f3 $file)
      trimmed_bai=$(cut -f4 $file)
      untrimmed_bam=$(cut -f5 $file)
      untrimmed_bai=$(cut -f6 $file)
      depth_windows=$(cut -f7 $file)
      depth_histograms=$(cut -f8 $file)
      alignment_stats=$(cut -f9 $file)
      amplicon_depths=$(cut -f10 $file)
      gene_depths=$(cut -f11 $file)
      sample_name=$(cut -f12 $file)

      echo $reference_name >> all_reference_names.txt
      echo $reference_fasta >> all_reference_fastas.txt
      echo $trimmed_bam >> ${reference_name}_trimmed_bam.txt
      echo $trimmed_bai >> ${reference_name}_trimmed_bai.txt
      echo $untrimmed_bam >> ${reference_name}_untrimmed_bam.txt
      echo $untrimmed_bai >> ${reference_name}_untrimmed_bai.txt
      echo $depth_windows >> ${reference_name}_depth_windows.txt
      echo $depth_histograms >> ${reference_name}_depth_histograms.txt
      echo $alignment_stats >> ${reference_name}_alignment_stats.txt
      echo $amplicon_depths >> ${reference_name}_amplicon_depths.txt
      echo $gene_depths >> ${reference_name}_gene_depths.txt
      echo $sample_name >> ${reference_name}_sample_name.txt

    done

    sort -u all_reference_names.txt > unique_reference_names.txt
    sort -u all_reference_fastas.txt > unique_reference_fastas.txt

    ls *_trimmed_bam.txt > all_trimmed_bams.txt
    ls *_trimmed_bai.txt > all_trimmed_bais.txt
    ls *_untrimmed_bam.txt > all_untrimmed_bams.txt
    ls *_untrimmed_bai.txt > all_untrimmed_bais.txt
    ls *_depth_windows.txt > all_depth_windows.txt
    ls *_depth_histograms.txt > all_depth_histograms.txt
    ls *_alignment_stats.txt > all_alignment_stats.txt
    ls *_amplicon_depths.txt > all_amplicon_depths.txt
    ls *_gene_depths.txt > all_gene_depths.txt
    ls *_sample_name.txt > all_sample_names.txt

  >>>
  output {
    Array[String] unique_reference_names = read_lines("unique_reference_names.txt")
    Array[File] unique_reference_fastas = read_lines("unique_reference_fastas.txt")
    Array[File] all_trimmed_bams = read_lines("all_trimmed_bams.txt")
    Array[File] all_trimmed_bais = read_lines("all_trimmed_bais.txt")
    Array[File] all_untrimmed_bams = read_lines("all_untrimmed_bams.txt")
    Array[File] all_untrimmed_bais = read_lines("all_untrimmed_bais.txt")
    Array[File] all_depth_windows = read_lines("all_depth_windows.txt")
    Array[File] all_depth_histograms = read_lines("all_depth_histograms.txt")
    Array[File] all_alignment_stats = read_lines("all_alignment_stats.txt")
    Array[File] all_amplicon_depths = read_lines("all_amplicon_depths.txt")
    Array[File] all_gene_depths = read_lines("all_gene_depths.txt")
    Array[File] all_sample_names = read_lines("all_sample_names.txt")
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}

task create_sample_output_group {
 input {
    String reference_name
    String reference_fasta
    String trimmed_bam
    String trimmed_bai
    String untrimmed_bam
    String untrimmed_bai
    String depth_windows
    String depth_histograms
    String alignment_stats
    String amplicon_depths
    String gene_depths
    String sample_name
    
    Int cpu = 1
    Int memory = 1
    Int disk_size = 10
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0"
  }
  command <<<
    echo -e "~{reference_name}\t~{reference_fasta}\t~{trimmed_bam}\t~{trimmed_bai}\t~{untrimmed_bam}\t~{untrimmed_bai}\t{~depth_windows}\t~{depth_histograms}\t~{alignment_stats}\t~{amplicon_depths}\t~{gene_depths}\t~{sample_name}" > reference_info.txt
  >>>
  output {
    File sample_output_group = "sample_reference_info.txt"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}


task concatenate_stats_by_reference {
  input {
    Array[File] depth_windows
    Array[File] depth_histograms
    Array[File] alignment_stats
    Array[File] amplicon_depths
    Array[File] gene_depths

    String reference_name
    
    Int cpu = 1
    Int memory = 2
    Int disk_size = 25
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/rsveillance:0.1"
  }
  command <<<
    echo "DEBUG: concatenating depth windows"
    cat ~{sep=" " depth_windows} > ~{reference_name}_depth_windows.txt

    echo "DEBUG: concatenating depth histograms"
    cat ~{sep=" " depth_histograms} > ~{reference_name}_depth_histograms.txt

    echo "DEBUG: concatenating alignment stats"
    cat ~{sep=" " alignment_stats} > ~{reference_name}_alignment_stats.txt

    echo "DEBUG: concatenating amplicon depths"
    cat ~{sep=" " amplicon_depths} > ~{reference_name}_amplicon_depths.txt

    echo "DEBUG: concatenating gene depths"
    cat ~{sep=" " gene_depths} > g~{reference_name}_gene_depths.txt

  >>>
  output {
    File concatenated_depth_windows = "~{reference_name}_depth_windows.txt"
    File concatenated_depth_histograms = "~{reference_name}_depth_histograms.txt"
    File concatenated_alignment_stats = "~{reference_name}_alignment_stats.txt"
    File concatenated_amplicon_depths = "~{reference_name}_amplicon_depths.txt"
    File concatenated_gene_depths = "~{reference_name}_gene_depths.txt"
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