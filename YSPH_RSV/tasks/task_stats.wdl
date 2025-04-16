version 1.0

task get_depth_histogram {
  input {
    File depth_file
    String reference_name
    String associated_sample

    Int window_size = 10
    Int downsample_fraction = 1

    Int cpu = 1
    Int memory = 4
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/rsveillance:0.1"
  }
  String output_prefix = "~{associated_sample}_~{reference_name}"
  command <<<
    echo "DEBUG: calculating the depth windows and depth histograph"
    get_depth_distribution.py \
      -d ~{depth_file} \
      -s ~{associated_sample} \
      -F ~{downsample_fraction} \
      -w ~{window_size} \
      -o ~{output_prefix}
  >>>
  output {
    String stats_docker = docker
    File depth_windows = " ~{output_prefix}_depthwins.txt"
    File depth_histogram = " ~{output_prefix}_depthhist.txt"
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

task get_alignment_stats {
  input {
    File alignment_flagstat
    File depth_histogram
    String reference_name
    String associated_sample

    Int downsample_fraction = 1
    Int min_depth = 10

    Int cpu = 1
    Int memory = 4
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/rsveillance:0.1"
  }
  String output_prefix = "~{associated_sample}_~{reference_name}"
  command <<<
    echo "DEBUG: calculating alignment stats"
    get_align_stats.py \
      -s ~{associated_sample} \
      -t ~{reference_name} \
      -F ~{downsample_fraction} \
      -i ~{alignment_flagstat} \
      -d ~{depth_histogram} \
      -o ~{output_prefix} \
      -m ~{min_depth}
  >>>
  output {
    File alignment_stats = "~{output_prefix}_alignstats.txt"
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

task get_depths {
  input {
    File depth_file    
    File amplicon_bed
    File reference_gff

    String reference_name
    String associated_sample

    Int downsample_fraction = 1

    Int cpu = 1
    Int memory = 4
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/rsveillance:0.1"
  }
  String output_prefix = "~{associated_sample}_~{reference_name}"
  command <<<
    echo "DEBUG: calculating amplicon depths"
    get_depth_windows.py \
      -b ~{amplicon_bed} \
      -d ~{depth_file} \
      -s ~{associated_sample} \
      -F ~{downsample_fraction} \
      -o ~{output_prefix}_ampdepth.txt

    echo "DEBUG: calculating gene depths"
    get_depth_windows.py \
      -g ~{reference_gff} \
      -d ~{depth_file} \
      -s ~{associated_sample} \
      -F ~{downsample_fraction} \
      -o ~{output_prefix}_genedepth.txt
  >>>
  output {
    File amplicon_depths = "~{output_prefix}_ampdepth.txt"
    File gene_depths = "~{output_prefix}_genedepth.txt"
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

task get_diversity {
  input {
    File bcftools_vcf
    File amplicon_bed
    File reference_gff
    String reference_name

    Int downsample_fraction = 1

    Int cpu = 1
    Int memory = 4
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/rsveillance:0.1"
  }
  command <<<
    echo "DEBUG: calculating amplicon nucleotide diversity (pi)"
    get_refdist_distribution.py \
      -v ~{bcftools_vcf} \
      -b ~{amplicon_bed} \
      -o ~{reference_name}_ampdiv.txt 

    echo "DEBUG: calculating amplicon nucleotide diversity (pi)"
    get_refdist_distribution.py \
      -v ~{bcftools_vcf} \
      -g ~{reference_gff} \
      -o ~{reference_name}_ampdiv.txt 
  >>>
  output {
    File amplicon_diversity = "~{reference_name}_ampdiv.txt"
    File gene_diversity = "~{reference_name}_genediv.txt"
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

task get_diversity_metrics {
  input {
    File bcftools_vcf
    File depth_file
    File primer_bed
    String reference_name

    Int cpu = 1
    Int memory = 4
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/rsveillance:0.1"
  }
  command <<<
    echo "DEBUG: calculating diversity metrics"
    get_pidp_ranges.py \
      --vcf ~{bcftools_vcf} \
      --depth ~{depth_file} \
      --bed ~{primer_bed} \
      --out ~{reference_name}_primer-stats
  >>>
  output {
    File mean_depth = "~{reference_name}_primer-stats_meandp.tsv"
    File mean_depth_plot = "~{reference_name}_primer-stats_meandp.png"
    File covpc = "~{reference_name}_primer-stats_covpc.tsv"
    File covpc_plot = "~{reference_name}_primer-stats_covpc.png"
    File pi = "~{reference_name}_primer-stats_pi.tsv"
    File pi_plot = "~{reference_name}_primer-stats_pi.png"
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

task get_genotyping_report {
  input {
    File trimmed_vcf
    File untrimmed_vcf
    File primer_bed
    String reference_name

    Int cpu = 1
    Int memory = 4
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/rsveillance:0.1"
  }
  command <<<
    echo "DEBUG: calculating genotyping report"


    get_pi_ranges.py \
      --vcf ~{trimmed_vcf} \
      --bed ~{primer_bed} \
      --out ~{reference_name}_primer

    python3 get_pi_ranges.py \
      --vcf ~{untrimmed_vcf} \
      --bed ~{primer_bed} \
      --out ~{reference_name}_untrim-primer
  >>>
  output {
    File trimmed_pi = "~{reference_name}_primer-pi.tsv"
    File untrimmed_pi = "~{reference_name}_untrim-primer-pi.tsv"
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

task call_rsvab {
 input {
    File concatenated_alignment_stats

    Float min_coverage = 0.8
    Float min_ratio = 0.95

    Int cpu = 1
    Int memory = 4
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/rsveillance:0.1"
  }
  command <<<
    echo "DEBUG: calling RSVA/RSVB"
    call_RSVAB.py \
      --alignstats ~{concatenated_alignment_stats} \
      --coverage ~{min_coverage} \
      --ratio ~{min_ratio} \
      --out "final"
  >>>
  output {
    File final_calls = "final_calls.txt"
    File final_alignment_stats = "final_alignstats.txt"
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