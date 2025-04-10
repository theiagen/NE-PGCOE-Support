version 1.0

task depth_histogram {
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
  command <<<
    echo "DEBUG: calculating the depth windows and depth histograph"
    python3 scripts/get_depth_distribution.py \
      -d ~{depth_file} \
      -s ~{associated_sample} \
      -F ~{downsample_fraction} \
      -w ~{window_size} \
      -o ~{associated_sample}_~{reference_name}
  >>>
  output {
    File depth_windows = "~{samplename}_depthwins.txt"
    File depth_histogram = "~{samplename}_depthhist.txt"
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