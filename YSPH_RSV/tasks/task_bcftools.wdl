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
    bcftools --version | head -n1 | tee VERSION

    echo "DEBUG: creating mpileup file"
    bcftools mpileup \
      -Ov \
      --threads ~{threads} \
      -o ~{output_prefix}_variants.vcf \
      -f ~{reference_fasta} \
      ~{sep=" " input_bam}

    echo "DEBUG: calling variants"
    bcftools call \
      --threads ~{threads} \
      --ploidy ~{ploidy} \
      -A -vcO z \
      -o ~{output_prefix}_all_unfiltered.vcf.gz \
      ~{output_prefix}_variants.vcf

    echo "DEBUG: filtering variants"
    bcftools filter \
      -O z \
      -o ~{output_prefix}_all.vcf.gz \
      -i 'QUAL>10 && DP>10' \
      ~{output_prefix}_all_unfiltered.vcf.gz

    # this line is done because of a malformed info field that has a 'Version="3"' column in the INFO field that causes the scikit-allel to fail
    #   /opt/conda/lib/python3.12/site-packages/allel/io/vcf_read.py:1732: UserWarning: invalid INFO header: '##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">\n'
    #   warnings.warn('invalid INFO header: %r' % header)
    # because of this, I'm making a duplicate version of the VCF that has all INFO fields removed just in case as they are not needed for the analysis
    bcftools view ~{output_prefix}_all.vcf.gz | sed '/^##INFO/d' | bcftools view -O z -o ~{output_prefix}_all_noinfo.vcf.gz
  >>>
  output {
    String bcftools_version = read_string("VERSION")
    File filtered_vcf = "~{output_prefix}_all.vcf.gz"
    File filtered_vcf_noinfo = "~{output_prefix}_all_noinfo.vcf.gz"
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

