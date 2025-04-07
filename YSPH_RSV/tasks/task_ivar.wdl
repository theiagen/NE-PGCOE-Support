version 1.0


### ivar consensus process:

## primer clip 
## input sorted bam & primers
## output: unsorted trimmed bam, sorted trimmed bam & bai
## ivar trim -i input.bam -b input.primers -o unsorted.bam
## samtools sort unsorted bam & samtools index on sorted bam

## samtools mpileup
## input: trimmed bam & bai
## output: mpileup
## samtools mpileup -aa -A -Q 0 -d maxdepth(10000) -f ref -o mpileup input.bam

## ivar consensus
## input mpileup
## output consensus.fa
## cat mpileup | ivar consensus -t threshold (0.75) -m depth (20) -p output.prefix -i samplename

## ivar variants
## input mpileup
## output variants.tsv
## cat mpileup | ivar variants -q qual (2) -r ref -t threshold (0.2) -m depth (20) -p prefix


task ivar_trim {
  input {
    File aligned_bam
    File primer_bed
    String samplename
    File reference_fasta

  }
  command <<<
    ivar trim -i ~{aligned_bam} -b ~{primer_bed} -o ~{samplename}_trimmed.bam -e 
    samtools sort -@ 4 -o ~{samplename}_trimmed_sorted.bam ~{samplename}_trimmed.bam
    samtools index ~{samplename}_trimmed_sorted.bam
    samtools mpileup -aa -A -Q 0 -d 10000 -f ~{reference_fasta} -o ~{samplename}_mpileup.txt ~{samplename}_trimmed_sorted.bam


  >>>
  output {
    File trimmed_bam = "~{samplename}_trimmed_sorted.bam"
    File trimmed_bai = "~{samplename}_trimmed_sorted.bam.bai"
    File mpileup = "~{samplename}_mpileup.txt"
  }
  runtime {
    # tbd
  }
}

task ivar_consensus {
  input {
    File mpileup
    String samplename
  }
  command <<<
    cat ~{mpileup} | ivar consensus -t 0.75 -m 20 -p ~{samplename} -i ~{samplename}
  >>>
  output {
    File consensus = "~{samplename}.fa"
  }
  runtime {
    # tbd
  }
}

task ivar_variants {
  input {
    File mpileup
    File reference_fasta
    String samplename
  }
  command <<<
    cat ~{mpileup} | ivar variants -q 2 -r ~{reference_fasta} -t 0.2 -m 20 -p ~{samplename}
  >>>
  output {
    File variants = "~{samplename}.tsv"
  }
  runtime {
    # tbd
  }
}