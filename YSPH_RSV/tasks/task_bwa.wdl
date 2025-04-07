version 1.0

### bwa alignment process:

## bwa index 
## input: reference fasta
## output: .bwt of fasta
## uses indexer.sh but uses the bwa index for it instead of mash index

## bwa align
## input: r1 & r2 and .bwt of ref fasta
## output: unsorted bam file
## bwa mem -t cpus ref r1 r2 -o output.bam

## samtools flagstat
## input unsorted bam file
## output flagstats
## samtools flagstat -@ cores. -O tsv input.bam > output.flagstat

## samtools remove unaligned from bam
## input: unsorted bam file
## otuput: only aligned bam
## samtools view -b -F 4 -F 2048 -o output.bam input.bam

## subsample bam
## input: only aligned bam
## output: subsampled bam and substat
## subfraction determined by dividing the fastq size by the max fq size of 1024 and returning the floor
## subprocess samtools view -b -s subfraction -o subsampled.bam input.bam

## samtools sort
## input unsorted only aligned bam
## output sorted bam & bai
## samtools sort -@ cores -o sorted.bam input.bam
## samtools index sorted.bam


## samtools depth
## input trimmed bam
## output depth file
## samtools depth -a -H input.bam -o depth


task bwa {
  input {
    File read1
    File read2
    String samplename
    File reference_fasta
    File reference_bwa_index
    
    Int cpu = 4

  }
  command <<<
    bwa mem -t ~{cpu} ~{reference_fasta} ~{read1} ~{read2} -o ~{samplename}.bam

    # generate flagstat
    samtools flagstat -O tsv ~{samplename}.bam > ~{samplename}.flagstat.txt

    # remove unaligned reads
    samtools view -b -F 4 -F 2048 -o ~{samplename}_aligned.bam ~{samplename}.bam

    # subsample bam??
    
    # sorting bam
    samtools sort -@ ~{cpu} -o ~{samplename}_sorted.bam ~{samplename}_aligned.bam
    samtools index ~{samplename}_sorted.bam

    # calculate depth
    samtools depth -a -H ~{samplename}_sorted.bam > ~{samplename}_depth.txt
  >>>
  output {
    File flagstat = "~{samplename}.flagstat.txt"
    File sorted_bam = "~{samplename}_sorted.bam"
    File sorted_bai = "~{samplename}_sorted.bam.bai"
    File depth = "~{samplename}_depth.txt"
  }
  runtime {
    # tbd
  }
}