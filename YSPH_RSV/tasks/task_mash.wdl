version 1.0

  
### mash process:

## index references w/ mash
## input: reference fastas
## output: msh (index_all.msh)
## uses indexer.sh w/ var parameters:
        # genome_size = "11k"
    # alternatively, there's an option to make a bwa index for the reference fasta
    # making mash sketch:  mash sketch -g ${genome_size} ${FASTAS[$i]}
    # merging mash indices: mash paste $mashout ${FASTAS[*]}

## mash input: r1 & r2
## mash output: mash calls and mash output txt
## uses masher.sh script w/ var parameters:
        # reads=10000, # compare top N reads to refs
        # bloom=10,    # bloom filter kmers with < N coverage (seq errors)
        # gsize="11k", # estimated genome size (for prob assignment)
        # prob=1e-50,  # max mash prob to call
        # dist=0.25,   # max mash dist to call
    # mash dist -m bloom -r -g gzise ref reads.fastq > mash.txt
    # extracting max prob / dist:
    # awk -v dist=dist -v prob=prob -v sample=sample '($3+0 < dist+0 && $4+0 < prob+0) {sub(".fasta","",$1); print sample, $1}' ${PREFIX}_mash.txt > ${PREFIX}_calls.txt

## mash calls concatenated together for all samples
## input: mash results from all samples
## output; concatenated mash outputs w/ cat


task mash {
  input {
    File read1
    File read2
    String samplename
    File reference_mash_index

    Int number_of_reads = 10000
    Int bloom_filter = 10
    String genome_size = "11k"
    Float max_mash_prob = 1e-50
    Float max_mash_dist = 0.25
  }
  command <<<
    # get the number_of_reads reads ?? seems odd - will ask

    # concatenated read1 and read2 is input into this task
    mash dist -m ~{bloom_filter} -r -g ~{genome_size} ~{reference_mash_index} {read1} {read2} > ~{samplename}_mash.txt

    awk -v dist=~{max_mash_dist} -v prob=~{max_mash_prob} -v sample=~{samplename} '($3+0 < dist+0 && $4+0 < prob+0) {sub(".fasta","",$1); print sample, $1}' ~{samplename}_mash.txt > ~{samplename}_calls.txt
  >>>
  output {
    File mash_calls = "~{samplename}_calls.txt"
    File mash_output = "~{samplename}_mash.txt"
  }
  runtime {
    # tbd
  }
}