version 1.0

import "../tasks/task_versioning.wdl" as versioning

## rsveillance

workflow ysph_rsv {
  input {
    Array[File] read1
    Array[File] read2
    Array[String] samplename

    # questions for seth: 
    #  does he still want it to run on the set or is individual okay for terra?
    #  does he want commented out tasks implemented?
    #  do you want to save the indexes for the references instead of regenerating them each time?

    ### references: RSVA.fasta & RSVB.fasta & assoc. gffs
    ### amplicons/primers: RSVA & RSVB
    ### metadata: r script to pull from gisaid and google sheets -- individual files exist in repo

  }
  call versioning.version_capture {
    input:
  }

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

### bcftools variant calling process

## bcftools call
## input: trimmed bam & bai, reference
## output: bamlist, bcf, vcf
## bcftools mpileup -Ov --threads 3 -o output.bcf -f ref -b bamlist
## bcftools call --threads 3 --ploidy 1 -A -vcO z -o output.vcf output.bcf
## tabix -p vcf output.vcf

## bcftools filter
## input: vcf
## output: filtered vcf & index
## bcftools filter -O z -o filtered.vcf -i 'QUAL>10 & DP>10' input.vcf
## tabix -p vcf filtered.vcf

## repeat previous two steps on UNTRIMMED bamfile -- same outputs but with untrimmed suffix

### stats process

## samtools depth
## input trimmed bam
## output depth file
## samtools depth -a -H input.bam -o depth

## make depth historgram
## input depth file & subsampled stats
## output depthwins.txt & depthhist.txt
## python scripts/get_depth_distribution.py -d depthfile -s samplename -F `cat {input.subfactor} | cut -f 2,2` \ -w {params.winsize} -o {params.prefix} 

## alignment stats
## input subsampled, flagstats, and depthhist.txt
## output alignment stats
## built-in python script

## get amplicon depths
## input bams, substats, depth
## output ampdepth.txt
## scripts/get_depth_windows (see also line 141 of stats.smk)

## get depth of genes
## input: bams, substat, depth & ref gff3
## output: genedepth.txt
## scripts/get_depth_windows but with a gff this time (see also line 165 of stats.smk)

## div_amplicons:
## input vcf & amplicon.bed
## output ampdiv.txt
## scripts/get_refdist_distribution (see line 189)

## div_genes
## input vcf & gff
## output genediv.txt
## scripts/get_refdist_distribution (line 207)

## cat all stats files
## input all stats files calculated (dhists, dwinds, alstats, ampdepths, genedepths)
## output concatenated versions
## cat all files




  output {
    String ysph_rsv_analysis_date = version_capture.date
    String ysph_rsv_wf_version = version_capture.version


# alignstats.txt ->
# consensus.fasta -> from get_ivar_report -> input: ivar (ivariants.tsv & consensus.fa)
# genotype-summary.txt -> from get_genotyping_report -> input: bcftools (all.vcf.gz & untrim_all.vcf.gz) & primer.bed
# final_calls.txt ->
# final_alignstats.tzt -> from call_rsvab -> input: alignstats.txt
# ampdepths -> 

  }
}