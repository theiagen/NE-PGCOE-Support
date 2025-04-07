version 1.0


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


task bcftools {
  input {
    File input_bam
    File input_bai
    String samplename
    File reference_fasta
  }
  command <<<
    bcftools mpileup -Ov --threads 3 -o ~{samplename}.mpileup -f ~{reference_fasta} ~{input_bam}
    bcftools call --threads 3 --ploidy 1 -A -vcO z -o ~{samplename}_trimmed.vcf ~{samplename}._mpileup
    bcftools filter -O z -o ~{samplename}_trimmed_filtered.vcf -i 'QUAL>10 & DP>10' ~{samplename}_trimmed.vcf
    tabix -p vcf ~{samplename}_trimmed_filtered.vcf
  >>>
  output {
    File filtered_vcf = "~{samplename}_trimmed.vcf"
    File filtered_vcf_index = "~{samplename}_trimmed.vcf.tbi"
  }
}