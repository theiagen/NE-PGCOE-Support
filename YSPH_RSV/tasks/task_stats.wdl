version 1.0

### stats process

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

task stats {
  input {
    File depth
    File subsampled_stats ??
    File amplicon_bed
    File reference_gff
    File vcf
    File bam
    String samplename

  }
  command <<<
    python3 scripts/get_depth_distribution.py -d ~{depth} -s ~{samplename} -F `cat ~{subsampled_stats} | cut -f 2,2` -w 1000 -o ~{samplename}
    python3 scripts/get_depth_windows.py -d ~{depth} -s ~{samplename} -F `cat ~{subsampled_stats} | cut -f 2,2` -w 1000 -o ~{samplename}
    python3 scripts/get_depth_windows.py -d ~{depth} -s ~{samplename} -F `cat ~{subsampled_stats} | cut -f 2,2` -w 1000 -o ~{samplename} --gff ~{reference_gff}
    python3 scripts/get_refdist_distribution.py -v ~{vcf} -b ~{amplicon_bed} -o ~{samplename}
    python3 scripts/get_refdist_distribution.py -v ~{vcf} -g ~{reference_gff} -o ~{samplename}
    # this one is the built-in python script
    python3 scripts/get_align_stats.py -d ~{depth} -s ~{samplename} -F `cat ~{subsampled_stats} | cut -f 2,2` -w 1000 -o ~{samplename}
  >>>
}