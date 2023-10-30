#!/bin/sh
query=$1
output_prefix=$2
nucmer -p ${output_prefix} /gpfs1/wucc/Streptococcus_suis-2022-07-22/refs/Streptococcus_suis_BM407.fasta ${query}
delta-filter -1 -q -r ${output_prefix}.delta > ${output_prefix}.filter
show-snps -Clr  ${output_prefix}.filter >  ${output_prefix}.snps
