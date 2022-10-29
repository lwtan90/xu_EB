#!/bin/bash

## ENSURE THAT demuxlet can be found in $PATH

## Variable
BAM=possorted_genome_bam.bam
VCF=/labs/joewu/wlwtan/joewu/xu/merged/final.vcf.gz
OUTPUT=demux.newgenome

demuxlet --sam $BAM --vcf $VCF --field GT --out $OUTPUT --group-list barcode.txt




