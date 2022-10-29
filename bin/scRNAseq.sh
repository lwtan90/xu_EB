#!/bin/bash


#### This pipeline will take FASTQ file ####

#### Modules
module load cellranger

S1=$1
mkdir $S1

#### Global
REF=/labs/joewu/wlwtan/annotation/hg38/10x/hg38
FQ=/labs/joewu/wlwtan/joewu/xu/EB_nov2021/X2/usftp21.novogene.com/raw_data/$S1

S1=$1
mkdir $S1
cellranger count --id=$S1 --include-introns --transcriptome=$REF --fastqs=$FQ --localcores=5 --localmem=100 --sample=$S1

