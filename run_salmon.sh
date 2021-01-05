#!/bin/bash

## Run Salmon to quantify transcripts on Paired end RNA-seq data. 

#SBATCH -c 4
#SBATCH -t 0-06:00:00
#SBATCH --mem=30G
#SBATCH -p short

fastq_pfx=$1
fastq_dir=$2
outdir=$3

REF_DIR="/n/data1/bch/genetics/lee/eam63/projects/lncRNA_project/RNA_seq/reference/" # this path is hard coded, but it can be updated for future runs. 

salmon quant -i ${REF_DIR} \
 -l A \
 -1 ${fastq_dir}/${fastq_pfx}_1.fastq.gz \
 -2 ${fastq_dir}/${fastq_pfx}_2.fastq.gz \
 -p 4 \
 -o ${outdir}/${fastq_pfx}_quant \
 --seqBias \
 --gcBias \
 --useVBOpt \
 --numBootstraps 30 \
 --validateMappings

