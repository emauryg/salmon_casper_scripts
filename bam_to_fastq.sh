#!/bin/bash

## Script to convert bam files to fastq files

bam_dir=$1
out_dir=$2

mkdir -p $out_dir

for fname in $(ls ${bam_dir}/*.bam)
do
	temp=$(basename ${fname})
	sra_id=${temp%.bam}
	sbatch -p short --mem=30G -t 0-6:00:00 -o ${sra_id}.out --wrap="java -jar $PICARD SamToFastq I=${bam_dir}/${sra_id}.bam FASTQ=${out_dir}/${sra_id}_1.fastq.gz F2=${out_dir}/${sra_id}_2.fastq.gz"
done


