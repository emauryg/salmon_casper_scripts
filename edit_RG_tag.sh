#!/bin/bash


## Script to convert bam files to Add RG information

file_list=$1
bam_dir=$2
out_dir=$3

mkdir -p $out_dir

for sra_id in $(cat ${file_list})
do
	sbatch -p short --mem=30G -t 0-6:00:00 -o ${sra_id}.out --wrap="samtools view -H ${bam_dir}/${sra_id}.bam | sed 's,^@RG.*,@RG\tID:None\tSM:None\tLB:None\tPL:Illumina,g' |  samtools reheader - ${bam_dir}/${sra_id}.bam > ${bam_dir}/${sra_id}_RG.bam && samtools index ${bam_dir}/${sra_id}_RG.bam"
done
