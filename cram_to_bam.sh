#!/bin/bash

## Script to combert cram files to bam files after downloading from EGA

cram_dir=$1
output_dir=$2
fnames=$(find ${cram_dir} -mindepth 1 -maxdepth 1 -type d  -printf '%f\n')

rm -r $output_dir
mkdir -p $output_dir

for f in ${fnames}
do
	cram_name=$(ls ${cram_dir}/${f}/*.cram)
	base_name=$(basename $cram_name)
	pfx=${base_name%.cram}
	sbatch -o ${pfx}.out -p short -t 0-3:00:00 --mem 30G --wrap="samtools view -b  -o ${output_dir}/${pfx}.bam $cram_name && \
	samtools index ${output_dir}/${pfx}.bam"
done
