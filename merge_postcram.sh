#!/bin/bash

## Code to merge bam files after being extracted from a cram
## This code was used for the anaplistic samples fom collord et al

bam_dir=$1
out_dir=$2

sample_list=$(ls ${bam_dir}/ | cut -d"_" -f1 | uniq)
subsamp=(1 2)

mkdir -p $out_dir

for i in ${sample_list[@]}
do
	echo $i
	for j in ${subsamp[@]}
		do
			echo $j
			sbatch -o ${i}_pair1.out -p short -t 0-3:00:00 --mem 20G --wrap="sambamba merge ${out_dir}/${i}_1.merged.bam ${bam_dir}/${i}_1#1.bam ${bam_dir}/${i}_2#1.bam"
			sbatch -o ${i}_pair2.out -p short -t 0-3:00:00 --mem 20G --wrap="sambamba merge ${out_dir}/${i}_2.merged.bam ${bam_dir}/${i}_1#2.bam ${bam_dir}/${i}_2#2.bam"
		done	
done

