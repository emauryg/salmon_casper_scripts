#!/bin/bash

## Sort with samtools
## Example: ./sort_bams.sh ../bams/Agnihotri_etal/ ../EGA_downloads/Agnihotri_etal/sample_names.txt

REFERENCE=/n/data1/bch/genetics/lee/eam63/reference/Homo_sapiens_assembly19.fasta
BAM_dir=$1
sample_names=$2

for i in `cat ${sample_names} | uniq`
do
	sbatch -c 1 -p short -t 0-6:00:00 --mem 10G -o ${i}.sort.out \
	 --wrap="sambamba index 6_Aligned.sortedByCoord.out.bam"
done