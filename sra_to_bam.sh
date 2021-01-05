#!/bin/bash

## Run fastq-dump to download multiple sra files

module load perl
### accession list

### accession list
ncbi_dir=$1
output_dir=$2



mkdir -p $output_dir

for sra in $(ls ${ncbi_dir}/*.sra)
do
	base=$(basename $sra)
	sra_id=${base%.sra}
	sbatch -p short -c 4 --mem=60G -t 0-6:00:00 -o ${sra_id}.out --wrap="sam-dump -s -= ${ncbi_dir}/${sra_id}.sra | samtools view -bS - > ${output_dir}/${sra_id}.bam && samtools index ${output_dir}/${sra_id}.bam"
done
