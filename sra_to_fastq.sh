#!/bin/bash

## Run fastq-dump to download multiple sra files

module load perl
### accession list

### accession list
sample_list=$1 
output_dir=$2
ncbi_dir=$3



mkdir -p $output_dir

for sra_id in $(cat ${sample_list})
do
	sbatch -p short -c 4 --mem=60G -t 0-6:00:00 -o ${sra_id}.out --wrap="fasterq-dump ${ncbi_dir}/${sra_id}.sra -t /tmp -e 4 -O ${output_dir} --split-3 --skip-technical -M 1 && gzip ${output_dir}/${sra_id}*fastq*"
done


