#!/bin/bash

## Trying samdump with sratoolkits 2.10.8

#!/bin/bash

## Run prefect to download mutliple sra files using sratoolkit 2.10.8
#module load sratoolkit/2.9.0


### accession list
sample_list=$1 

#rm -r -f /n/scratch2/eam63/ncbi/ncbi_downloads/

for sra_id in $(cat ${sample_list})
do
	sbatch -p short -c 1 --mem=60G -t 0-12:00:00 -o ${sra_id}.out --wrap="sam-dump ${sra_id}  | samtools view -bS - > ${sra_id}.bam && samtools index ${sra_id}.bam"
done
