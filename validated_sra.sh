#!/bin/bash

## Check that all .sra files were downloaded correctly
module load sratoolkit/2.9.0


### accession list
sample_list=$1 


for sra_id in $(cat ${sample_list})
do
	sbatch -p short -c 1 --mem=60G -t 0-1:00:00 -o ${sra_id}.out --wrap="vdb-validate ${sra_id}"
done

