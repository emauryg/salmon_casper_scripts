#!/bin/bash

## Run prefect to download mutliple sra files using sratoolkit 2.10.8
#module load sratoolkit/2.9.0


### accession list
sample_list=$1 

#rm -r -f /n/scratch2/eam63/ncbi/ncbi_downloads/

for sra_id in $(cat ${sample_list})
do
	sbatch -p short -c 1 --mem=60G -t 0-12:00:00 -o ${sra_id}.out --wrap="prefetch ${sra_id}"
done


# ## Code to run failed downloads
# for sra_id in $(grep "err" -io -l *.out | cut -d"." -f1 )
# do
# 	sbatch -p short -c 1 --mem=60G -t 0-12:00:00 -o ${sra_id}.out --wrap="prefetch -v ${sra_id}"
# done

# ## Download via https if you get this "error: file not found while copying file - cannot download 'SRR3996005' using requested transport"

# prefetch --option-file ../sra_files/clark_etal/SRR_Acc_List.txt
