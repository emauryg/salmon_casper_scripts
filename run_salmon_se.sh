#!/bin/bash

### Run salmon for single end RNA-seq data
## Written by Eduardo Maury (eduardo_maury@hms.harvard.edu)

fastq_dir=$1
outdir=$2
REF_DIR="/n/data1/bch/genetics/lee/eam63/reference/salmon_index/salmon_index/" # this path is hard coded, but it can be updated for future runs. 
mkdir -p $outdir


for fname in $(ls ${fastq_dir}/*.fastq.gz | uniq)
do
	temp=$(basename ${fname})
	fastq_pfx=${temp%.fastq.gz}	

	sbatch -p short -c 4 -t 6:00:00 --mem=30G --wrap="salmon quant -i ${REF_DIR} \
	 -l A \
	 -r ${fastq_dir}/${fastq_pfx}.fastq.gz \
	 -p 4 \
	 -o ${outdir}/${fastq_pfx}_quant \
	 --seqBias \
	 --gcBias \	
	 --useVBOpt \
	 --numBootstraps 30 \
	 --validateMappings"

done