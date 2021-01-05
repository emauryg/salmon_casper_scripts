#!/bin/bash

## Run salmon on a batch of fastq files from RNA-seq experiments

fastq_dir=$1
outdir=$2
### accession list
#sample_list=$3

mkdir -p $outdir

for fname in $(ls ${fastq_dir}/*_1.fastq.gz | uniq)
do
	temp=$(basename ${fname})
	sra_id=${temp%_1.fastq.gz}	
	sbatch -o ${sra_id}.out run_salmon.sh $sra_id $fastq_dir $outdir
done