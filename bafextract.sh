#!/bin/bash

GENOME_LIST=/n/scratch3/users/e/eam63/clinical_projects/meningioma/RNAseq/references/hg19.list
FASTA_pileup=/n/scratch3/users/e/eam63/clinical_projects/meningioma/RNAseq/references/hg19/

bam_dir=$1
output_dir=$2

rm -rf $output_dir
mkdir -p $output_dir


cmd=/n/data1/bch/genetics/lee/eam63/tools/BAFExtract-master/bin/BAFExtract



for bam_file in `find ${bam_dir} -type f -name "*.bam"`
do
	base=$(basename ${bam_file})
	file=${base%_Aligned.sortedByCoord.out.bam}
	sbatch -o ${file}.out -p short -t 0-3:00:00 --mem 30G --wrap="samtools view ${bam_dir}/${file}_Aligned.sortedByCoord.out.bam | ${cmd}  -generate_compressed_pileup_per_SAM stdin ${GENOME_LIST} ${output_dir} 50 0; ${cmd} -get_SNVs_per_pileup  ${GENOME_LIST} ${output_dir} ${FASTA_pileup} 20 4 0.1 ${output_dir}/${file}.bcf"
done

