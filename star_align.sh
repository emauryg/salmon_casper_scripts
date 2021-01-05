#!/bin/bash

## Running STAR alignment using STAR 2.7

REFERENCE=/n/data1/bch/genetics/lee/eam63/reference/Homo_sapiens_assembly19.fasta
GENOMEDIR=/n/data1/bch/genetics/lee/eam63/reference/STAR_genomedir/
GTF=/n/data1/bch/genetics/lee/eam63/reference/gencode.v33lift37.annotation_fixed.gtf
FASTQ_dic=$1
output_dir=$2
	
rm -r -f $output_dir
mkdir -p $output_dir

for file in `ls ${FASTQ_dic}*_1.fastq* | uniq`
do
	base_name=$(basename ${file})
	i=${base_name%_1.fastq.gz}
	sbatch -c 4 -p short -t 0-6:00:00 --mem 70G -o ${i}.star.out --wrap="STAR --twopassMode Basic \
	--twopass1readsN -1 --runThreadN 4 \
	--runMode alignReads --genomeDir ${GENOMEDIR} \
	--sjdbGTFfile ${GTF} --sjdbOverhang 75 \
	--readFilesIn ${FASTQ_dic}/${i}_1.fastq.gz ${FASTQ_dic}${i}_2.fastq.gz \
	--outFileNamePrefix ${output_dir}/${i}_  \
	--outSAMtype BAM SortedByCoordinate \
	--limitBAMsortRAM 45000000000 \
	--limitOutSJcollapsed 1000000 \
	--outTmpKeep None \
	--readFilesCommand gunzip -c"
done


