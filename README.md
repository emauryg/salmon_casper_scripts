# Scripts to run Salmon and Casper on bulk RNA-seq data

The workflow assumes either single-end or paired-end sequencing data. Scripts should work on most HPC slurm systems, but it is worthwhile to check commands for compatability.
We employed several hard-coded paths in these scripts that should be adapted to each user's needs. 

A typical workflow will look something like this:

1. Run `run_salmon_se.sh` or `run_salmon.sh` for single end and paired end respectively
2. 
