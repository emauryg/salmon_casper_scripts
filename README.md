# Scripts to run Salmon and Casper on bulk RNA-seq data

The workflow assumes either single-end or paired-end sequencing data. Scripts should work on most HPC slurm systems, but it is worthwhile to check commands for compatability.
We employed several hard-coded paths in these scripts that should be adapted to each user's needs. 

A typical workflow will look something like this:

1. Run `run_salmon_se.sh` or `run_salmon.sh` for single end and paired end respectively. Make sure you have the correct index libraries
2. Run `star_align_se.sh` or `star_align.sh` make sure you have the correct index libraries.
3. Run `bafextract.sh`
4. Run `casper_running_script.R`, make sure to change the paths to where your salmon and bafextract outputs are located. 

**Note:** controls samples were generated with a combination of samples to be used across tumors. You might need to build your own control reference TPM values and bafextract outputs following a protocol similar to how we did it for tumors.