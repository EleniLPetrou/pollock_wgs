#!/bin/bash
#SBATCH --job-name=pollock_bwaindex
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=2-12:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu


##### ENVIRONMENT SETUP ##########
MYCONDA=/gscratch/merlab/software/miniconda3/etc/profile.d/conda.sh # path to conda installation on our Klone node. Do NOT change this.
MYENV=bwa_env #name of the conda environment containing samtools software. 

## Specify directories with data
REFGENOME=/gscratch/merlab/genomes/atlantic_cod/GCF_902167405.1_gadMor3.0_genomic.fna #path to fasta genome
MYPREFIX=GCF_902167405.1_gadMor3.0_genomic
##################################################################
## Activate the conda environment:
## start with clean slate
module purge

## This is the filepath to our conda installation on Klone. Source command will allow us to execute commands from a file in the current shell
source $MYCONDA

## activate the conda environment
conda activate $MYENV

## Index the genome using bwa

bwa index -a bwtsw $REFGENOME $MYPREFIX