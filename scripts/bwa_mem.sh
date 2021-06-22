#!/bin/bash
#SBATCH --job-name=pollock_bwamem
#SBATCH --account=merlab
#SBATCH --partition=ckpt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=2:00:00
## Memory per node
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/elpetrou/pollock/fastq_trimmed/fastq_PWS20


##### ENVIRONMENT SETUP ##########
MYCONDA=/gscratch/merlab/software/miniconda3/etc/profile.d/conda.sh # path to conda installation on our Klone node. Do NOT change this.
MYENV=bwa_env #name of the conda environment containing samtools software. 

## Specify directories with data
REFGENOME=/gscratch/merlab/genomes/atlantic_cod/GCF_902167405.1_gadMor3.0_genomic.fna #path to ref genome that was indexed by bwa
SUFFIX1=_R1_001.trim.fastq # Suffix to trimmed fastq files. The forward reads with paired-end data.
SUFFIX2=_R2_001.trim.fastq # Suffix to trimmed fastq files. The reverse reads with paired-end data.

##################################################################
## Activate the conda environment:
## start with clean slate
module purge

## This is the filepath to our conda installation on Klone. Source command will allow us to execute commands from a file in the current shell
source $MYCONDA

## activate the conda environment
conda activate $MYENV

###################################

# Save the base name of each input file
MYBASE=$(basename --suffix=$SUFFIX1 "$1")

# Sanity check
echo "$1"
echo $MYBASE

## options explained
#-M : mark shorter split hits as secondary
# -t INT        number of threads [1]
#<index_prefix> is the index for the reference genome generated from bwa index
# input_reads_pair_1.fastq, input_reads_pair_2.fastq are the input files of sequencing data that is paired-end
# -o FILE       sam file to output results to [stdout]

#bwa mem -M -t 1 <idxbase> <in1.fq> [in2.fq]

bwa mem -M -t ${SLURM_JOB_CPUS_PER_NODE} $REFGENOME ${MYBASE}$SUFFIX1 ${MYBASE}$SUFFIX2 -o ${MYBASE}.sam

