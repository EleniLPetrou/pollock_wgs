#!/bin/bash
#SBATCH --job-name=pollock_multiqc
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=3-12:00:00
## Memory per node
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu

##### ENVIRONMENT SETUP ##########
DATADIR1=/gscratch/scrubbed/elpetrou/pollock/fastqc/Plate3
DATADIR2=/gscratch/scrubbed/elpetrou/pollock/fastqc/Lane1
DATADIR3=/gscratch/scrubbed/elpetrou/pollock/fastqc/Lane2

MYCONDA=/gscratch/merlab/software/miniconda3/etc/profile.d/conda.sh # path to conda installation on our Klone node. Do NOT change this.
MYENV=multiqc_env #name of the conda environment containing samtools software.


################################################################################
## Activate the conda environment:

## start with clean slate
module purge

## This is the filepath to our conda installation on Klone. Source command will allow us to execute commands from a file in the current shell
source $MYCONDA

## activate the conda environment
conda activate $MYENV

# Run multiqc
cd $DATADIR1
multiqc .

cd $DATADIR2
multiqc .

cd $DATADIR3
multiqc .
