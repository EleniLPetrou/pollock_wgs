# The purpose of this script is to calculate a covariance matrix for PCA.
# It uses a beagle file created from angsd to do this.
# For some odd reason I could only get pcangsd to work on an interactive node.

# Request an interactive node on Klone
srun -p compute-hugemem -A merlab --nodes=1 --ntasks-per-node=16 --time=12:00:00 --mem=500G --pty /bin/bash

#####################################################################
##### ENVIRONMENT SETUP ##########

## software information
MYCONDA=/gscratch/merlab/software/miniconda3/etc/profile.d/conda.sh # path to conda installation on our Klone node. Do NOT change this.
MYENV=pcangsd_env #name of the conda environment containing samtools software. 
PCANGSD=/gscratch/merlab/software/miniconda3/envs/pcangsd/pcangsd.py #path to pcangsd software on Klone

## data information
DATADIR=/gscratch/scrubbed/elpetrou/pollock/angsd #path to input files (beagle file)
OUTDIR=/gscratch/scrubbed/elpetrou/pollock/pcangsd #path to output files
MYINFILE=samples617_maf0.01_miss0.3.beagle.gz #name of beagle file
OUTPREFIX=results_pca_samples617_maf0.01_miss0.3 #prefix for output files

##################################################################
## Activate the conda environment:
## start with clean slate
module purge

## This is the filepath to our conda installation on Klone. Source command will allow us to execute commands from a file in the current shell
source $MYCONDA

## activate the conda environment
conda activate $MYENV

## run pcangsd
python $PCANGSD -beagle $DATADIR'/'$MYINFILE -o $OUTDIR'/'$OUTPREFIX -threads 16