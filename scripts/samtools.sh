#!/bin/bash
#SBATCH --job-name=pollock_samtools
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=6-12:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu


##### ENVIRONMENT SETUP ##########
## Specify the directory containing data
DATADIR=/gscratch/scrubbed/elpetrou/pollock/sam #directory with sam files
SUFFIX1=_RG.sam #file suffix
MYCONDA=/gscratch/merlab/software/miniconda3/etc/profile.d/conda.sh # path to conda installation on our Klone node. Do NOT change this.
MYENV=samtools_env #name of the conda environment containing samtools software. 

## Activate the conda environment:
## start with clean slate
module purge

## This is the filepath to our conda installation on Klone. Source command will allow us to execute commands from a file in the current shell
source $MYCONDA

## activate the conda environment
conda activate $MYENV


###################################################################################################################
## Move into the working directory and run script
cd $DATADIR

## Run samtools commands. This takes about 5 min per sample (so like 2 days total for whole data set?)
for MYSAMPLEFILE in *$SUFFIX1
do
    echo $MYSAMPLEFILE
    MYBASE=`basename --suffix=$SUFFIX1 $MYSAMPLEFILE`
    samtools view -bS -F 4 $MYBASE'_RG.sam' > $MYBASE'.bam'
    samtools view -h -q 20 $MYBASE'.bam' | samtools view -buS - | samtools sort -o $MYBASE'_minq20_sorted.bam'
    samtools index $MYBASE'_minq20_sorted.bam'
done

## Flag explanations for samtools view:
## -b       output BAM
## -h       include header in SAM output
## -q INT   only include reads with mapping quality >= INT [0]
##-F INT   only include reads with none of the bits set in INT set in FLAG [0] (aka when this is set to 4, you remove unmapped reads)

## deactivate the conda environment
conda deactivate

