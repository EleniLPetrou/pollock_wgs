#!/bin/bash
#SBATCH --job-name=angsd_window_fst
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=6:00:00
## Memory per node
#SBATCH --mem=300G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu

##### ENVIRONMENT SETUP ##########
MYCONDA=/gscratch/merlab/software/miniconda3/etc/profile.d/conda.sh # path to conda installation on our Klone node. Do NOT change this.
MYENV=angsd_env_0.921 #name of the conda environment containing samtools software. 

## Specify directories with data
DATADIR=/gscratch/scrubbed/elpetrou/pollock/angsd_sfs # Path to directory containing the input data to realSFS program
POP_FILE=pairwise_population_comparisons.txt # text file containing all possible (non-redundant) pairwise population comparisons, in the format pop1 tab pop2

## Specify variables for sliding window analysis in angsd
WINDOW=50000
WINDOW_STEP=10000

##################################################################
## Activate the conda environment:
## start with clean slate
module purge

## This is the filepath to our conda installation on Klone. Source command will allow us to execute commands from a file in the current shell
source $MYCONDA

## activate the conda environment
conda activate $MYENV

####################################################################
## Calculate the 2-D SFS for each pair of populations
cd $DATADIR

while IFS=$'\t' read POP1 POP2 REST; 
do echo $POP1 $POP2; 
realSFS fst stats2 $POP1'.'$POP2'.fst.idx' -win $WINDOW -step $WINDOW_STEP -P ${SLURM_JOB_CPUS_PER_NODE} > $POP1'.'$POP2'.slidingwindow.fst'	
done < $POP_FILE

# Leave conda environment
conda deactivate

#Notes from Merot: Used realSFS functions to compute FST in sliding windows of 25 KB with a step of 5 KB.
#Notes from Clucas: Created Manhattan plots of pairwise FST in non-overlapping 15 KB windows(I am interpreting this to mean that the step size was also 15Kb). Then they annotated SNPs that were in the upper 99.9th percentile of the FST distribution.  
#Notes from Pinsky: They calculated FST using sliding windows of 50 KB and a step of 10KB

# Notes about sliding window output: first column is the specific info on the genomic region, second column is chromosome, third column is the center of the window, fourth column is the number of sites (SNPs), and last column is the average of fst for that window.