#!/bin/bash
#SBATCH --job-name=pollock_angsd_2dsfs
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=9-12:00:00
## Memory per node
#SBATCH --mem=300G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu

##### ENVIRONMENT SETUP ##########
MYCONDA=/gscratch/merlab/software/miniconda3/etc/profile.d/conda.sh # path to conda installation on our Klone node. Do NOT change this.
MYENV=angsd_env_0.921 #name of the conda environment containing samtools software. 

## Specify directories with data
DATADIR=/gscratch/scrubbed/elpetrou/pollock/angsd_sfs # Path to directory containing the input data to realSFS program
POP_FILE=pairwise_population_comparisons.txt # text file containing all possible (non-redundant) pairwise population comparisons, in the format pop1 <tab> pop2

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


## read in the tab-delimited text file that has format pop1 <tab> pop2,  line by line
## split each line by a tab, first field = POP1, second field = POP2
## run realSFS module on angsd for POP1 and POP2 to estimate the 2-dimensional frequency spectrum from the site allele frequency likelihoods

while IFS=$'\t' read POP1 POP2 REST; 
do echo $POP1 $POP2; 
realSFS $POP1'.saf.idx' $POP2'.saf.idx' -P ${SLURM_JOB_CPUS_PER_NODE} -maxIter 30 > $POP1'.'$POP2'.ml'	
done < $POP_FILE

## Estimate the pairwise Fst between each pair of populations
## First we will index the sample so that the same sites are analysed for each population
## Then we will get the global estimate of FST between each population pair

while IFS=$'\t' read POP1 POP2 REST; 
do echo $POP1 $POP2; 
realSFS fst index $POP1'.saf.idx' $POP2'.saf.idx' -sfs $POP1'.'$POP2'.ml' -fstout $POP1'.'$POP2
realSFS fst stats $POP1'.'$POP2'.fst.idx' > $POP1'.'$POP2'.global.fst'	
done < $POP_FILE

# Leave conda environment
conda deactivate

