# Estimating population-specific allele frequencies and site frequency spectra

## Install stable version of angsd (0.921) that can handle -sites flag

My goal is to calculate population-specific allele frequencies and site frequency spectra for each population separately, so I can estimate pairwise FST values and look for signatures of natural selection. I want to do this using the 607K SNPs that I discovered across all individuals. Unfortunately, the most recent version of angsd has an error (related to compatibility with htslib), such that the -sites argument yields 0 sites at the end of the analysis. After much searching on the issues page of angsd on GitHub, I determined that the most recent version of angsd in which the -sites flag works is version 0.921 (with htslib v1.9). Here is how I installed this software on klone:

``` bash
# Try installing a different version that can handle the -sites filter for the sfs analyses
cd /gscratch/merlab/software/miniconda3/envs
conda create -n angsd_env_0.921
conda activate angsd_env_0.921
conda install -c bioconda angsd=0.921 htslib=1.9

```

### Create a sites file for angsd that specifies the SNPs you want to use in the downstream angsd analysis

``` bash
srun -p compute-hugemem -A merlab --nodes=1 --ntasks-per-node=1 --time=01:00:00 --mem=80G --pty /bin/bash

# activate conda angsd

conda activate angsd_env

# Specify paths and file names
DATADIR=/gscratch/scrubbed/elpetrou/pollock/angsd
MAFS_FILE=samples617_maf0.05_miss0.3.nuclear.mafs
SITES_FILE=$MAFS_FILE'.sites'

# Make the sites file
cd $DATADIR
gunzip $MAFS_FILE'.gz'

cut -f 1,2,3,4 $MAFS_FILE > $SITES_FILE 

# index the sites file
angsd sites index $SITES_FILE

```

## Generate site frequency likelihoods for each population and estimate the site frequency spectrum (1-dimensional) for each population

``` bash
#!/bin/bash
#SBATCH --job-name=pollock_angsd_saf
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=6-6:00:00
## Memory per node
#SBATCH --mem=350G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu


## The purpose of the script is to estimate the 1D SFS from angsd

##### ENVIRONMENT SETUP ##########
MYCONDA=/gscratch/merlab/software/miniconda3/etc/profile.d/conda.sh # path to conda installation on our Klone node. Do NOT change this.
MYENV=angsd_env_0.921 #name of the conda environment containing samtools software. 
REALSFS=/gscratch/merlab/software/miniconda3/envs/angsd_env_0.921/bin/realSFS #full path to realSFS program in angsd

## Specify directories with data
REFGENOME=/gscratch/merlab/genomes/atlantic_cod/GCF_902167405.1_gadMor3.0_genomic.fna #path to fasta genome
BAMDIR=/gscratch/scrubbed/elpetrou/pollock/realigned_bam #path to directory containing bam files
SITES_FILE=/gscratch/scrubbed/elpetrou/pollock/angsd/samples617_maf0.05_miss0.3.nuclear.mafs.sites #path to sites file for angsd
OUTDIR=/gscratch/scrubbed/elpetrou/pollock/angsd_sfs #directory for output files

## Specify filtering values for angsd
FILTERS="-minMapQ 30 -minQ 20 -minInd 10 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1"

## Specify output options for angsd
OPTIONS="-GL 1 -doMaf 1 -doSaf 1 -doMajorMinor 3"

## Specify the base name of each of the populations

LIST_POP=(
JPN20
CHU19
NBS17
ZHE17
PRI17
ADK07
ATK07
BOG02
BOG18
MOR17
SAN17
SHU17
SHE17
RES17
PWS17
PWS20
) 

##################################################################
## Activate the conda environment:
## start with clean slate
module purge

## This is the filepath to our conda installation on Klone. Source command will allow us to execute commands from a file in the current shell
source $MYCONDA

## activate the conda environment
conda activate $MYENV

####################################################################
## Calculate the Site Frequency Spectrum for each population separately

## Part1: Create a text file containing the names of the bam files for each population
# try to make an SFS file for one population

## Make a list of bams for each population
cd $BAMDIR

for POP in ${LIST_POP[@]}
do
	ls $POP*'.bam' > $POP'_bams.txt'
done

#Generate site frequency likelihoods using ANGSD
for POP in ${LIST_POP[@]}
do
	angsd -b $BAMDIR/$POP'_bams.txt' -ref $REFGENOME -anc $REFGENOME -out $OUTDIR/$POP \
	-sites $SITES_FILE  \
	$FILTERS \
	$OPTIONS \
	-nThreads ${SLURM_JOB_CPUS_PER_NODE}
done

# Estimate the site frequency spectrum for each population without having to call genotypes or variable sites directly from the site frequency likelihoods

cd $OUTDIR

for POP in ${LIST_POP[@]}
do
	$REALSFS $POP'.saf.idx' > $POP'.sfs'
done

# Deactivate conda environment
conda deactivate
	

```
## Download the results of these analyses to local computer
``` bash
scp elpetrou@klone.hyak.uw.edu:/gscratch/scrubbed/elpetrou/pollock/angsd_sfs/*.global.fst /mnt/hgfs/D
```


## Concantenate results of pairwise population FST and save them to a text file


``` r
# The purpose of this script is to compile data on pairwise population fst 
# that was output by angsd

################################################################################
# Load libraries
library(tidyverse)

# To run this code, put all of your angsd output files in a single directory
DATADIR <- "/gscratch/scrubbed/elpetrou/pollock/angsd_sfs/"

# set working directory
setwd(DATADIR)
#list.files()

# Specify the names of data files used
fileNames <- Sys.glob("*.global.fst") #this is R's version of a wildcard


################################################################################
# Part 1: Create a concatenated dataframe and save it as a text file
# read in the files and start data processing

temp_df <- map(fileNames, read.table, sep = '', header = FALSE) %>%
  set_names(fileNames) %>%
  bind_rows(.id = 'comparison')

output_df <- temp_df %>%
  separate(comparison, c("Pop1", "Pop2"), remove = FALSE) %>%
  rename(unweighted_fst = V1, weighted_fst = V2 )

# save the dataframe as a text file
write.table(output_df, file = "pairwise_population_FST_concatenated_results.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE)

```
