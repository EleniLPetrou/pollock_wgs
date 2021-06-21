# WGS data

## Data download and backup

Sequencing data were downloaded from the NW Genomics Center using Globus software. The raw data (.tar files) have been saved to three locations:
  1. An external hard drive that is in Ellie's posession
  2. The HauserLab fireproof external hard drive: ioSafe_HauserLab/
  3. Klone: /gscratch/scrubbed/elpetrou/pollock
  4. LOLO Archive: /archive/merlab/pollock_wgs/raw_sequences


For each of these data backups, I used the MD5sum file to verify that the data were not corrupted. 


## Unzip data
I unzipped the data on Klone using the folowing scipt:


``` bash
#!/bin/bash
#SBATCH --job-name=pollock_gunzip
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

# Move the folders into a new directory
PLATE3=/gscratch/scrubbed/elpetrou/pollock/fastq/Plate3
LANE1=/gscratch/scrubbed/elpetrou/pollock/fastq/Lane1
LANE2=/gscratch/scrubbed/elpetrou/pollock/fastq/Lane2

# unzip the fastq files in each directory
cd $PLATE3
gunzip *fastq.gz

cd $LANE1
gunzip *fastq.gz

cd $LANE2
gunzip *fastq.gz

```

## Plot distribution of raw sequences per sample

I plotted the distribution of raw sequencing reads for each sample using the R script below and data provided to me by the NW Genomics Sequencing Center in .csv format (these were included in the tar files containing raw seqencing data). 

Here is a summary of the raw sequencing data :
  - We sequenced 556 herring (number of WA samples = 281; number of AK samples = 275)
  - Average number of reads per sample = 12.88 million reads
  - Range = 0 to 80.06 million reads
  - Standard deviation = 4.41 million reads
  - Seven samples sequenced poorly (reads per sample less than 2 sd from the mean)
  -     SMBY15_012  SMBY15_018  SQUA14_045   BERN16_004  BERN16_010  BERN16_032 BERN16_031
  - Nine samples sequenced very deeply (reads per sample more than 2 sd from mean)
  -     CHPT16_041  ELBY15_118  ELBY15_129  ELBY15_177  QLBY19_080  QLBY19_097  BERN16_033  OLGA19_003  SITKA17_043


## Run FastQC

To check the quality of the raw sequence data I ran the software FastQC on Klone

``` bash
#!/bin/bash
#SBATCH --job-name=pollock_fastqc
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
DATADIR=/gscratch/scrubbed/elpetrou/pollock/fastq/Plate3
OUTDIR=/mmfs1/gscratch/scrubbed/elpetrou/pollock/fastqc/Plate3
MYSUFFIX=.fastq

#### CODE FOR JOB #####
## make output directory
mkdir $OUTDIR

## navigate into directory holding fastq files
cd $DATADIR 

# Run fastqc on all files that end in a particular suffix ( in this case, .fastq)
for SAMPLEFILE in *$MYSUFFIX
do
    fastqc -f fastq --extract \
    $SAMPLEFILE
done

## Move all the fastqc results files (ending in _fastqc, .zip, .html) to the output directory
mv *.html $OUTDIR
mv *.zip $OUTDIR
mv *_fastqc $OUTDIR

```

## Visualize FastQC output using MultiQC

I visualized the voluminous FastQC output using MultiQC software on Klone

First, I installed multiqc version 1.10.1 using conda:

``` bash

#Create a conda environment for multiqc
conda create -n multiqc_env

# To activate or enter the environments you just created simply type:
conda activate multiqc_env

# Once in the environment, install the versions of software that you want using conda:
conda install -c bioconda -c conda-forge multiqc

```

Next, I ran multiqc in the folder that contained the fastqc output:

``` bash
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

```


### Results 
Download the full [MultiQC html report here](https://github.com/EleniLPetrou/herring_whole_genome_sequencing/blob/main/Markdown/multiqc_report.html)

The sequencing quality looks really great for almost all samples, hurray!

![phred_plot](https://github.com/EleniLPetrou/herring_whole_genome_sequencing/blob/main/Markdown/plots/plot_fastqc_mean_qual_scores_raw.png)



