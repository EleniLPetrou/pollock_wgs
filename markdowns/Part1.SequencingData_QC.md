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

## Merge Fastq files for samples that were run on NovaseqSP lane ("Plate3" samples)

According to Erica from the UW Genomics Core, the samples run on a Novaseq XP lane we run without a lane splitter, so they have some fastq files that need to be concatenated. These samples came fom the following populations: BOG02, NBS17, PWS17, and SHE17.

This is how I concatenated the two R1 and R2 fastq files for each sample ("mergelanes.sh")

``` bash
#!/bin/bash
#SBATCH --job-name=elp_mergelanes
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=12:00:00
## Memory per node
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu

##### ENVIRONMENT SETUP ##########
DATADIR=/gscratch/scrubbed/elpetrou/pollock/fastq/Plate3 #directory with fastq files
OUTDIR=/gscratch/scrubbed/elpetrou/pollock/fastq/Plate3_merged_fastq

# Navigate to directory with fastq files and create a text file that contains of the unique individual id basenames (naming convention = POP_INDIV_stuff_LaneID_ReadID_001.fastq)
# cut -d specifies that the delimiter (in our case, and underscore)
# -f 1,2 will retain the first two fields (POP_INDIV) 

cd $DATADIR
ls *R1* | cut -d _ -f 1,2 | sort | uniq > sample_basenames.txt

# sanity check
#cat sample_basenames

# Run loop to merge fastq files
# For each sample, print content to standard output for both lane 1 and lane 2 fastq files.
# Do this separately for the forward  reads (*R1) and the reverse reads (*R2).
# Save the standard output to a new file (thus merging the fastq from Lane 1 and Lane 2 for each sample and read orientation)


for sample in `cat sample_basenames.txt`
do
	cat ${sample}*R1_001.fastq > ${sample}_merged_R1_001.fastq
	cat ${sample}*R2_001.fastq > ${sample}_merged_R2_001.fastq
done

# Move merged files to new directory
mv *merged*.fastq $OUTDIR

```


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
DATADIR=/gscratch/scrubbed/elpetrou/pollock/fastq/Plate3_merged_fastq
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
## Plot distribution of raw sequences per sample

I plotted the distribution of raw sequencing reads for each sample using the R script below and the data summaries produced by multiqc in these files: multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt. 

R script for plotting:

``` R
# The purpose of this script is to plot the distribution of raw sequencing reads for each herring sample. The input data are data frames created by multiqc program.

# Load libraries
library(tidyverse)

# Specify path to input files
MYPATH <- "./sample_metadata/multiqc"
setwd(MYPATH)
list.files()

# Specify file names
file1 <- "Lane1_mqc_fastqc_sequence_counts_plot_1.txt" 
file2 <- "Lane2_mqc_fastqc_sequence_counts_plot_1.txt"
file3 <- "Plate3_mqc_fastqc_sequence_counts_plot_1.txt"

# Read in the data
Lane1_df <- read.delim(file1)
Lane2_df <- read.delim(file2)
Lane3_df <- read.delim(file3)

Lane1_df <- Lane1_df %>%
  separate(Sample, c("pop", "id"), sep = "_", remove = FALSE) %>%
  unite("sample_id", pop:id, sep = "_") %>%
  group_by(sample_id) %>%
  summarise(sum(Unique.Reads)) %>%
  mutate(batch = "Lane1")
  
Lane2_df <- Lane2_df %>%
  separate(Sample, c("pop", "id"), sep = "_", remove = FALSE) %>%
  unite("sample_id", pop:id, sep = "_") %>%
  group_by(sample_id) %>%
  summarise(sum(Unique.Reads)) %>%
  mutate(batch = "Lane2")

Lane3_df <- Lane3_df %>%
  separate(Sample, c("pop", "id"), sep = "_", remove = FALSE) %>%
  unite("sample_id", pop:id, sep = "_") %>%
  group_by(sample_id) %>%
  summarise(sum(Unique.Reads)) %>%
  mutate(batch = "Lane3")

plot_df <- rbind(Lane1_df, Lane2_df, Lane3_df)

plot_df <- plot_df %>%
  mutate(Unique.Reads.Millions = `sum(Unique.Reads)`/1000000) %>%
  filter(sample_id != "Undetermined_S0")
  

# plot the data

distro_plot <- ggplot(data = plot_df) +
  geom_histogram(aes(x = Unique.Reads.Millions), bins = 20) +
  facet_wrap(~batch) +
  ylab("Number of samples") +
  xlab("Number of raw sequences (millions)") +
  theme_bw()
  
distro_plot 

# Calculate some summary statistics
(mymean <- mean(plot_df$Unique.Reads.Millions)) # mean = 12.02 million reads
(mymedian <- median(plot_df$Unique.Reads.Millions)) # mean = 10.45 million reads
(myrange <- range(plot_df$Unique.Reads.Millions)) #range = 0.02 to 49.89 million reads
(mysd <- sd(plot_df$Unique.Reads.Millions)) # sd = 6.22 million reads

# Flad some outlier samples
(lower_threshold <- mymean - 2*mysd) 
(upper_threshold <- mymean + 1*mysd)
  
outlier_df <- plot_df %>%
  filter(Unique.Reads.Millions < lower_threshold | Unique.Reads.Millions > upper_threshold )


# save output of analyses
ggsave("plot_raw_seq_distro.pdf", distro_plot)

```

### Results 

We sequenced 633 pollock and obtained 14.30 million raw sequences per sample on average (median = 10.73 million, range = 0.04 – 61.28 million reads, sd = 9.91 million reads) 

Check out the distribution of raw sequences per sample: 

![raw_seq_plot](https://github.com/EleniLPetrou/pollock_wgs/blob/main/markdowns/plots/raw_seq_distro.png)

Check out the sequencing quality for one of the Lanes (Lane1):

![phred_plot](https://github.com/EleniLPetrou/pollock_wgs/blob/main/markdowns/plots/Lane1_seq_quality.png)



