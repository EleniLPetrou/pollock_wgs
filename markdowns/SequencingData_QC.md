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
## Plot distribution of raw sequences per sample

I plotted the distribution of raw sequencing reads for each sample using the R script below and the data summaries produced by multiqc in these files: multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt. 

Here is a summary of the raw sequencing data :

We sequenced 633 pollock and obtained 12.02 million raw sequences per sample on average (median = 10.45 million, range = 0.02 – 49.89 million reads, sd = 6.22 million reads) 

R script for plotting:

``` R
# The purpose of this script is to plot the distribution of raw sequencing reads for each herring sample. The input data are data frames specifying the raw number of reads per sample (saved as .csv files). These data were provided to me by the NW Genomics Sequencing Center.

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
Download the full [MultiQC html reports here](https://github.com/EleniLPetrou/herring_whole_genome_sequencing/blob/main/Markdown/multiqc_report.html)

The sequencing quality looks really great for almost all samples, hurray!

![phred_plot](https://github.com/EleniLPetrou/herring_whole_genome_sequencing/blob/main/Markdown/plots/plot_fastqc_mean_qual_scores_raw.png)



