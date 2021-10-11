# Estimate read depth bam files after indel realignment

Using samtools, I estimated the read depth of each bam file after indel realignment using the script below ("samtools_depth.sh"). As far as I could tell, samtools depth command does not support multithreading, so it took about ~2 min to run for each sample. 

``` bash
#!/bin/bash
#SBATCH --job-name=pollock_depth
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=3-12:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu


##### ENVIRONMENT SETUP ##########
## Specify the directory containing data
DATADIR=/mmfs1/gscratch/scrubbed/elpetrou/pollock/realigned_bam #directory with sam files
SUFFIX1=_realigned.bam #file suffix
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
    samtools depth -aa $MYSAMPLEFILE | cut -f 3 | gzip > $MYSAMPLEFILE'.depth.gz'
done

## Flag explanations for samtools depth:
## -aa: output absolutely all positions, including unused ref. sequences

## deactivate the conda environment
conda deactivate

```

# Parse the output and calculate summary statistics

Next, I used the Rscript "estimate_realigned_sequencing_depth.R" to take the sequencing depth output of samtools and summarize it. 

``` R
# The purpose of this script is to compile and plot data on sequencing read depth
# after bam files have been filtered and indel realignment has taken place. 
# As input is takes the .depth files that are output from samtools depth.

################################################################################
# Load libraries
library(tidyverse)

################################################################################
# read in the files and start data processing

command_args <- commandArgs(trailingOnly = TRUE)
fileName <- command_args[1]

# For each file, read in the depth data and calculate summary statistics
# Compute sequencing depth summary statistics

depth <- read_tsv(fileName, col_names = F)$X1
mean_depth <- mean(depth)
sd_depth <- sd(depth)
mean_depth_nonzero <- mean(depth[depth > 0])
mean_depth_within2sd <- mean(depth[depth < mean_depth + 2 * sd_depth])
median <- median(depth)
presence <- as.logical(depth)
proportion_of_reference_covered <- mean(presence)
  
# save these results to a dataframe
output <- data.frame(fileName, mean_depth, sd_depth, 
                       mean_depth_nonzero, mean_depth_within2sd, 
                       median, proportion_of_reference_covered)
  
output_df <- output %>%
    mutate(across(where(is.numeric), round, 3)) %>%
    separate(fileName, "population", sep = "_", remove = FALSE,  extra = "drop")
  
# write the output_df to a text file
write.table(output_df, file = paste0(fileName,"depth_results.txt"), 
              append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = TRUE)
  
#remove the giant intermediate files that take up a lot of memory
remove(depth)
remove(presence)
  
#clean up the workspace memory
gc()
  
```

I ran the above R script on klone, using the following sbatch script ("sbatch_estimate_realigned_sequencing_depth.sh"):

``` bash
#!/bin/bash
#SBATCH --job-name=pollock_R_seq_depth
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=4-12:00:00
## Memory per node
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu


##### ENVIRONMENT SETUP ##########
## Specify the directory containing data
DATADIR=/mmfs1/gscratch/scrubbed/elpetrou/pollock/realigned_bam #directory with depth files created by samtools
MYSCRIPT=/mmfs1/home/elpetrou/pollock_scripts/plot_realigned_sequencing_depth.R #path to R script
SUFFIX1=.bam.depth.gz #suffix of files that you would like to analyze

## start with clean slate
module purge

###################################################################################################################
## Move into the working directory and run script
cd $DATADIR


## Run R script
for MYSAMPLEFILE in *$SUFFIX1
do
    Rscript --vanilla $MYSCRIPT $MYSAMPLEFILE
done

```

Following this, I plotted the output data using another R script ("plot_realigned_sequencing_depth.R").

``` R
# The purpose of this script is to compile and plot data on sequencing read depth
# after bam files have been filtered and indel realignment has taken place. 
# As input is takes the .depth files that are output from samtools depth.

################################################################################
# Load libraries
library(tidyverse)

# To run this code, put all of your depth files in a single directory

DATADIR <- "/gscratch/scrubbed/elpetrou/pollock/realigned_bam/"

# set working directory
setwd(DATADIR)
#list.files()

# Specify the names of data files used
fileNames <- Sys.glob("*.gzdepth_results.txt") #this is R's version of a wildcard


################################################################################
# Part 1: Create a concatenated dataframe and save it as a text file
# read in the files and start data processing

output_df = data.frame() #initialize empty dataframe


for (fileName in fileNames) {
  print(fileName) #counter
  df <- read.delim(fileName)
  output_df <- rbind(output_df, df) # add each individual dataframe to a big dataframe
}


# Mean depth and standard deviation over all individuals
output_df$mean_depth <- as.numeric(output_df$mean_depth)
output_df$sd_depth <- as.numeric(output_df$sd_depth)

# save the dataframe as a text file
write.table(output_df, file = "sequencing_depth_after_realignment_concatenated_results.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE)


################################################################################
# Part2 : Plot the depth distribution

# read in the concatenated dataframe

#output_df <- read.delim("sequencing_depth_after_realignment_concatenated_results.txt")

#allow scientific notation
options(scipen = 0) 

plot1 <- ggplot(output_df, aes(x = population, color = population)) +
  geom_point(aes(x = population, y = mean_depth)) +
  geom_boxplot(aes(x = population, y = mean_depth)) +
  ylab("sequencing depth") +
  xlab("population") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

plot1


# save the plot as a pdf
ggsave("sequencing_depth_after_realignment.pdf", plot1)

```
# Results

The mean sequencing depth across the entire reference genome was 0.8X (median = 0.79, range = 0.24 – 4.73X).

Information to add: 1) plot 2) how low-depth individuals were identified and removed from downstream analyses 3) list of removed infividuals
