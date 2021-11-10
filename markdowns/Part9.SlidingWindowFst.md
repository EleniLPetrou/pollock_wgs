# Sliding window analysis in angsd

Once I had computed the 2D site frequency spectrum (2D SFS) for each pair of populations in angsd (see Step8.md), I wanted to conduct a sliding window analysis to visualize genetic differentiation across the genome. To do this, I used the realSFS program in angsd. The code for this is below. Here are some additional notes from other studies that have done this analysis:

- Notes from Merot et al. 2021: Used realSFS functions to compute FST in sliding windows of 25 KB with a step of 5 KB.
- Notes from Clucas et al. 2019: Created Manhattan plots of pairwise FST in non-overlapping 15 KB windows(I am interpreting this to mean that the step size was also 15Kb). Then they annotated SNPs that were in the upper 99.9th percentile of the FST distribution.  
- Notes from Pinsky et al. 2020: They calculated FST using sliding windows of 50 KB and a step of 10KB

Here is the code I ran on Klone:

``` bash
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

# Notes about sliding window output: first column is the specific info on the genomic region, second column is chromosome, third column is the center of the window, fourth column is the number of sites (SNPs), and last column is the average of fst for that window.
```

# Download the data to local computer
``` 
scp elpetrou@klone.hyak.uw.edu:/gscratch/scrubbed/elpetrou/pollock/angsd_sfs/*slidingwindow.fst /hgfs/mnt/D
```


# Plot the results of the sliding window analysis in R

``` r
# The purpose of this script is to take sliding window output from angsd and plot it. The plots compare one focal population to all others. To run this script you need: 1) angsd output files in a single directory, 2) a table containing information about how you would like the chromosomes to be named, 3) the name of the focal population, and 4) the name of all populations in your analysis.


# Load libraries
library(tidyverse)
library(viridis)
library(scales)

# Specify working directory and input files
DATADIR <- "E:/Dropbox (MERLAB)/Eleni/postdoc_wgs/results/angsd/all_samples_maf0.05_miss_0.3/WA_comparisons"

CHROMFILE <- "E:/Dropbox (MERLAB)/Eleni/postdoc_wgs/chrom_lookup_table.txt" #tab-delimited table that has refseq chrom name and simplified name (LG1, etc.)

# Specify the name of the outpt file (pdf)
OUTFILE <- "plot_fst_sliding_windows_grid_WASamples.pdf"

# Specify the population you would like to focus on, and a full list of population names
FOCALPOP <- "CHPT16"
POPLIST <- as.vector(c("SQUA14", "PORT14", "SMBY15", "QLBY19", "ELBY15", "CHPT16"))

#################################################################################
#################################################################################
# Read in files
setwd(DATADIR)

# Remove the FOCALPOP from the POPLIST
COMPARISONS <- subset(POPLIST, !(POPLIST %in% FOCALPOP))

# Create a search term (WILDCARD) that will help you find all angsd output files that contain the FOCALPOP in the file name.
WILDCARD <- paste0("*",FOCALPOP,"*","slidingwindow.fst")


FILENAMES <- Sys.glob(WILDCARD) #this is R's version of a wildcard
file_list <- as.list(FILENAMES)

chrom_df <- read.table(CHROMFILE, header = TRUE)

################################################################################
# Part 1: Create a concatenated dataframe and save it as a text file
# This part of the code is just parsing the angsd files, concatenating them, and adding some information that will make plotting easier. 

fst_df <- file_list %>%
  set_names(nm = FILENAMES) %>%
  map_dfr(
    ~ read_delim(.x, skip = 1, col_types = cols(), col_names = c("region", "chr", "midpos", "nsites", "fst"), delim = "\t"),
    .id = "comparison"
  )

head(fst_df)


# Append information about simplified chromosome names (from RefSeq to numeric)
fst_df <- left_join(fst_df, chrom_df, by = "chr")
nCHR <- length(unique(fst_df$chr))

# Calculate the distance along chromosomes in Mb
fst_df <- fst_df %>%
  mutate(midpos_Mb = midpos/1000000)

head(fst_df)

# Strip the original file name such that only the population that FOCALPOP is being compared to is retained (for plot labelling)

final_df <- fst_df %>%
  mutate(temp1 = gsub(".slidingwindow.fst", "", comparison)) %>%
  mutate(temp2 = gsub(FOCALPOP, "", temp1)) %>%
  mutate(temp3 = gsub("\\.", "", temp2))

# Make the POP2 column in the dataframe into a factor with a specific order
final_df$POP2 <- factor(final_df$temp3, levels = COMPARISONS)

# Make the linkage_group column in the dataframe into a factor with a specific order
final_df$linkage_group <- factor(final_df$linkage_group, levels = chrom_df$linkage_group)

#################################################################
# Part 2: Plot the data

# Pick a color scheme and labelling scheme
show_col(viridis_pal(option = "plasma")(20))
mypalette <- c("#0D0887FF", "#0D0887FF", "#0D0887FF", "#D14E72FF", "#D14E72FF", "#FA9E3BFF") #specify viridis

color_df <- as.data.frame(cbind(POPLIST, mypalette)) %>%
  filter(POPLIST %in% COMPARISONS)

mycols <- as.vector(color_df$mypalette)

manhplot <- ggplot() +
  geom_point(data = final_df, aes(x = midpos_Mb, y = fst, color = POP2), 
             alpha = 1, size = 0.5) +
  scale_color_manual(values = mycols) +
  facet_grid(POP2 ~ linkage_group, scales = "free_x") +
  ylab(expression(italic(F[ST]))) +
  xlab("Chromosome position (Mb)") +
  theme_bw() +
  ggtitle(FOCALPOP) +
  # set tick mark spacing
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4)) +
  scale_x_continuous(breaks = c(15,30)) +
  theme( 
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, size = 7, vjust = 0.5),
    panel.background = element_rect(fill = "white"), 
    panel.spacing = unit(0,"lines")
  )


# Save the plot to a pdf file

ggsave(paste0(FOCALPOP,"_", OUTFILE), plot = manhplot, width = 11, height = 5, units = "in")

```
