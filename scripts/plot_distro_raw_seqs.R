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
