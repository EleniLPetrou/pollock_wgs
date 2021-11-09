# The purpose of this script is to compile data on pairwise population fst 
# that was output by angsd

################################################################################
# Load libraries
library(tidyverse)


# To run this code, put all of your angsd output files in a single directory
DATADIR <- "E:/Dropbox (MERLAB)/Eleni/pollock/results/angsd/samples617_miss0.3_maf0.05_nuclear"

# set working directory
setwd(DATADIR)
list.files()

# Specify the names of data files used
fileNames <- Sys.glob("*.global.fst") #this is R's version of a wildcard

fileNames
################################################################################
# Part 1: Create a concatenated dataframe and save it as a text file
# read in the files and start data processing

temp_df <- map(fileNames, read.table, sep = ' ', header = FALSE) %>%
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



### extra code for debugging
output_df = data.frame() #initialize empty dataframe


for (fileName in fileNames) {
  print(fileName) #counter
  df <- read.delim(fileName)
  output_df <- rbind(output_df, df) # add each individual dataframe to a big dataframe
}


