# Load libraries
library(tidyverse)
library(reshape2)

# Specify data directory containing input file
DATADIR <- "~/pollock"

# Specify name of input file
INFILE <- "./results/pairwise_population_FST_concatenated_results.txt" #tab-delimited text file containing fst data
METAFILE <- "./sample_metadata/collection_site_metadata.txt" # full path to sampling location metadata

OUTFILE <- "pairwise_population_FST_samples617_maf0.05_miss0.3.nuclear.pdf" # Specify name of output file

# Specify a custom order for the populations in the heatmap (WA by spawn time, then AK by spawn time)

my_levels <- c("JPN20", "CHU19", "NBS17", "ZHE17", "PRI17", "ADK07", 
               "ATK07", "BOG02", "BOG18", "MOR17", "SAN17", "SHU17",
               "SHE17", "RES17", "PWS17", "PWS20")

my_levels2 <- c("Japan", "Chukchi", "N. Bering", "Zhemchug", 
                "Pribilof", "Adak", "Atka", 
                "Bogoslof.02","Bogoslof.18", "Morzhovoi", 
                "Sanak", "Shumigan", "Shelikof", "Resurrection",
                "PWS.17", "PWS.20")

##############################################################################
# set working directory
setwd(DATADIR)

# Read in the data and manipulate it for plotting
fst_df <- read.table(INFILE, header = TRUE)
meta_df <- read.delim(METAFILE)

mini_df <- meta_df %>%
  select(pop, location)


# make a temporary df with the population names joined and bind them together, 
# to make the full pairwise matrix.

temp_df <- fst_df %>%
  rename(Pop1 = Pop2, Pop2 = Pop1)

full_df <- rbind(fst_df, temp_df)

full_df$weighted_fst <- round(full_df$weighted_fst, digits = 3)

# Order the levels according to a custom order  

full_df$Pop1 <- factor(x = full_df$Pop1,
                       levels = my_levels, 
                       ordered = TRUE)

full_df$Pop2 <- factor(x = full_df$Pop2,
                       levels = my_levels, 
                       ordered = TRUE)

###Part 2: remove duplicate pairwise-columns

# Turn the dataframe into a matrix

my_mat <- acast(full_df, Pop1~Pop2, value.var = "weighted_fst")

## Specify some functions to retrieve upper part of matrix
# Get lower triangle of the correlation matrix

get_lower_tri <- function(Fstmat){
  Fstmat[upper.tri(Fstmat)] <- NA
  return(Fstmat)
}

## subset the matrix
lower_tri <- get_lower_tri(my_mat)
lower_tri

##Use the package reshape to melt the matrix into a df again:
final_df <- melt(lower_tri, value.name = "weighted_fst")

plotting_df <- left_join(final_df, mini_df, by = c("Var1" = "pop")) %>%
  rename(Pop1 = location) %>%
  left_join(mini_df, by = c("Var2" = "pop")) %>%
  rename(Pop2 = location)


plotting_df$Pop1 <- factor(plotting_df$Pop1, levels = my_levels2)
plotting_df$Pop2 <- factor(plotting_df$Pop2, levels = my_levels2)

# Make a heatmap and visualize the FST values

heatmap_plot <- ggplot(data = plotting_df, aes(Pop1, Pop2, fill = weighted_fst)) +
  geom_raster() +
  geom_text(aes(label = weighted_fst), size = 2) +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(angle = 0, vjust = 1, size = 10, hjust = 1)) +
  ylab("Population A") +
  xlab("Population B") +
  labs(fill = expression(italic(F[ST]))) +
  coord_fixed() 

heatmap_plot

# save pdf to file

ggsave(OUTFILE, heatmap_plot)
