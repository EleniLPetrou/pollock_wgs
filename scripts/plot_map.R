# This R  script was written by Eleni Petrou on 20191003 using  "R version 3.6.1 (2019-07-05)"
# Its purpose is to make a map of sampling locations using GPS data

###########################################################
# Load R libraries

# Libraries
library(tidyverse)
library(ggrepel)
library(viridis)
library(lubridate)
library(grid)
library(cowplot)
library(maps)

# Specify directory names and file names
BASEDIR <- "E:/Dropbox (MERLAB)/Eleni/pollock"

#tab-delimited text file with sampling metadata
INFILE <- "./sample_metadata/collection_site_metadata.txt" 

# name of pdf to output
OUTFILE <- "./results/plot_map.pdf"

setwd(BASEDIR)
###########################################################

# Read in your data frame with longitude, latitude, and other metadata.
coord_df <- read.delim(INFILE)
head(coord_df)

# Specify the order of some factors in coord_df for plotting later
coord_df$pop <- factor(coord_df$pop, levels = coord_df$pop)
coord_df$region <- factor(coord_df$region, levels = coord_df$region)
coord_df$location <- factor(coord_df$location, levels = coord_df$location)


# Specify some countries that you want in the map
MYREGIONS <- c("Japan", "Russia", "China", "North Korea", "South Korea", "Taiwan", 
               "Canada", "USA", "Mexico")

# Re-center the map on the Pacific Ocean
coord_df$long2 <- ifelse(coord_df$long < -25, coord_df$long + 360, coord_df$long)

mapWorld <- map_data('world', wrap = c(0, 360)) %>%
  filter(region %in% MYREGIONS)


# Set the ggplot theme
theme_set(
  theme(panel.background = element_rect(fill = "aliceblue"), 
        panel.grid.major = element_line(colour = NA),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(), 
        #axis.text = element_text(size = 12),
        #axis.title = element_text(size = 14),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))  
)

# Specify some limits for the plot
my_xlim <- c(135, 220)
my_ylim <- c(32,73)

# Specify a color palette
mypalette <- c("#d94801", # Japan - orange
               "#08306b", # Chukchi - deep blue
               "#08519c", "#4292c6", "#9ecae1", # Bering Sea - blues
               "#006d2c", "#31a354", "#74c476", "#bae4b3", # Aleutians - greens
               "#fcc5c0", "#fa9fb5", "#f768a1", # AK Pen - pinks
               "#dd3497", "#ae017e", "#7a0177", "#49006a") # GOA - purples

plotA <- ggplot() +
  geom_polygon(data = mapWorld, 
               aes(x = long, y = lat, group = group), 
               fill = " grey47", alpha = 0.5) +
  geom_point(data = coord_df, aes(x = long2, y = lat, 
                                  color = location), size = 4) +
  coord_map(xlim = my_xlim, ylim = my_ylim) +
  scale_color_manual(values = mypalette) +
  guides(color = guide_legend(title = "Collection Site")) +
  annotate("text", x = 140, y = 38, label = "Japan", size = 6) +
  annotate("text", x = 148.8, y = 55, label = "Sea of Okhotsk") +
  annotate("text", x = 156, y = 66, label = "Russia", size = 6) +
  annotate("text", x = 189, y = 71, label = "Chukchi Sea") +
  annotate("text", x = 182, y = 56.99, label = "Bering Sea") +
  annotate("text", x = 188, y = 51, label = "Aleutian Islands") +
  annotate("text", x = 216, y = 58.8, label = "Gulf of Alaska") +
  annotate("text", x = 214, y = 66, label = "Alaska", size = 6) +
  annotate("text", x = 182, y = 42, label = "Pacific Ocean", size = 7)

plotA

# Save plot to file

ggsave(OUTFILE, plot = plotA) 






