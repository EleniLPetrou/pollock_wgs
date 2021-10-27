# The purpose of this script is to plot the results of a pca analysis, using a 
# covariance matrix produced by pcangsd v1.02 and the bamlist used in angsd.


library(tidyverse)
library(RColorBrewer)


# Specify data directory and file names
DATADIR <- "~/pollock/results/pcangsd/samples580_maf0.05_miss0.3_noJPN20"
COV_FILE <- "results_pca_samples580_maf0.05_miss0.3_noJPN20.cov"
BAM_LIST <- "samples580_noJPN20_bam.filelist" 
METADATA <- "~/pollock/sample_metadata/collection_site_metadata.txt"

OUTFILE <- "pcangsd_samples580_maf0.05_miss0.3_noJPN20.pdf"

# Set working directory

setwd(DATADIR)

# Load the covariance matrix
cov_mat <- as.matrix(read.table(COV_FILE, header = F))

bam_df <- read.delim(BAM_LIST, header = F, col.names = "sample_id")
meta_df <- read.delim(METADATA, header = TRUE)
meta_df$pop <- factor(meta_df$pop, levels = meta_df$pop)

# Add a population information column to bam_df

bam_df <-  bam_df %>%
  separate(sample_id , c("pop", "id"), sep = "_", remove = FALSE, extra = "drop")

# Perform the pca on the covariance matrix using the eigen function.
mme.pca <- eigen(cov_mat)  

# Extract Eigenvectors 
eigenvectors = mme.pca$vectors 
head(eigenvectors)
eigen_df <- as.data.frame(eigenvectors)


# Combine PCA results with sampling location information
pca_df = cbind(bam_df, eigen_df)

plotting_df <- left_join(pca_df, meta_df, by = "pop" )

plotting_df$pop <- factor(plotting_df$pop, levels = meta_df$pop)
plotting_df$region <- factor(plotting_df$region, levels = unique(meta_df$region))


# Summary statistics about the PCA

pca.eigenval.sum = sum(mme.pca$values) #sum of eigenvalues
varPC1 <- (mme.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (mme.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (mme.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (mme.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4

PC1 <- round(varPC1, 2)
PC2 <- round(varPC2, 2)

### Plot the PCA results

mypalette <- c(
               "#08306b", # Chukchi - deep blue
               "#08519c", "#4292c6", "#9ecae1", # Bering Sea - blues
               "#006d2c", "#31a354", "#74c476", "#bae4b3", # Aleutians - greens
               "#fcc5c0", "#fa9fb5", "#f768a1", # AK Pen - pinks
               "#dd3497", "#ae017e", "#7a0177", "#49006a") # GOA - purples

pca_plot <- ggplot(data = plotting_df) +
  geom_point(aes(x = V1 , y = V2, color = interaction(region, pop, sep = " - ")), 
             alpha = 0.7, pch = 16, size = 2 ) +
  scale_color_manual(values = mypalette) +
  guides(color = guide_legend(title = "Sample")) +
  xlab(paste0("PC1: ", PC1, "% variation")) +
  ylab(paste0("PC2: ", PC2, "% variation")) +
  theme_classic()

pca_plot

# Save the plot
ggsave(OUTFILE, plot = pca_plot)
