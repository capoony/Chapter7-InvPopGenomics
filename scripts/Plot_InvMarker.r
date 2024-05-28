library(tidyverse)
library(stringr)
library(ggpubr)

# Define command line arguments
args <- commandArgs(trailingOnly = TRUE)
INV <- args[1]
WD <- args[2]

# Set working directory
setwd(WD)

# Load data
data_file <- paste0("results/SNPs_", INV, "/", INV, "_pos.af")
data <- read.table(data_file, header = TRUE)

# Create output directory if it doesn't exist
output_dir <- paste0("results/SNPs_", INV, "/", INV, "_plots/")
dir.create(output_dir, showWarnings = FALSE)

# Function to create and save plots for each ID
create_and_save_plots <- function(data, id, output_dir) {
    # Filter data for the given ID
    data_at <- data %>% filter(ID == id)

    # Calculate median frequency for the given ID
    medians_df <- data_at %>%
        group_by(ID) %>%
        summarize(median_value = median(Freq))

    # Histogram plot with median line
    plot1 <- ggplot(data_at, aes(x = Freq)) +
        geom_histogram(bins = 10, fill = "blue", color = "black") +
        facet_grid(. ~ ID) +
        theme_bw() +
        geom_vline(
            data = medians_df,
            aes(xintercept = median_value),
            color = "red",
            linetype = "dashed",
            size = 1
        )

    # Scatter plot of frequency by position
    plot2 <- ggplot(data_at, aes(x = Pos, y = Freq, col = Freq)) +
        geom_point() +
        facet_grid(ID ~ .) +
        theme_bw() +
        theme(legend.position = "none")

    # Combine the two plots
    plot_final <- ggarrange(plot1, plot2, nrow = 2)

    # Save the combined plot
    output_file <- paste0(output_dir, id, ".png")
    ggsave(file = output_file, plot_final, width = 6, height = 6)
}

# Create and save plots for each ID in the data
unique_ids <- unique(data$ID)
for (id in unique_ids) {
    create_and_save_plots(data, id, output_dir)
}
