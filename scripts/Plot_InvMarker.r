library(tidyverse)
library(stringr)
library(ggpubr)
library(ggExtra)

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

    # Scatter plot of frequency by position
    PLOT <- ggplot(data_at, aes(x = Pos, y = Freq, col = Freq)) +
        geom_point() +
        # facet_grid(ID ~ .) +
        theme_bw() +
        theme(legend.position = "none") +
        geom_hline(
            data = medians_df,
            aes(yintercept = median_value),
            color = "red",
            linetype = "dashed",
            size = 1
        ) +
        ggtitle(id)

    # Combine the two plots
    plot_final <- ggMarginal(PLOT,
        type = "histogram",
        margins = "y",
        color = "black",
        fill = "blue",
        bins = 10
    )

    # Save the combined plot
    output_file <- paste0(output_dir, id, ".png")
    ggsave(file = output_file, plot_final, width = 6, height = 4)
    output_file <- paste0(output_dir, id, ".pdf")
    ggsave(file = output_file, plot_final, width = 6, height = 4)
}

# Create and save plots for each ID in the data
unique_ids <- unique(data$ID)
for (id in unique_ids) {
    create_and_save_plots(data, id, output_dir)
}
