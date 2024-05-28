library(tidyverse)

# Define command line arguments
args <- commandArgs(trailingOnly = TRUE)
INV <- args[1]
WD <- args[2]

# Set working directory
setwd(WD)

# Function to load and clean FST data
load_and_clean_data <- function(file_path) {
    data <- read.table(file_path, na.strings = "-nan", header = TRUE) %>%
        na.omit()
    data$WEIR_AND_COCKERHAM_FST[data$WEIR_AND_COCKERHAM_FST < 0] <- 0
    return(data)
}

# Function to filter data by chromosome
filter_by_chromosome <- function(data, chromosomes) {
    filtered_data <- data %>% filter(CHROM %in% chromosomes)
    return(filtered_data)
}

# Function to create FST plot
create_fst_plot <- function(data) {
    plot <- ggplot(data, aes(x = POS / 1000000, y = WEIR_AND_COCKERHAM_FST)) +
        geom_point(col = rgb(0, 0, 0, 0.1), pch = 16) +
        facet_grid(. ~ CHROM, scales = "free_x", space = "free") +
        theme_bw() +
        xlab("Position (Mbp)")
    return(plot)
}

# Function to save plot to file
save_plot <- function(plot, file_path, width = 16, height = 5) {
    ggsave(file = file_path, plot, width = width, height = height)
}

# Main script
file_path <- paste0("results/SNPs_", INV, "/", INV, ".fst.weir.fst")
data <- load_and_clean_data(file_path)
chromosomes <- c("X", "2L", "2R", "3L", "3R", "4")
filtered_data <- filter_by_chromosome(data, chromosomes)
plot <- create_fst_plot(filtered_data)
output_file_path <- paste0("results/SNPs_", INV, "/", INV, ".fst.weir.fst.png")
save_plot(plot, output_file_path)
