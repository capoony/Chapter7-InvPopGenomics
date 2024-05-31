library(tidyverse)

# Define command line arguments
args <- commandArgs(trailingOnly = TRUE)
INV <- args[1]
Chr <- args[2]
Start <- args[3]
End <- args[4]
WD <- args[5]

# Set working directory
setwd(WD)

# Function to load and clean FST data
load_and_clean_data <- function(file_path) {
    data <- read.table(file_path, na.strings = "-nan", header = TRUE) %>%
        na.omit()
    data$PI[data$PI < 0] <- 0
    return(data)
}

# Function to filter data by chromosome
filter_by_chromosome <- function(data, chromosomes) {
    filtered_data <- data %>% filter(CHROM %in% chromosomes)
    return(filtered_data)
}

# Create data frame for the rectangle (highlight region) in the plot
rect_data <- data.frame(
    xmin = as.numeric(Start) / 1000000,
    xmax = as.numeric(End) / 1000000,
    ymin = 0,
    ymax = Inf,
    CHROM = Chr
)

# Function to create FST plot
create_fst_plot <- function(data) {
    plot <- ggplot(data, aes(x = BIN_START / 1000000, y = PI, col = Type)) +
        geom_line() +
        # geom_point(col = rgb(0, 0, 0, 0.1), pch = 16, size = 0.1) +
        facet_grid(. ~ CHROM, scales = "free_x", space = "free") +
        theme_bw() +
        geom_rect(
            data = rect_data,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "blue", alpha = 0.2, inherit.aes = FALSE
        ) +
        xlab("Position (Mbp)")
    return(plot)
}

# Function to save plot to file
save_plot <- function(plot, file_path, width = 16, height = 5) {
    ggsave(file = file_path, plot, width = width, height = height)
}

# Main script
file_path <- paste0("results/SNPs_", INV, "/", INV, "_pi.tsv")
data <- load_and_clean_data(file_path)
chromosomes <- c("X", "2L", "2R", "3L", "3R", "4")
filtered_data <- filter_by_chromosome(data, chromosomes)
plot <- create_fst_plot(filtered_data)
output_file_path <- paste0("results/SNPs_", INV, "/", INV, "_pi.png")
save_plot(plot, output_file_path)
