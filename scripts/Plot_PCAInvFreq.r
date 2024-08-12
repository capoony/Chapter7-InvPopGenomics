# Load necessary libraries
library(tidyverse)
library(ggpubr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
INV <- args[1]
WD <- args[2]

# Set working directory
setwd(WD)

# Read the PCA and inversion frequency data
PCA.inv.na <- read.table(paste0("results/SNPs_", INV, "/PCA_", INV, "_NorthAmerica_inside.txt"), header = TRUE, sep = ",")
PCA.non.na <- read.table(paste0("results/SNPs_", INV, "/PCA_", INV, "_NorthAmerica_outside.txt"), header = TRUE, sep = ",")
PCA.inv.eur <- read.table(paste0("results/SNPs_", INV, "/PCA_", INV, "_Europe_inside.txt"), header = TRUE, sep = ",")
PCA.non.eur <- read.table(paste0("results/SNPs_", INV, "/PCA_", INV, "_Europe_outside.txt"), header = TRUE, sep = ",")

# Read the inversion frequency data and clean sample IDs
InvFreq <- read.table(paste0("results/SNPs_", INV, "/", INV, ".af"), header = TRUE)
InvFreq$sampleId <- gsub("-", ".", InvFreq$Sample)

# Join the PCA data with inversion frequency data
PCA.inv.na <- PCA.inv.na %>% inner_join(InvFreq, by = "sampleId")
PCA.non.na <- PCA.non.na %>% inner_join(InvFreq, by = "sampleId")
PCA.inv.eur <- PCA.inv.eur %>% inner_join(InvFreq, by = "sampleId")
PCA.non.eur <- PCA.non.eur %>% inner_join(InvFreq, by = "sampleId")

# Function to create plot and calculate R²
create_plot <- function(data, INV, title_suffix) {
    model <- lm(data[[INV]] ~ data[["Dim.1"]])
    r_squared <- summary(model)$r.squared

    plot <- ggplot(data, aes_string(x = INV, y = "Dim.1")) +
        geom_point() +
        geom_smooth(method = "lm", col = "blue") +
        annotate("text",
            x = Inf, y = Inf, label = paste("R² = ", round(r_squared, 2)),
            hjust = 1.1, vjust = 2, size = 5, color = "black"
        ) +
        labs(
            title = paste0(title_suffix),
            x = paste0(INV, " Frequency"),
            y = "Dim.1"
        ) +
        theme_bw()

    return(plot)
}

# Create plots for each condition
p1 <- create_plot(PCA.inv.na, INV, paste0("Inside ", INV, " North America"))
p2 <- create_plot(PCA.non.na, INV, paste0("Outside ", INV, " North America"))
p3 <- create_plot(PCA.inv.eur, INV, paste0("Inside ", INV, " Europe"))
p4 <- create_plot(PCA.non.eur, INV, paste0("Outside ", INV, " Europe"))

# Combine and save the plots
PLOT <- ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
FILE <- paste0("results/SNPs_", INV, "/PCA-InvFreq_", INV, ".png")
ggsave(
    file = FILE,
    plot = PLOT,
    width = 12,
    height = 6
)
FILE <- paste0("results/SNPs_", INV, "/PCA-InvFreq_", INV, ".pdf")
ggsave(
    file = FILE,
    plot = PLOT,
    width = 12,
    height = 6
)
