library(tidyverse)
library(ggpubr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
INV <- args[1]
WD <- args[2]

# Set working directory
setwd(WD)

# Load metadata and adjust sample IDs
meta <- read.csv("data/meta.csv", header = TRUE) %>%
    select(sampleId, continent, country, lat, long) %>%
    mutate(sampleId = gsub("-", ".", sampleId))

# Load inversion frequency data and adjust sample IDs
inv_freq <- read.table(paste0("results/SNPs_", INV, "/", INV, ".af"), header = TRUE) %>%
    mutate(sampleId = gsub("-", ".", Sample)) %>%
    inner_join(meta, by = "sampleId")

# Function to create plots
create_plot <- function(data, x_var, y_var, title_suffix, output_prefix) {
    model <- lm(asin(sqrt(data[[y_var]])) ~ data[[x_var]])
    p_value <- summary(model)$coefficients[2, 4]

    plot <- ggplot(data, aes_string(x = x_var, y = y_var)) +
        geom_point() +
        geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE, col = "blue") +
        annotate("text", x = Inf, y = Inf, label = paste("p-Value =", format(p_value, digits = 2)), hjust = 1.1, vjust = 2, size = 5, color = "black") +
        labs(title = paste0(x_var, " vs. ", INV, " ", title_suffix), y = paste0(INV, " Frequency"), x = x_var) +
        theme_bw()

    return(plot)
}

# Subset data for Europe and create plots
inv_freq_eur <- inv_freq %>% filter(continent == "Europe")
p_eur_lat <- create_plot(inv_freq_eur, "lat", INV, "Europe", "eur_lat")
p_eur_lon <- create_plot(inv_freq_eur, "long", INV, "Europe", "eur_lon")

# Subset data for North America (excluding Guadeloupe and Panama) and create plots
inv_freq_na <- inv_freq %>% filter(continent == "North_America" & country != "Guadeloupe" & country != "Panama")
p_na_lat <- create_plot(inv_freq_na, "lat", INV, "North America", "na_lat")
p_na_lon <- create_plot(inv_freq_na, "long", INV, "North America", "na_lon")

# Arrange plots into a grid and save the combined plot
combined_plot <- ggarrange(p_eur_lat, p_na_lat, p_eur_lon, p_na_lon, nrow = 2, ncol = 2)
output_file <- paste0("results/SNPs_", INV, "/Clines_", INV, ".png")
ggsave(file = output_file, plot = combined_plot, width = 12, height = 6)
output_file <- paste0("results/SNPs_", INV, "/Clines_", INV, ".pdf")
ggsave(file = output_file, plot = combined_plot, width = 12, height = 6)
