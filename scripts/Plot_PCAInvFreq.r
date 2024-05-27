library(tidyverse)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)

INV <- args[1]
WD <- args[2]

setwd(WD)

PCA.inv.na <- read.table(paste0("results/SNPs_", INV, "/PCA_", INV, "_NorthAmerica_inside.txt"), header = T, sep = ",")
PCA.non.na <- read.table(paste0("results/SNPs_", INV, "/PCA_", INV, "_NorthAmerica_outside.txt"), header = T, sep = ",")
PCA.inv.eur <- read.table(paste0("results/SNPs_", INV, "/PCA_", INV, "_Europe_inside.txt"), header = T, sep = ",")
PCA.non.eur <- read.table(paste0("results/SNPs_", INV, "/PCA_", INV, "_Europe_outside.txt"), header = T, sep = ",")

InvFreq <- read.table(paste0("results/SNPs_", INV, "/", INV, ".af"), header = T)
InvFreq$sampleId <- gsub("-", ".", InvFreq$Sample)

PCA.inv.na <- PCA.inv.na %>%
    inner_join(InvFreq, by = "sampleId")
PCA.non.na <- PCA.non.na %>%
    inner_join(InvFreq, by = "sampleId")
PCA.inv.eur <- PCA.inv.eur %>%
    inner_join(InvFreq, by = "sampleId")
PCA.non.eur <- PCA.non.eur %>%
    inner_join(InvFreq, by = "sampleId")

# calculate linear model
model <- lm(PCA.inv.na[[INV]] ~ PCA.inv.na[["Dim.1"]])

# calculate R²
r_squared <- summary(model)$r.squared

# Plot with ggplot2
p1 <- ggplot(PCA.inv.na, aes_string(x = INV, y = "Dim.1")) +
    geom_point() +
    geom_smooth(method = "lm", col = "blue") +
    annotate("text",
        x = Inf, y = Inf, label = paste("R² = ", round(r_squared, 2)),
        hjust = 1.1, vjust = 2, size = 5, color = "black"
    ) +
    labs(
        title = paste0("Inside ", INV, " North America"),
        x = paste0(INV, " Frequency"),
        y = "Dim.1"
    ) +
    theme_bw()
p1

# calculate linear model
model <- lm(PCA.non.na[[INV]] ~ PCA.non.na[["Dim.1"]])

# calculate R²
r_squared <- summary(model)$r.squared

# Plot with ggplot2
p2 <- ggplot(PCA.non.na, aes_string(x = INV, y = "Dim.1")) +
    geom_point() +
    geom_smooth(method = "lm", col = "blue") +
    annotate("text",
        x = Inf, y = Inf, label = paste("R² = ", round(r_squared, 2)),
        hjust = 1.1, vjust = 2, size = 5, color = "black"
    ) +
    labs(
        title = paste0("Outside ", INV, " North America"),
        x = paste0(INV, " Frequency"),
        y = "Dim.1"
    ) +
    theme_bw()
p2

# calculate linear model
model <- lm(PCA.inv.eur[[INV]] ~ PCA.inv.eur[["Dim.1"]])

# calculate R²
r_squared <- summary(model)$r.squared

# Plot with ggplot2
p3 <- ggplot(PCA.inv.eur, aes_string(x = INV, y = "Dim.1")) +
    geom_point() +
    geom_smooth(method = "lm", col = "blue") +
    annotate("text",
        x = Inf, y = Inf, label = paste("R² = ", round(r_squared, 2)),
        hjust = 1.1, vjust = 2, size = 5, color = "black"
    ) +
    labs(
        title = paste0("Inside ", INV, " Europe"),
        x = paste0(INV, " Frequency"),
        y = "Dim.1"
    ) +
    theme_bw()
p3

# calculate linear model
model <- lm(PCA.non.eur[[INV]] ~ PCA.non.eur[["Dim.1"]])

# calculate R²
r_squared <- summary(model)$r.squared

# Plot with ggplot2
p4 <- ggplot(PCA.non.eur, aes_string(x = INV, y = "Dim.1")) +
    geom_point() +
    geom_smooth(method = "lm", col = "blue") +
    annotate("text",
        x = Inf, y = Inf, label = paste("R² = ", round(r_squared, 2)),
        hjust = 1.1, vjust = 2, size = 5, color = "black"
    ) +
    labs(
        title = paste0("Outside ", INV, " Europe"),
        x = paste0(INV, " Frequency"),
        y = "Dim.1"
    ) +
    theme_bw()
p4

# save combined plot
PLOT <- ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
FILE <- paste0("results/SNPs_", INV, "/PCA-InvFreq_", INV, ".png")
ggsave(
    file = FILE,
    PLOT,
    width = 12,
    height = 6
)
