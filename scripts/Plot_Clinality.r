library(ggpubr)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

INV <- args[1]
WD <- args[2]

setwd(WD)

## Get meta data
meta <- read.csv("data/meta.csv",
    header = TRUE
)
meta.sub <- meta %>%
    select(sampleId, continent, country, lat, long)
meta.sub$sampleId <- gsub("-", ".", meta.sub$sampleId)

InvFreq <- read.table(paste0("results/SNPs_", INV, "/", INV, ".af"), header = T)
InvFreq$sampleId <- gsub("-", ".", InvFreq$Sample)

InvFreq <- InvFreq %>%
    inner_join(meta.sub, by = "sampleId")

InvFreq.eur <- InvFreq %>%
    filter(continent == "Europe")

# Plot with ggplot2

# Logistic regression model
model <- glm(InvFreq.eur[[INV]] ~ InvFreq.eur[["lat"]], family = binomial)

# p-value
p_value <- summary(model)$coefficients[2, 4]
p.eur.lat <- ggplot(InvFreq.eur, aes_string(x = "lat", y = INV)) +
    geom_point() +
    geom_smooth(
        method = "glm",
        method.args = list(family = "binomial"),
        se = FALSE,
        col = "blue"
    ) +
    annotate("text",
        x = Inf, y = Inf, label = paste("p-Value = ", format(p_value, digits = 2)),
        hjust = 1.1, vjust = 2, size = 5, color = "black"
    ) +
    labs(
        title = paste0("Latitude vs. ", INV, " Europe"),
        y = paste0(INV, " Frequency"),
        x = "Latitude"
    ) +
    theme_bw()
p.eur.lat

# Logistic regression model
model <- glm(InvFreq.eur[[INV]] ~ InvFreq.eur[["long"]], family = binomial)

# p-value
p_value <- summary(model)$coefficients[2, 4]

p.eur.lon <- ggplot(InvFreq.eur, aes_string(x = "long", y = INV)) +
    geom_point() +
    geom_smooth(
        method = "glm",
        method.args = list(family = "binomial"),
        se = FALSE,
        col = "blue"
    ) +
    annotate("text",
        x = Inf, y = Inf, label = paste("p-Value = ", format(p_value, digits = 2)),
        hjust = 1.1, vjust = 2, size = 5, color = "black"
    ) +
    labs(
        title = paste0("Longitude vs. ", INV, " Europe"),
        y = paste0(INV, " Frequency"),
        x = "Longitude"
    ) +
    theme_bw()
p.eur.lon


InvFreq.na <- InvFreq %>%
    filter(continent == "North_America" & country != "Guadeloupe" & country != "Panama")

# Logistic regression model
model <- glm(InvFreq.na[[INV]] ~ InvFreq.na[["lat"]], family = binomial)

# p-value
p_value <- summary(model)$coefficients[2, 4]

# Plot with ggplot2
p.na.lat <- ggplot(InvFreq.na, aes_string(x = "lat", y = INV)) +
    geom_point() +
    geom_smooth(
        method = "glm",
        method.args = list(family = "binomial"),
        se = FALSE,
        col = "blue"
    ) +
    annotate("text",
        x = Inf, y = Inf, label = paste("p-Value = ", format(p_value, digits = 2)),
        hjust = 1.1, vjust = 2, size = 5, color = "black"
    ) +
    labs(
        title = paste0("Latitude vs. ", INV, " North_America"),
        y = paste0(INV, " Frequency"),
        x = "Latitude"
    ) +
    theme_bw()
p.na.lat

# Logistic regression model
model <- glm(InvFreq.na[[INV]] ~ InvFreq.na[["long"]], family = binomial)

# p-value
p_value <- summary(model)$coefficients[2, 4]


p.na.lon <- ggplot(InvFreq.na, aes_string(x = "long", y = INV)) +
    geom_point() +
    geom_smooth(
        method = "glm",
        method.args = list(family = "binomial"),
        se = FALSE,
        col = "blue"
    ) +
    annotate("text",
        x = Inf, y = Inf, label = paste("p-Value = ", format(p_value, digits = 2)),
        hjust = 1.1, vjust = 2, size = 5, color = "black"
    ) +
    labs(
        title = paste0("Longitude vs. ", INV, " North_America"),
        y = paste0(INV, " Frequency"),
        x = "Longitude"
    ) +
    theme_bw()
p.na.lon

PLOT <- ggarrange(p.eur.lat, p.na.lat, p.eur.lon, p.na.lon, nrow = 2, ncol = 2)
FILE <- paste0("results/SNPs_", INV, "/PCA-Clines_", INV, ".png")
ggsave(
    file = FILE,
    PLOT,
    width = 12,
    height = 6
)
