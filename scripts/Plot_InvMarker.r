library(tidyverse)
library(stringr)
library(ggpubr)
args <- commandArgs(trailingOnly = TRUE)

INV <- args[1]
WD <- args[2]

setwd(WD)
DATA <- read.table(paste0("results/SNPs_", INV, "/", INV, "_pos.af"),
    header = TRUE
)

dir.create(paste0("results/SNPs_", INV, "/", INV, "_plots/"))

for (i in DATA$ID) {
    DATA.AT <- DATA %>%
        filter(ID == i)

    medians_df <- DATA.AT %>%
        group_by(ID) %>%
        summarize(median_value = median(Freq))

    PLOT <- ggplot(DATA.AT, aes(x = Freq)) +
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
    # PLOT

    PLOT2 <- ggplot(DATA.AT, aes(
        x = Pos,
        y = Freq,
        col = Freq
    )) +
        geom_point() +
        facet_grid(ID ~ .) +
        theme_bw() +
        theme(legend.position = "none")
    # PLOT2

    PLOT.final <- ggarrange(PLOT, PLOT2, nrow = 2)
    FILE <- paste0("results/SNPs_", INV, "/", INV, "_plots/", i, ".png")
    ggsave(
        file = FILE,
        PLOT.final,
        width = 6,
        height = 6
    )
}
