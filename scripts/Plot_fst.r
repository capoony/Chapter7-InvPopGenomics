library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)

INV <- args[1]
WD <- args[2]

setwd(WD)

DATA <- na.omit(read.table(paste0("results/SNPs_", INV, "/", INV, ".fst.weir.fst"),
    na.strings = "-nan",
    header = TRUE
))

DATA$WEIR_AND_COCKERHAM_FST[DATA$WEIR_AND_COCKERHAM_FST < 0] <- 0

DATA1 <- DATA %>%
    filter(CHROM %in% c("X", "2L", "2R", "3L", "3R", "4"))

PLOT <- ggplot(
    DATA1,
    aes(
        x = POS / 1000000,
        y = WEIR_AND_COCKERHAM_FST
    )
) +
    geom_point(col = rgb(0, 0, 0, 0.1), pch = 16) +
    facet_grid(. ~ CHROM,
        scales = "free_x",
        space = "free"
    ) +
    theme_bw() +
    xlab("Position (Mbp)")

ggsave(
    file = paset0("results/SNPs_", INV, "/", INV, ".fst.weir.fst.png"),
    PLOT,
    width = 16,
    height = 5
)
