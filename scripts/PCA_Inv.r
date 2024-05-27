library(tidyverse)
library(factoextra)
library(FactoMineR)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)

INV <- args[1]
Chr <- args[2]
Start <- args[3]
End <- args[4]
WD <- args[5]

setwd(WD)

## Get meta data
meta <- read.csv("data/meta.csv",
    header = TRUE
)
meta.sub <- meta %>%
    select(sampleId, country, province)
meta.sub$sampleId <- gsub("-", ".", meta.sub$sampleId)


for (i in c("Europe", "NorthAmerica")) {
    ### read frequency data
    DATA <- read.table(gzfile(paste0("results/SNPs/", i, "_freq.csv.gz")),
        header = TRUE,
        comment.char = ""
    )
    ### subset to SNPs inside the inversion
    DATA.inv <- as.data.frame(t(DATA %>%
        filter(X.CHROM == Chr & POS > as.numeric(Start) & POS < as.numeric(End)) %>%
        select(3:ncol(DATA))))
    ### get sampleIds and link to metadata
    DATA.inv$sampleId <- rownames(DATA.inv)
    DATA.inv <- DATA.inv %>%
        inner_join(meta.sub, by = "sampleId") %>%
        filter(country != "Panama" & country != "Guadeloupe")
    ### do PCA

    PCA.inv <- PCA(DATA.inv[, 1:(ncol(DATA.inv) - 3)],
        graph = FALSE
    )
    ### make output table
    write.table(
        file = paste0("results/SNPs_", INV, "/PCA_", INV, "_", i, "_inside.txt"),
        cbind(DATA.inv[, (ncol(DATA.inv) - 2):ncol(DATA.inv)], PCA.inv$ind$coord),
        sep = ",",
        quote = FALSE,
        row.names = FALSE
    )
    ### make scatterplot
    if (i == "Europe") {
        COLOR <- DATA.inv$country
    } else {
        COLOR <- DATA.inv$province
    }
    PLOT.inv <- fviz_pca_ind(PCA.inv,
        col.ind = COLOR,
        pointshape = 16,
        pointsize = 3,
        alpha = 0.7,
        mean.point = FALSE,
        label = COLOR,
        repel = TRUE
    ) + theme_bw() + ggtitle(paste0("PCA - ", i, " inside ", INV))

    ### repeat for SNPs outside the inversion
    DATA.non <- as.data.frame(t(DATA %>%
        filter(!(X.CHROM == Chr & POS > as.numeric(Start) & POS < as.numeric(End))) %>%
        select(3:ncol(DATA))))
    DATA.non$sampleId <- rownames(DATA.non)
    DATA.non <- DATA.non %>%
        inner_join(meta.sub, by = "sampleId")
    PCA.non <- PCA(DATA.non[, 1:(ncol(DATA.non) - 3)],
        graph = FALSE
    )
    ### make output table
    write.table(
        file = paste0("results/SNPs_", INV, "/PCA_", INV, "_", i, "_outside.txt"),
        cbind(DATA.non[, (ncol(DATA.non) - 2):ncol(DATA.non)], PCA.non$ind$coord),
        sep = ",",
        quote = FALSE,
        row.names = FALSE
    )
    if (i == "Europe") {
        COLOR <- DATA.non$country
    } else {
        COLOR <- DATA.non$province
    }
    PLOT.non <- fviz_pca_ind(PCA.non,
        col.ind = COLOR,
        pointshape = 16,
        pointsize = 3,
        alpha = 0.7,
        mean.point = FALSE,
        label = COLOR,
        repel = TRUE
    ) + theme_bw() + ggtitle(paste0("PCA - ", i, " outside ", INV))

    PLOT <- ggarrange(PLOT.inv, PLOT.non, common.legend = T)
    FILE <- paste0("results/SNPs_", INV, "/PCA_", INV, "_", i, ".png")
    ggsave(
        file = FILE,
        PLOT,
        width = 12,
        height = 6
    )
}
