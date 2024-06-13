library(tidyverse)
library(LEA)
library(ggpubr)
library(geodata)

### get paramters from commandline
args <- commandArgs(trailingOnly = TRUE)

INV <- args[1]
Chr <- args[2]
Start <- as.numeric(args[3])
End <- as.numeric(args[4])
WD <- args[5]

### ste working directory
setwd(WD)

## function to find the elbow point
find_elbow <- function(eigenvalues) {
    n <- length(eigenvalues)
    x <- 1:n
    y <- eigenvalues

    # Normalize the data
    x_norm <- (x - min(x)) / (max(x) - min(x))
    y_norm <- (y - min(y)) / (max(y) - min(y))

    # Find the distances to the line connecting the first and last points
    line <- lm(y_norm ~ x_norm)
    a <- line$coefficients[2]
    b <- -1
    c <- line$coefficients[1]

    distances <- abs(a * x_norm + b * y_norm + c) / sqrt(a^2 + b^2)
    elbow <- which.max(distances)

    return(elbow)
}

# Functions
get_meta_data <- function() {
    ### get metadata-table
    meta <- read.csv("data/meta.csv", header = TRUE)
    meta.sub <- meta %>% select(sampleId, continent, country, province, lat, long)
    meta.sub$sampleId <- gsub("-", ".", meta.sub$sampleId)
    return(meta.sub)
}

get_inv_freq <- function(meta.sub) {
    ### get estimated inversion frequencies
    InvFreq <- read.table(paste0("results/SNPs_", INV, "/", INV, ".af"), header = TRUE)
    InvFreq$sampleId <- gsub("-", ".", InvFreq$Sample)
    ### combine with metadata based on sampleId
    InvFreq <- InvFreq %>% inner_join(meta.sub, by = "sampleId")
    return(InvFreq)
}


plot_correlation <- function(data, var, axis, title) {
    ### plot correlations between climatic and geographic variables
    ggplot(data, aes_string(x = axis, y = var)) +
        geom_point(col = rgb(0, 0, 0, 0.1)) +
        ggtitle(title) +
        geom_smooth(method = "lm", col = "blue") +
        theme_bw()
}

save_plot <- function(file, plot, width, height) {
    ggsave(file, plot, width = width, height = height)
}

run_lfmm <- function(data, env, K) {
    mod.lfmm2 <- lfmm2(data, env, K = K)
    return(mod.lfmm2)
}

process_continent_data <- function(continent, inv_freq, meta_sub, chr, start, end) {
    ### do the magic
    dir.create(paste0("results/SNPs_", INV, "/LFMM_", continent))
    data <- read.table(gzfile(paste0("results/SNPs/", continent, "_freq.csv.gz")), header = TRUE, comment.char = "")
    data.pos <- paste(data[, 1], data[, 2], sep = ":")
    rownames(data) <- data.pos

    data.af <- data %>%
        filter(!(X.CHROM == chr & POS > start & POS < end)) %>%
        select(3:ncol(data)) %>%
        t() %>%
        as.data.frame() %>%
        select(where(~ sum(.) != 0))

    data.af$sampleId <- rownames(data.af)
    data.af <- na.omit(inv_freq %>% inner_join(data.af, by = "sampleId")) %>%
        filter(country != "Panama" & country != "Guadeloupe" & sampleId != "US_Lou_Bat_0_2013.09.15")

    rownames(data.af) <- data.af$sampleId
    data.af.meta <- data.af %>% select(sampleId, continent, country, lat, long)
    data.af.lfmm <- data.af %>% select(-Sample, -sampleId, -continent, -country, -province, -lat, -long)
    data.af.lfmm <- data.af.lfmm[, colSums(data.af.lfmm) != 0]
    data.env <- data.af.meta %>% select(lat, long)

    pca_result <- prcomp(data.af.lfmm)
    eigenvalues <- pca_result$sdev^2

    explained_variance <- eigenvalues / sum(eigenvalues)

    # Calculate cumulative explained variance
    cumulative_variance <- cumsum(explained_variance)

    # Find the number of components that explain at least 90% of the variance
    threshold <- 0.25
    optimal_k <- which(cumulative_variance >= threshold)[1]

    write.lfmm(data.af.lfmm, paste0("results/SNPs_", INV, "/LFMM_", continent, "/", continent, ".lfmm"))

    mod.lfmm2 <- run_lfmm(paste0("results/SNPs_", INV, "/LFMM_", continent, "/", continent, ".lfmm"), data.env, optimal_k)

    png(paste0("results/SNPs_", INV, "/LFMM_", continent, "/", continent, "_LatentFactors.png"), width = 500, height = 500)
    plot(mod.lfmm2@U, col = "blue", pch = 19, cex = 1.2)
    text(mod.lfmm2@U, data.af.meta$sampleId, cex = 0.65, pos = 3, col = "red")
    dev.off()

    scree_data <- data.frame(Principal_Component = seq_along(eigenvalues), Eigenvalue = eigenvalues)

    SCREE <- ggplot(scree_data, aes(x = Principal_Component, y = Eigenvalue)) +
        geom_point() +
        geom_line() +
        geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
        labs(title = "Scree Plot with Elbow Point", x = "Principal Component", y = "Eigenvalue") +
        theme_minimal() +
        ggtitle(paste0("Optimal K: ", optimal_k))

    ggsave(
        file = paste0("results/SNPs_", INV, "/LFMM_", continent, "/", continent, "_Scree.png"),
        SCREE
    )
    pv <- lfmm2.test(object = mod.lfmm2, input = data.af.lfmm, env = data.env, linear = TRUE)
    pval <- cbind(ID = colnames(data.af.lfmm), t(-log10(pv$pvalues)))

    inversion <- pval[1, ]
    hline_data <- data.frame(
        ylat = as.numeric(inversion[2]),
        ylon = as.numeric(inversion[3]),
        xmin = as.numeric(start) / 1000000,
        xmax = as.numeric(end) / 1000000,
        Chrom = c(chr)
    )

    others <- pval[2:nrow(pval), ] %>%
        as.data.frame() %>%
        separate(ID, into = c("Chrom", "Pos"), sep = ":")

    plot_list <- list(
        lat = ggplot(others, aes(x = as.numeric(Pos) / 1000000, y = as.numeric(lat))) +
            geom_point(col = rgb(0, 0, 0, 0.1), pch = 16) +
            facet_grid(. ~ Chrom, scales = "free_x", space = "free") +
            theme_bw() +
            geom_hline(yintercept = -log10(0.05 / (nrow(others) + 1)), colour = "blue", linetype = "dashed") +
            geom_segment(aes(x = xmin, y = ylat, xend = xmax, yend = ylat, colour = "segment"), data = hline_data) +
            xlab("Position (Mbp)") +
            ylab("-log10(p-value)") +
            ggtitle(paste0("LFMM: Latitude ", INV, " for ", continent)) +
            theme(legend.position = "none"),
        lon = ggplot(others, aes(x = as.numeric(Pos) / 1000000, y = as.numeric(long))) +
            geom_point(col = rgb(0, 0, 0, 0.1), pch = 16) +
            facet_grid(. ~ Chrom, scales = "free_x", space = "free") +
            theme_bw() +
            geom_hline(yintercept = -log10(0.05 / (nrow(others) + 1)), colour = "blue", linetype = "dashed") +
            geom_segment(aes(x = xmin, y = ylon, xend = xmax, yend = ylon, colour = "segment"), data = hline_data) +
            xlab("Position (Mbp)") +
            ylab("-log10(p-value)") +
            ggtitle(paste0("LFMM: Longitude ", INV, " for ", continent)) +
            theme(legend.position = "none")
    )

    return(plot_list)
}

main <- function() {
    meta.sub <- get_meta_data()
    inv_freq <- get_inv_freq(meta.sub)

    plot_list <- list(
        lat = list(),
        lon = list()
    )

    for (continent in c("Europe", "NorthAmerica")) {
        plot_continent <- process_continent_data(continent, inv_freq, meta.sub, Chr, Start, End)
        plot_list$lat[[continent]] <- plot_continent$lat
        plot_list$lon[[continent]] <- plot_continent$lon
    }

    combined_plots <- list(
        lat = ggarrange(plotlist = plot_list$lat, nrow = 2),
        lon = ggarrange(plotlist = plot_list$lon, nrow = 2)
    )

    save_plot(paste0("results/SNPs_", INV, "/LFMM_", INV, "_Latitude.png"), combined_plots$lat, 10, 6)
    save_plot(paste0("results/SNPs_", INV, "/LFMM_", INV, "_Longitude.png"), combined_plots$lon, 10, 6)
}

# Run main function
main()
