library(tidyverse)
library(LEA)
library(ggpubr)
library(geodata)

args <- commandArgs(trailingOnly = TRUE)

INV <- args[1]
Chr <- args[2]
Start <- as.numeric(args[3])
End <- as.numeric(args[4])
WD <- args[5]

setwd(WD)

# Functions
get_meta_data <- function() {
    meta <- read.csv("data/meta.csv", header = TRUE)
    meta.sub <- meta %>% select(sampleId, continent, country, province, lat, long)
    meta.sub$sampleId <- gsub("-", ".", meta.sub$sampleId)
    return(meta.sub)
}

get_inv_freq <- function(meta.sub) {
    InvFreq <- read.table(paste0("results/SNPs_", INV, "/", INV, ".af"), header = TRUE)
    InvFreq$sampleId <- gsub("-", ".", InvFreq$Sample)
    InvFreq <- InvFreq %>% inner_join(meta.sub, by = "sampleId")
    return(InvFreq)
}

get_worldclim_data <- function(meta.sub) {
    biod <- worldclim_global(var = "bio", 2.5, "data")
    bio <- raster::extract(biod, cbind(meta.sub$long, meta.sub$lat))
    bio.sub <- as.data.frame(bio) %>% select(`wc2.1_2.5m_bio_1`, `wc2.1_2.5m_bio_12`)
    colnames(bio.sub) <- c("AvTemp", "AvPrec")
    return(bio.sub)
}

create_plot <- function(data, var, title, colors, ggtitle_text, xlab_text, ylab_text) {
    world_coordinates <- map_data("world")
    WORLD <- ggplot(data, aes(x = long, y = lat, col = !!sym(var))) +
        geom_map(
            data = world_coordinates, map = world_coordinates,
            aes(long, lat, map_id = region),
            color = "black", fill = "lightgrey"
        ) +
        geom_point() +
        ggtitle(ggtitle_text) +
        scale_colour_gradientn(colours = colors) +
        theme_bw() +
        theme(legend.position = "bottom")
    return(WORLD)
}

plot_correlation <- function(data, var, axis, title) {
    ggplot(data, aes_string(x = axis, y = var)) +
        geom_point(col = rgb(0, 0, 0, 0.1)) +
        ggtitle(title) +
        geom_smooth(method = "lm", col = "blue") +
        theme_bw()
}

save_plot <- function(file, plot, width, height) {
    ggsave(file, plot, width = width, height = height)
}

run_lfmm <- function(data, env, K_range = 1:10) {
    project <- snmf(data, K = K_range, entropy = TRUE, repetitions = 2, project = "new")
    CE <- sapply(K_range, function(K) mean(cross.entropy(project, K = K)))
    K <- which.min(CE)
    mod.lfmm2 <- lfmm2(data, env, K = K)
    return(mod.lfmm2)
}

process_continent_data <- function(continent, inv_freq, meta_sub, chr, start, end) {
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
        filter(country != "Panama" & country != "Guadeloupe" & sampleId != "US_Lou_Bat_0_2013-09-15")

    world_temp_plot <- create_plot(data.af, "AvTemp", "Average Temperature (°C)", terrain.colors(10), "Average Temperature (°C)", "Longitude", "Latitude")
    corr_lat_temp_plot <- plot_correlation(data.af, "AvTemp", "lat", "Corr. Latitude & Average Temp")
    corr_lon_temp_plot <- plot_correlation(data.af, "AvTemp", "long", "Corr. Longitude & Average Temp")
    plot_temp <- ggarrange(world_temp_plot, ggarrange(corr_lat_temp_plot, corr_lon_temp_plot, ncol = 2), nrow = 2, heights = c(2, 1))
    save_plot(paste0("results/SNPs_", INV, "/LFMM_", INV, "_", continent, "_WorldTemp.png"), plot_temp, 10, 10)

    world_prec_plot <- create_plot(data.af, "AvPrec", "Average Precipitation (mm)", terrain.colors(10), "Average Precipitation (mm)", "Longitude", "Latitude")
    corr_lat_prec_plot <- plot_correlation(data.af, "AvPrec", "lat", "Corr. Latitude & Average Prec")
    corr_lon_prec_plot <- plot_correlation(data.af, "AvPrec", "long", "Corr. Longitude & Average Prec")
    plot_prec <- ggarrange(world_prec_plot, ggarrange(corr_lat_prec_plot, corr_lon_prec_plot, ncol = 2), nrow = 2, heights = c(2, 1))
    save_plot(paste0("results/SNPs_", INV, "/LFMM_", INV, "_", continent, "_WorldPrecipitation.png"), plot_prec, 10, 10)

    rownames(data.af) <- data.af$sampleId
    data.af.meta <- data.af %>% select(sampleId, continent, country, lat, long, AvTemp, AvPrec)
    data.af.lfmm <- data.af %>% select(-Sample, -sampleId, -continent, -country, -province, -lat, -long, -AvTemp, -AvPrec)
    data.af.lfmm <- data.af.lfmm[, colSums(data.af.lfmm) != 0]
    data.env <- data.af.meta %>% select(AvTemp, AvPrec)

    write.lfmm(data.af.lfmm, paste0("results/SNPs_", INV, "/LFMM_", continent, "/", continent, ".lfmm"))

    mod.lfmm2 <- run_lfmm(paste0("results/SNPs_", INV, "/LFMM_", continent, "/", continent, ".lfmm"), data.env)

    png(paste0("results/SNPs_", INV, "/LFMM_", continent, "/", continent, "_LatentFactors.png"), width = 500, height = 500)
    plot(mod.lfmm2@U, col = "blue", pch = 19, cex = 1.2)
    text(mod.lfmm2@U, data.af.meta$sampleId, cex = 0.65, pos = 3, col = "red")
    dev.off()

    pv <- lfmm2.test(object = mod.lfmm2, input = data.af.lfmm, env = data.env, linear = TRUE)
    pval <- cbind(ID = colnames(data.af.lfmm), t(-log10(pv$pvalues)))

    inversion <- pval[1, ]
    hline_data <- data.frame(
        yat = as.numeric(inversion[2]),
        yap = as.numeric(inversion[3]),
        xmin = as.numeric(start) / 1000000,
        xmax = as.numeric(end) / 1000000,
        Chrom = c(chr)
    )

    others <- pval[2:nrow(pval), ] %>%
        as.data.frame() %>%
        separate(ID, into = c("Chrom", "Pos"), sep = ":")

    write.table(pval,
        file = paste0("results/SNPs_", INV, "/LFMM_", INV, "_", continent, ".txt"),
        quote = FALSE,
        row.names = FALSE
    )

    plot_list <- list(
        temp = ggplot(others, aes(x = as.numeric(Pos) / 1000000, y = as.numeric(AvTemp))) +
            geom_point(col = rgb(0, 0, 0, 0.1), pch = 16) +
            facet_grid(. ~ Chrom, scales = "free_x", space = "free") +
            theme_bw() +
            geom_hline(yintercept = -log10(0.05 / (nrow(others) + 1)), colour = "blue") +
            geom_segment(aes(x = xmin, y = yat, xend = xmax, yend = yat, colour = "segment"), data = hline_data) +
            xlab("Position (Mbp)") +
            ylab("-log10(p-value)") +
            ggtitle(paste0("LFMM: Temp ", INV, " for ", continent)) +
            theme(legend.position = "none"),
        prec = ggplot(others, aes(x = as.numeric(Pos) / 1000000, y = as.numeric(AvPrec))) +
            geom_point(col = rgb(0, 0, 0, 0.1), pch = 16) +
            facet_grid(. ~ Chrom, scales = "free_x", space = "free") +
            theme_bw() +
            geom_hline(yintercept = -log10(0.05 / (nrow(others) + 1)), colour = "blue") +
            geom_segment(aes(x = xmin, y = yap, xend = xmax, yend = yap, colour = "segment"), data = hline_data) +
            xlab("Position (Mbp)") +
            ylab("-log10(p-value)") +
            ggtitle(paste0("LFMM: Precipitation ", INV, " for ", continent)) +
            theme(legend.position = "none")
    )

    return(plot_list)
}

main <- function() {
    meta.sub <- get_meta_data()
    inv_freq <- get_inv_freq(meta.sub)
    bio.sub <- get_worldclim_data(meta.sub)
    inv_freq <- cbind(inv_freq, bio.sub)

    plot_list <- list(
        temp = list(),
        prec = list()
    )

    for (continent in c("Europe", "NorthAmerica")) {
        plot_continent <- process_continent_data(continent, inv_freq, meta.sub, Chr, Start, End)
        plot_list$temp[[continent]] <- plot_continent$temp
        plot_list$prec[[continent]] <- plot_continent$prec
    }

    combined_plots <- list(
        temp = ggarrange(plotlist = plot_list$temp, nrow = 2),
        prec = ggarrange(plotlist = plot_list$prec, nrow = 2)
    )

    save_plot(paste0("results/SNPs_", INV, "/LFMM_", INV, "_AvTemp.png"), combined_plots$temp, 10, 8)
    save_plot(paste0("results/SNPs_", INV, "/LFMM_", INV, "_AvPrec.png"), combined_plots$prec, 10, 8)
}

# Run main function
main()
