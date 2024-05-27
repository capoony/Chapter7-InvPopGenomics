library(tidyverse)
# devtools::install_github("bcm-uga/LEA")
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

## Get meta data
meta <- read.csv("data/meta.csv",
    header = TRUE
)
meta.sub <- meta %>%
    select(sampleId, continent, country, province, lat, long)
meta.sub$sampleId <- gsub("-", ".", meta.sub$sampleId)

InvFreq <- read.table(paste0("results/SNPs_", INV, "/", INV, ".af"), header = T)
InvFreq$sampleId <- gsub("-", ".", InvFreq$Sample)

InvFreq <- InvFreq %>%
    inner_join(meta.sub, by = "sampleId")


# get WorldClim data
biod <- worldclim_global(var = "bio", 2.5, "data")

# extact for each coordinate bio clim variables
bio <- raster::extract(biod, cbind(meta.sub$long, meta.sub$lat))
bio.sub <- as.data.frame(bio) %>%
    select(`wc2.1_2.5m_bio_1`, `wc2.1_2.5m_bio_12`)

colnames(bio.sub) <- c("AvTemp", "AvPrec")

InvFreq <- cbind(InvFreq, bio.sub)

PlotList.lat <- list()
PlotList.lon <- list()
PlotList.at <- list()
PlotList.ap <- list()

for (i in c("Europe", "NorthAmerica")) {
    dir.create(paste0("results/SNPs_", INV, "/LFMM_", i))
    ### read frequency data
    DATA <- read.table(gzfile(paste0("results/SNPs/", i, "_freq.csv.gz")),
        header = TRUE,
        comment.char = ""
    )
    ## get positions
    DATA.pos <- paste(DATA[, 1], DATA[, 2], sep = ":")
    rownames(DATA) <- DATA.pos

    ### subset to SNPs inside the inversion
    DATA.af <- as.data.frame(t(DATA %>%
        filter(!(X.CHROM == Chr & POS > Start & POS < End)) %>%
        select(3:ncol(DATA)))) %>%
        select(where(~ sum(.) != 0))

    ### get sampleIds and link to metadata
    DATA.af$sampleId <- rownames(DATA.af)
    DATA.af <- na.omit(InvFreq %>%
        inner_join(DATA.af, by = "sampleId"))

    ###
    world_coordinates <- map_data("world")
    ggplot() +
        # geom_map() function takes world coordinates
        # as input to plot world map
        geom_map(
            data = world_coordinates, map = world_coordinates,
            aes(long, lat, map_id = region),
            color = "black",
            fill = "lightgrey",
        ) +
        theme_bw() +
        rownames(DATA.af) <- DATA.af$sampleId
    DATA.af.meta <- DATA.af %>%
        select(sampleId, continent, country, lat, long, AvTemp, AvPrec)
    DATA.af.lfmm <- DATA.af %>%
        select(-Sample, -sampleId, -continent, -country, -province, -lat, -long, -AvTemp, -AvPrec)
    DATA.af.lfmm <- DATA.af.lfmm[, colSums(DATA.af.lfmm) != 0]
    DATA.env <- DATA.af.meta %>%
        select(lat, long, AvTemp, AvPrec)
    write.lfmm(
        DATA.af.lfmm,
        paste0("results/SNPs_", INV, "/LFMM_", i, "/", i, ".lfmm")
    )

    ### get optimal numbers of K (more repetitions are better)
    project <- snmf(paste0("results/SNPs_", INV, "/LFMM_", i, "/", i, ".lfmm"),
        K = 1:10,
        entropy = TRUE,
        repetitions = 2,
        project = "new"
    )

    ### get cross.entropy scores and average
    CE <- c()
    for (j in seq(1, 10, 1)) {
        CE <- cbind(CE, mean(cross.entropy(project, K = j)))
    }
    K <- which.min(CE)

    ### run LFMM2
    mod.lfmm2 <- lfmm2(DATA.af.lfmm,
        DATA.env,
        K = K
    )

    ### Plot Latent Factors
    png(paste0("results/SNPs_", INV, "/LFMM_", i, "/", i, "_LatentFactors.png"), width = 500, height = 500)
    plot(mod.lfmm2@U, col = "blue", pch = 19, cex = 1.2)
    text(mod.lfmm2@U, DATA.af.meta$sampleId,
        cex = 0.65, pos = 3, col = "red"
    )
    dev.off()

    ### calculate P-values
    pv <- lfmm2.test(
        object = mod.lfmm2,
        input = DATA.af.lfmm,
        env = DATA.env, ,
        linear = TRUE
    )

    ### split into Inversion and all other SNPs
    PVAL <- cbind(ID = colnames(DATA.af.lfmm), t(-log10(pv$pvalues)))
    INVERSION <- PVAL[1, ]

    ### make input to plot Inversion pval separately
    hline_data <- data.frame(
        ylat = as.numeric(INVERSION[2]),
        ylon = as.numeric(INVERSION[3]),
        yat = as.numeric(INVERSION[4]),
        yap = as.numeric(INVERSION[5]),
        xmin = as.numeric(Start) / 1000000, # specify x-coordinate range for the lines
        xmax = as.numeric(End) / 1000000,
        Chrom = c(Chr)
    )

    ### prepare all SNPs for Manhattan plot
    OTHERS <- PVAL[2:nrow(PVAL), ]
    OTHERS <- as.data.frame(OTHERS) %>%
        separate(ID, into = c("Chrom", "Pos"), sep = ":")

    #### Latitude
    PLOT <- ggplot(
        OTHERS,
        aes(
            x = as.numeric(Pos) / 1000000,
            y = as.numeric(lat)
        )
    ) +
        geom_point(col = rgb(0, 0, 0, 0.1), pch = 16) +
        facet_grid(. ~ Chrom,
            scales = "free_x",
            space = "free"
        ) +
        theme_bw() +
        ## add Bonferroni-corrected p-value threshold
        geom_hline(yintercept = -log10(0.05 / (nrow(OTHERS) + 1)), colour = "blue") +
        ## add p-value of inversion in corresponding region
        geom_segment(aes(x = xmin, y = ylat, xend = xmax, yend = ylat, colour = "segment"),
            data = hline_data
        ) +
        xlab("Position (Mbp)") +
        ylab("-log10(p-value)") +
        ggtitle(paste0("LFMM: Latitude ", INV, " for ", i)) +
        theme(legend.position = "none")

    PlotList.lat[[i]] <- PLOT

    #### Longitude
    PLOT <- ggplot(
        OTHERS,
        aes(
            x = as.numeric(Pos) / 1000000,
            y = as.numeric(long)
        )
    ) +
        geom_point(col = rgb(0, 0, 0, 0.1), pch = 16) +
        facet_grid(. ~ Chrom,
            scales = "free_x",
            space = "free"
        ) +
        theme_bw() +
        ## add Bonferroni-corrected p-value threshold
        geom_hline(yintercept = -log10(0.05 / (nrow(OTHERS) + 1)), colour = "blue") +
        ## add p-value of inversion in corresponding region
        geom_segment(aes(x = xmin, y = ylon, xend = xmax, yend = ylon, colour = "segment"), data = hline_data) +
        xlab("Position (Mbp)") +
        ylab("-log10(p-value)") +
        ggtitle(paste0("LFMM: Longitude ", INV, " for ", i)) +
        theme(legend.position = "none")
    PlotList.lon[[i]] <- PLOT

    #### Average Temp
    PLOT <- ggplot(
        OTHERS,
        aes(
            x = as.numeric(Pos) / 1000000,
            y = as.numeric(AvTemp)
        )
    ) +
        geom_point(col = rgb(0, 0, 0, 0.1), pch = 16) +
        facet_grid(. ~ Chrom,
            scales = "free_x",
            space = "free"
        ) +
        theme_bw() +
        ## add Bonferroni-corrected p-value threshold
        geom_hline(yintercept = -log10(0.05 / (nrow(OTHERS) + 1)), colour = "blue") +
        ## add p-value of inversion in corresponding region
        geom_segment(aes(x = xmin, y = yat, xend = xmax, yend = yat, colour = "segment"),
            data = hline_data
        ) +
        xlab("Position (Mbp)") +
        ylab("-log10(p-value)") +
        ggtitle(paste0("LFMM: Average Temperature ", INV, " for ", i)) +
        theme(legend.position = "none")

    PlotList.at[[i]] <- PLOT

    #### Average Temp
    PLOT <- ggplot(
        OTHERS,
        aes(
            x = as.numeric(Pos) / 1000000,
            y = as.numeric(AvPrec)
        )
    ) +
        geom_point(col = rgb(0, 0, 0, 0.1), pch = 16) +
        facet_grid(. ~ Chrom,
            scales = "free_x",
            space = "free"
        ) +
        theme_bw() +
        ## add Bonferroni-corrected p-value threshold
        geom_hline(yintercept = -log10(0.05 / (nrow(OTHERS) + 1)), colour = "blue") +
        ## add p-value of inversion in corresponding region
        geom_segment(aes(x = xmin, y = yap, xend = xmax, yend = yap, colour = "segment"),
            data = hline_data
        ) +
        xlab("Position (Mbp)") +
        ylab("-log10(p-value)") +
        ggtitle(paste0("LFMM: Average Precipitation ", INV, " for ", i)) +
        theme(legend.position = "none")

    PlotList.ap[[i]] <- PLOT
}

## combine Plots of both continents per predictor (Lat/Long)
PLOT.lat <- ggarrange(plotlist = PlotList.lat, nrow = 2)
PLOT.lon <- ggarrange(plotlist = PlotList.lon, nrow = 2)
PLOT.ap <- ggarrange(plotlist = PlotList.ap, nrow = 2)
PLOT.at <- ggarrange(plotlist = PlotList.at, nrow = 2)

## save combined Plots
FILE <- paste0("results/SNPs_", INV, "/LFMM_", INV, "_Latitude.png")
ggsave(FILE,
    PLOT.lat,
    width = 10,
    height = 8
)

FILE <- paste0("results/SNPs_", INV, "/LFMM_", INV, "_Longitude.png")
ggsave(FILE,
    PLOT.lon,
    width = 10,
    height = 8
)

FILE <- paste0("results/SNPs_", INV, "/LFMM_", INV, "_AvTemp.png")
ggsave(FILE,
    PLOT.at,
    width = 10,
    height = 8
)

FILE <- paste0("results/SNPs_", INV, "/LFMM_", INV, "_Avprec.png")
ggsave(FILE,
    PLOT.ap,
    width = 10,
    height = 8
)
