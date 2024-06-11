library(tidyverse)
library(factoextra)
library(FactoMineR)
library(ggpubr)

# Define command line arguments
args <- commandArgs(trailingOnly = TRUE)
INV <- args[1]
Chr <- args[2]
Start <- as.numeric(args[3])
End <- as.numeric(args[4])
WD <- args[5]

# Set working directory
setwd(WD)

map_levels_to_numbers <- function(vec, start, end) {
    # Get the unique levels in the vector
    unique_levels <- unique(vec)

    # Create a mapping of unique levels to numbers ranging from 1 to start, repeated
    level_mapping <- setNames(rep(start:end, length.out = length(unique_levels)), unique_levels)

    # Replace the character levels with corresponding numbers
    num_vec <- level_mapping[vec]

    # Return the resulting vector
    return(num_vec)
}

# Function to load and process metadata
load_metadata <- function() {
    meta <- read.csv("data/meta.csv", header = TRUE)
    ## subset to only contain the sampleId, country and province columns
    meta.sub <- meta %>% select(sampleId, country, province)
    ## modify the sampleId to be consistent with allele frequency datase
    meta.sub$sampleId <- gsub("-", ".", meta.sub$sampleId)
    return(meta.sub)
}

# Function to read frequency data
read_frequency_data <- function(region) {
    ## read the allele frequency dataset
    DATA <- read.table(gzfile(paste0("results/SNPs/", region, "_freq.csv.gz")), header = TRUE, comment.char = "")
    return(DATA)
}

# Function to process inversion data
process_inversion_data <- function(DATA, meta.sub, Chr, Start, End) {
    ## transpose the allele frequeny matrix and only retain SNPs, within the inverted genomic region.
    DATA.inv <- as.data.frame(t(DATA %>%
        filter(X.CHROM == Chr & POS > Start & POS < End) %>%
        select(3:ncol(DATA))))
    ## make new column from rownames (sampleIds)
    DATA.inv$sampleId <- rownames(DATA.inv)
    ## merge the allele frequency data and the country and province information from the metadata set
    DATA.inv <- DATA.inv %>%
        inner_join(meta.sub, by = "sampleId") %>%
        ## get rid of Panama and Guadeloupe since they are Central America and not North America
        filter(country != "Panama" & country != "Guadeloupe")
    return(DATA.inv)
}

# Function to perform PCA
perform_pca <- function(DATA) {
    PCA.result <- PCA(DATA, graph = FALSE)
    return(PCA.result)
}

# Function to create PCA plot
create_pca_plot <- function(PCA.result, DATA, region, INV, inside = TRUE) {
    plot_type <- ifelse(inside, "inside", "outside")
    if (region == "Europe") {
        COLOR <- DATA$country
    } else {
        COLOR <- DATA$province
    }
    COLOR2 <- map_levels_to_numbers(COLOR, 1, 8)
    SHAPE <- map_levels_to_numbers(COLOR, 15, 19)
    PLOT <- fviz_pca_ind(PCA.result,
        col.ind = COLOR,
        pointsize = 3,
        alpha = 0.7,
        mean.point = FALSE,
        label = COLOR,
        repel = TRUE
    ) +
        theme_bw() +
        ggtitle(paste0("PCA - ", region, " ", plot_type, " ", INV)) +
        labs(color = "Count(r)y") +
        labs(shape = "Count(r)y") +
        guides(fill = guide_legend(nrow = 2)) +
        scale_shape_manual(values = SHAPE) +
        scale_color_manual(values = COLOR2)
    return(PLOT)
}

# Function to save PCA results
save_pca_results <- function(PCA.result, DATA, region, INV, inside = TRUE) {
    plot_type <- ifelse(inside, "inside", "outside")
    write.table(
        file = paste0("results/SNPs_", INV, "/PCA_", INV, "_", region, "_", plot_type, ".txt"),
        cbind(DATA[, (ncol(DATA) - 2):ncol(DATA)], PCA.result$ind$coord),
        sep = ",",
        quote = FALSE,
        row.names = FALSE
    )
}

# Main function to process data for each region
process_region <- function(region, meta.sub, Chr, Start, End, INV) {
    DATA <- read_frequency_data(region)

    # Process SNP data inside the inversion
    DATA.inv <- process_inversion_data(DATA, meta.sub, Chr, Start, End)
    PCA.inv <- perform_pca(DATA.inv[, 1:(ncol(DATA.inv) - 3)])
    save_pca_results(PCA.inv, DATA.inv, region, INV, inside = TRUE)
    PLOT.inv <- create_pca_plot(PCA.inv, DATA.inv, region, INV, inside = TRUE)

    # Process SNP data outside the inversion
    DATA.non <- as.data.frame(t(DATA %>%
        filter(!(X.CHROM == Chr & POS > Start & POS < End)) %>%
        select(3:ncol(DATA))))
    DATA.non$sampleId <- rownames(DATA.non)
    DATA.non <- DATA.non %>%
        inner_join(meta.sub, by = "sampleId") %>%
        filter(country != "Panama" & country != "Guadeloupe")
    PCA.non <- perform_pca(DATA.non[, 1:(ncol(DATA.non) - 3)])
    save_pca_results(PCA.non, DATA.non, region, INV, inside = FALSE)
    PLOT.non <- create_pca_plot(PCA.non, DATA.non, region, INV, inside = FALSE)

    # Combine and save plots
    PLOT <- ggarrange(PLOT.inv, PLOT.non,
        common.legend = TRUE,
        legend = "bottom"
    )
    FILE <- paste0("results/SNPs_", INV, "/PCA_", INV, "_", region, ".png")
    ggsave(file = FILE, PLOT, width = 10, height = 5)
}

# Load metadata
meta.sub <- load_metadata()

# Process each region
for (region in c("Europe", "NorthAmerica")) {
    process_region(region, meta.sub, Chr, Start, End, INV)
}
