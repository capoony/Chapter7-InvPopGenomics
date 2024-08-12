# Load necessary libraries
library(readr)
library(tidyverse)
library(lme4)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
INV <- args[1]
Chr <- args[2]
Start <- args[3]
End <- args[4]
WD <- args[5]

# Set working directory
setwd(WD)

# Function to load and preprocess data
load_and_preprocess_data <- function(data_matrix_path, weights_matrix_path, dependent_variable_path) {
    # Load the data matrix, weights matrix, and dependent variable DataFrame
    data_matrix <- read.table(gzfile(data_matrix_path), header = TRUE, comment.char = "")
    weights_matrix <- read.table(gzfile(weights_matrix_path), header = TRUE, comment.char = "")
    dependent_variable <- read.table(dependent_variable_path, header = TRUE)
    colnames(dependent_variable) <- c("Label", "Value")
    dependent_variable$Label <- gsub("-", ".", dependent_variable$Label)

    # Extract the labels from the data matrix header, excluding the first two columns
    labels <- colnames(data_matrix)[-c(1, 2)]

    # Subset the dependent variable DataFrame to match the labels in the data matrix
    dependent_variable_subset <- dependent_variable %>%
        filter(Label %in% labels) %>%
        arrange(factor(Label, levels = labels))

    # Subset the data matrix and weights matrix to match the labels in the dependent variable
    data_matrix_subset <- data_matrix %>%
        select(c(1, 2, all_of(dependent_variable_subset$Label)))

    weights_matrix_subset <- weights_matrix %>%
        select(c(1, 2, all_of(dependent_variable_subset$Label)))

    # Create the independent variable vector in the correct order
    X <- dependent_variable_subset$Value

    # Return the processed data
    list(
        data_matrix_subset = data_matrix_subset,
        weights_matrix_subset = weights_matrix_subset,
        X = X
    )
}

# Function to perform logistic regression for each row and compute p-values
logistic_regression_function <- function(row_data, row_weights, X) {
    # Create a data frame for the model
    model_data <- data.frame(
        y = row_data,
        X = X
    )

    # Fit the logistic regression model using glm, incorporating weights
    model <- glm(y ~ X, data = model_data, family = binomial, weights = row_weights)

    # Compute p-value using anova
    anova_results <- anova(model, test = "Chisq")
    p_value <- anova_results$`Pr(>Chi)`[2]

    return(p_value)
}

# Main function to perform the analysis
perform_analysis <- function(data_matrix_path, weights_matrix_path, dependent_variable_path) {
    # Load and preprocess the data
    processed_data <- load_and_preprocess_data(data_matrix_path, weights_matrix_path, dependent_variable_path)

    data_matrix_subset <- processed_data$data_matrix_subset
    weights_matrix_subset <- processed_data$weights_matrix_subset
    X <- processed_data$X

    # Apply the logistic regression function to each row in the subsetted data matrix
    p_values <- sapply(1:nrow(data_matrix_subset), function(row_index) {
        row_data <- as.numeric(data_matrix_subset[row_index, -c(1, 2)])
        row_weights <- as.numeric(weights_matrix_subset[row_index, -c(1, 2)])
        logistic_regression_function(row_data, row_weights, X)
    })

    # Create a dataframe with p-values and label columns from the data matrix
    p_values_df <- data.frame(
        Chrom = data_matrix_subset[, 1],
        Pos = data_matrix_subset[, 2],
        p_value = p_values
    )

    return(p_values_df)
}

# Create directory for results
dir.create(paste0("results/SNPs_", INV, "/LDwithSNPs"))

# Perform analysis for Europe and North America
p_values_df.Eur <- perform_analysis(
    "results/SNPs/Europe_freq.csv.gz",
    "results/SNPs/Europe_weight.csv.gz",
    paste0("results/SNPs_", INV, "/", INV, ".af")
)
p_values_df.NA <- perform_analysis(
    "results/SNPs/NorthAmerica_freq.csv.gz",
    "results/SNPs/NorthAmerica_weight.csv.gz",
    paste0("results/SNPs_", INV, "/", INV, ".af")
)

# Add continent labels
p_values_df.Eur$Continent <- "Europe"
p_values_df.NA$Continent <- "North America"

# Combine results
p_values_df <- rbind(p_values_df.Eur, p_values_df.NA)

# Save p-values to a file
write.table(p_values_df,
    file = paste0("results/SNPs_", INV, "/LDwithSNPs/SNPs.pval"),
    quote = FALSE,
    row.names = FALSE
)

# Create data frame for the rectangle (highlight region) in the plot
rect_data <- data.frame(
    xmin = as.numeric(Start) / 1000000,
    xmax = as.numeric(End) / 1000000,
    ymin = 0,
    ymax = Inf,
    Chrom = Chr
)

# Create Manhattan plot
PLOT <- ggplot(p_values_df, aes(x = Pos / 1000000, y = -log10(p_value))) +
    geom_point(color = rgb(0, 0, 0, 0.5), pch = 16) +
    facet_grid(Continent ~ Chrom, scales = "free", space = "free_x") +
    geom_rect(
        data = rect_data,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "blue", alpha = 0.2, inherit.aes = FALSE
    ) +
    theme_bw() +
    xlab("Position (Mbp)") +
    ggtitle(paste0("Logistic regression with ", INV))

# Save the plot
FILE <- paste0("results/SNPs_", INV, "/LDwithSNPs/", INV, "_LD.png")
ggsave(
    file = FILE,
    plot = PLOT,
    width = 10,
    height = 6
)

# Save the plot
FILE <- paste0("results/SNPs_", INV, "/LDwithSNPs/", INV, "_LD.pdf")
ggsave(
    file = FILE,
    plot = PLOT,
    width = 10,
    height = 6
)
