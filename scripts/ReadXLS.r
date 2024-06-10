# Load necessary libraries
library(readxl)
library(tidyverse)

# Parse command line arguments and set working directory
args <- commandArgs(trailingOnly = TRUE)
WD <- args[1]
setwd(WD)

# Read the Excel file, skipping the first 6 rows
myDF <- read_excel("data/TableS1_individuals.xls", skip = 6)

# Filter the dataframe to keep only haploid embryos with single SRA accessions that start with 'SRR'
myDF1 <- myDF %>%
    filter(`Genome Type` == "haploid_embryo" & !str_detect(`SRA Accession`, ",") & str_detect(`SRA Accession`, "SRR"))

# Get the first 20 In(2L)t inverted and standard individuals
IN2LT.yes <- myDF1 %>%
    filter(`In(2L)t` == "INV" & `In(3R)P` == "ST") %>%
    select(`Stock ID`, `SRA Accession`, `In(2L)t`) %>%
    head(20)

IN2LT.no <- myDF1 %>%
    filter(`In(2L)t` == "ST" & `In(3R)P` == "ST") %>%
    select(`Stock ID`, `SRA Accession`, `In(2L)t`) %>%
    head(20)

# Write the combined data to a file
write.table(rbind(IN2LT.yes, IN2LT.no),
    file = "data/IN2LT.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = ","
)

# Get the first 20 In(3R)Payne inverted and standard individuals
IN3RP.yes <- myDF1 %>%
    filter(`In(3R)P` == "INV" & `In(2L)t` == "ST") %>%
    select(`Stock ID`, `SRA Accession`, `In(3R)P`) %>%
    head(20)

IN3RP.no <- myDF1 %>%
    filter(`In(3R)P` == "ST" & `In(2L)t` == "ST") %>%
    select(`Stock ID`, `SRA Accession`, `In(3R)P`) %>%
    head(20)

# Write the combined data to a file
write.table(rbind(IN3RP.yes, IN3RP.no),
    file = "data/IN3RP.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = ","
)
