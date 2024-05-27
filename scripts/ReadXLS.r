## load libraries
library(readxl)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)

WD <- args[1]

setwd(WD)

## read XLS file with SRA IDs and Inversion status
myDF <- read_excel("TableS1_individuals.xls",
    skip = 6
)

## only keep certain columns
myDF1 <- myDF[myDF[["Genome Type"]] == "haploid_embryo" & !grepl(",", myDF[["SRA Accession"]]) & grepl("SRR", myDF[["SRA Accession"]]), ]

## get In(2L)t inverted and standard individuals
IN2LT.yes <- head(myDF1[myDF1[["In(2L)t"]] == "INV" & myDF1[["In(3R)P"]] == "ST", ]
%>% select(`Stock ID`, `SRA Accession`, `In(2L)t`), 20)
IN2LT.no <- head(myDF1[myDF1[["In(2L)t"]] == "ST" & myDF1[["In(3R)P"]] == "ST", ]
%>% select(`Stock ID`, `SRA Accession`, `In(2L)t`), 20)
write.table(rbind(IN2LT.yes, IN2LT.no),
    file = "IN2LT.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = ","
)

# ## get In(3R)Payne inverted and standard individuals
IN3RP.yes <- head(myDF1[myDF1[["In(3R)P"]] == "INV" & myDF1[["In(2L)t"]] == "ST", ]
%>% select(`Stock ID`, `SRA Accession`, `In(3R)P`), 20)
IN3RP.no <- head(myDF1[myDF1[["In(3R)P"]] == "ST" & myDF1[["In(2L)t"]] == "ST", ]
%>% select(`Stock ID`, `SRA Accession`, `In(3R)P`), 20)
write.table(rbind(IN3RP.yes, IN3RP.no),
    file = "IN3RP.txt",
    quote = FALSE,
    row.names = FALSE,
    sep = ","
)
