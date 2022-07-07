options(stringsAsFactors=FALSE)
require(gdata)
require(dplyr)
require(ggplot2)
#library(XLConnect)
library(xlsx)
library(data.table)
library(stringr)
library(tidyr)
library(pracma)
library(MASS)

#read the input data in to begin the analysis
data_path <- "C:/Research_Assistant/work/data/TCR_data/NS087/"
input_data<-read.csv(paste0(data_path,"result.csv"))

#only keep columns with productive=TRUE and stop_codon=FALSE
input_data_prod <- subset(input_data, subset = productive == TRUE)
input_data_prod_cod <- subset(input_data_prod, subset = stop_codon == FALSE)

# keep only the columns that are necessary in the data frame
myvars <- c("filename", "junction_aa", "duplicate_count")
input_data_chop <- input_data_prod_cod[myvars]

# Adding column of day based on string in filename:
input_data_chop %>%
  mutate(day = case_when(
    grepl("pre_1", filename, fixed = TRUE) ~ 0,
    grepl("post_1", filename, fixed = TRUE) ~ 0,
    grepl("day14", filename, fixed = TRUE) ~ 14,
    grepl("day7", filename, fixed = TRUE) ~ 7,
    grepl("day13", filename, fixed = TRUE) ~ 13,
    grepl("day10", filename, fixed = TRUE) ~ 10,
  ))

# Adding column of subject based on string in filename:
input_data_chop %>%
  mutate(day = case_when(
    grepl("550630", filename, fixed = TRUE) ~ 550630,
    grepl("627506", filename, fixed = TRUE) ~ 627506,
    grepl("634105", filename, fixed = TRUE) ~ 634105,
    grepl("635331", filename, fixed = TRUE) ~ 635331,
    grepl("635729", filename, fixed = TRUE) ~ 635729,
    grepl("635779", filename, fixed = TRUE) ~ 635779,
    grepl("637340", filename, fixed = TRUE) ~ 637340,
    grepl("647785", filename, fixed = TRUE) ~ 647785,
    grepl("635495", filename, fixed = TRUE) ~ 635495,
    grepl("651806", filename, fixed = TRUE) ~ 651806,
    grepl("666427", filename, fixed = TRUE) ~ 666427,
    grepl("666482", filename, fixed = TRUE) ~ 666482,
    grepl("666660", filename, fixed = TRUE) ~ 666660,
    grepl("667186", filename, fixed = TRUE) ~ 667186,
    grepl("439679", filename, fixed = TRUE) ~ 439679,
    grepl("635742", filename, fixed = TRUE) ~ 635742,
    grepl("636163", filename, fixed = TRUE) ~ 636163,
    grepl("637269", filename, fixed = TRUE) ~ 637269,
    grepl("643925", filename, fixed = TRUE) ~ 643925,
    grepl("634418", filename, fixed = TRUE) ~ 634418,
    grepl("645438", filename, fixed = TRUE) ~ 645438,
    grepl("655401", filename, fixed = TRUE) ~ 655401,
    grepl("666475", filename, fixed = TRUE) ~ 666475,
    grepl("667082", filename, fixed = TRUE) ~ 667082,
    grepl("672533", filename, fixed = TRUE) ~ 672533,
    grepl("673067", filename, fixed = TRUE) ~ 673067,
  ))

