# R version 3.5.0
# gdata version 2.18.0
# dplyr version 0.8.3
# ggplot2 version 3.2.1

options(stringsAsFactors=FALSE)
require(gdata)
require(dplyr)
require(ggplot2)
#library(XLConnect)
library(xlsx)
library(data.table)
library(stringr)
library(tidyr)

#read in the 3 timepoints for 1st patient alpha chain
data_path<-"C:/Research_Assistant/work/data/TCR_data/NS085/"
pat_439679_day_0_alpha<-read.csv(paste0(data_path,"dcr_HVO_439679_pre_1_alpha.csv"))
pat_439679_day_7_alpha<-read.csv(paste0(data_path,"dcr_HVO_439679_day7_1_alpha.csv"))
pat_439679_day_14_alpha<-read.csv(paste0(data_path,"dcr_HVO_439679_day14_1_alpha.csv"))

#only keep columns with productive=TRUE
pat_439679_day_0_alpha_prod <- subset(pat_439679_day_0_alpha, subset = productive == TRUE) 
pat_439679_day_7_alpha_prod <- subset(pat_439679_day_7_alpha, subset = productive == TRUE)
pat_439679_day_14_alpha_prod <- subset(pat_439679_day_14_alpha, subset = productive == TRUE)

# add column to dataframe of the day
pat_439679_day_0_alpha_prod$day <- 0
pat_439679_day_7_alpha_prod$day <- 7
pat_439679_day_14_alpha_prod$day <- 14

#combine the 3 data frames for patient 439679 into 1 big data frame (to get the long format)
pat_439679_alpha_prod <- bind_rows(pat_439679_day_0_alpha_prod, pat_439679_day_7_alpha_prod, pat_439679_day_14_alpha_prod)

#make the data frame wide
pat_439679_alpha_prod_wide <- reshape(pat_439679_alpha_prod, idvar = "junction_aa", timevar = "day", direction = "wide")

#count how many ones there are in duplicate_count
day_0_over_1 <- nrow(pat_439679_alpha_prod_wide[pat_439679_alpha_prod_wide$duplicate_count.0>1, ])
day_0_over_0 <- nrow(pat_439679_alpha_prod_wide[pat_439679_alpha_prod_wide$duplicate_count.0>0, ])

# select variables junction_aa, duplicate_count.0, duplicate_count.7, duplicate_count.14
myvars <- c("junction_aa", "duplicate_count.0", "duplicate_count.7", "duplicate_count.14")
subset_pat_439679_alpha_prod_wide <- pat_439679_alpha_prod_wide[myvars]

#replace duplicate count NA with duplicate count = 0
subset_pat_439679_alpha_prod_wide[is.na(subset_pat_439679_alpha_prod_wide)] <- 0

#add column for day 7 significant?
subset_pat_439679_alpha_prod_wide <- subset_pat_439679_alpha_prod_wide %>% 
  mutate(sig_day_7_from_day_0 = if_else(duplicate_count.7 <= duplicate_count.0 - sqrt(duplicate_count.0) | duplicate_count.7 >= duplicate_count.0 + sqrt(duplicate_count.0) , 1, 0))
