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

#read in the 3 timepoints for 1st patient alpha chain
data_path<-"C:/Research_Assistant/work/data/TCR_data/NS085/"
pat_439679_day_0_alpha<-read.csv(paste0(data_path,"dcr_HVO_439679_pre_1_alpha.csv"))
pat_439679_day_7_alpha<-read.csv(paste0(data_path,"dcr_HVO_439679_day7_1_alpha.csv"))
pat_439679_day_14_alpha<-read.csv(paste0(data_path,"dcr_HVO_439679_day14_1_alpha.csv"))

#only keep columns with productive=TRUE
pat_439679_day_0_alpha_prod <- subset(pat_439679_day_0_alpha, subset = productive == TRUE) 
pat_439679_day_7_alpha_prod <- subset(pat_439679_day_7_alpha, subset = productive == TRUE)
pat_439679_day_14_alpha_prod <- subset(pat_439679_day_14_alpha, subset = productive == TRUE)

###make a list of all the TCRs at day 0 and check if they are in the data at day 7 and day 14
pat_439679_day_0_alpha_prod$present_in_day_7 <- as.integer(pat_439679_day_0_alpha_prod$junction_aa %in% pat_439679_day_7_alpha_prod$junction_aa)
pat_439679_day_0_alpha_prod$present_in_day_14 <- as.integer(pat_439679_day_0_alpha_prod$junction_aa %in% pat_439679_day_14_alpha_prod$junction_aa)

###if there is a zero in the present_in_day_7 column want to inflate the day 7 data frame with a row of that juntion_aa but a 0 as duplicate count
N_TCRs <- nrow(pat_439679_day_0_alpha_prod) #find the number of rows for the first patient
k<-1
for (k in 1:N_TCRs){
  if pat_439679_day_0_alpha_prod$present_in_day_7[k] == 0{
    #create a new row in the data frame with that TCR but abundance zero
