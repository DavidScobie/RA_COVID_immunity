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
library(pracma)
library(MASS)

#read in the 3 timepoints for 1st patient alpha chain
data_path<-"C:/Research_Assistant/work/data/TCR_data/"
input_data_chop_all<-read.csv(paste0(data_path,"input_data_chop_all_reg_days_no_filename.csv"))

# #exclude all rows which contain day10, day13 and post as these make data frame too big
# input_data_chop_all_without_10 <- input_data_chop_all[!grepl("day10", input_data_chop_all$filename),]
# input_data_chop_all_without_10_13 <- input_data_chop_all_without_10[!grepl("day13", input_data_chop_all_without_10$filename),]
# input_data_chop_all_without_reg_days <- input_data_chop_all_without_10_13[!grepl("post", input_data_chop_all_without_10_13$filename),]
# 
# #write the input_data_chop to a new file
# data_path <- "C:/Research_Assistant/work/data/TCR_data/"
# write.csv(input_data_chop_all_reg_days_no_filename,paste0(data_path,'input_data_chop_all_reg_days_no_filename.csv'), row.names = FALSE)
# 
# myvars <- c("junction_aa", "duplicate_count", "day", "subject", "chain", "PCR_positive")
# input_data_chop_all_reg_days_no_filename <- input_data_chop_all_without_reg_days[myvars]
# 
# input_data_chop_all_wide <- reshape(input_data_chop_all, idvar=c('junction_aa','subject','chain','PCR_positive'),timevar='day',direction='wide')

#chop out all the data from day 14 from the data set, as this will give us too many counts which are (0,0)
input_data_chop_all_without_14 <- input_data_chop_all[!grepl(14, input_data_chop_all$day),]

#make a list containing all of the possible subjects for subsetting the data
list_of_subjects <- as.list(unique(input_data_chop_all_without_14$subject))

#filter out individual patients
all_pats_wide <- data.frame()
#for (n in 1:length(list_of_subjects)) {
for (n in 1:2) {
  a_subject <- input_data_chop_all_without_14 %>% filter(subject == list_of_subjects[n])
  print(paste("patient",list_of_subjects[n]))
  #next make the data frame wide
  a_subject_wide <- reshape(a_subject, idvar=c('junction_aa','subject','chain','PCR_positive'),timevar='day',direction='wide')
  print(paste("number of rows of a_subject_wide",nrow(a_subject_wide)))
  #bind this dataframe onto the dataframes for the other subjects
  all_pats_wide <- bind_rows(all_pats_wide,a_subject_wide)
  print(paste("number of rows of all_pats_wide",nrow(all_pats_wide)))
}






