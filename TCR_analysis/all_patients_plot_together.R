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

##To start with we want to look at all PCR+ve patients from day 0 to day 7, on both alpha and beta chains

#read in the data for all of the patients with no filename, on days 0, 7 and 14, with only necessary data inside
data_path<-"C:/Research_Assistant/work/data/TCR_data/"
input_data_chop_all<-read.csv(paste0(data_path,"input_data_chop_all_reg_days_no_filename.csv"))

#chop out all the patients who are PCR-ve, as we just want PCR+ve people in this case
input_data_chop_all_PCR_pos <- input_data_chop_all[!grepl(0, input_data_chop_all$PCR_positive),]

#chop out all the data from day 14 from the data set, as this will give us too many counts which are (0,0)
input_data_chop_all_without_14 <- input_data_chop_all_PCR_pos[!grepl(14, input_data_chop_all_PCR_pos$day),]

#make a list containing all of the possible subjects for subsetting the data
list_of_subjects <- as.list(unique(input_data_chop_all_without_14$subject))

#filter out individual patients and make individual wide dataframes for them. Then bind all these dataframes together
all_pats_wide <- data.frame()
for (n in 1:length(list_of_subjects)) {  #from 1 to the number of patients
  a_subject <- input_data_chop_all_without_14 %>% filter(subject == list_of_subjects[n])  #subsetting the full dataframe for just 1 subject
  #make the wide dataframe. Keeping 'junction_aa','subject','chain','PCR_positive'. And making 'duplicate_count.x', where x are the days in question (0 and 7 to begin with)
  a_subject_wide <- reshape(a_subject, idvar=c('junction_aa','subject','chain','PCR_positive'),timevar='day',direction='wide')
  #bind the dataframe for this subject onto the dataframes for the other subjects
  all_pats_wide <- bind_rows(all_pats_wide,a_subject_wide)
}






