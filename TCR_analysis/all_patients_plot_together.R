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

all_pats_wide[is.na(all_pats_wide)] <- 0  #replace duplicate count = NA with duplicate count = 0

rownames(all_pats_wide) <- NULL  #This resets the index of the rows of the dataframe

#######the significance level

p_value <- 0.0001

############# This function is summing poisson distributions of means (lam_vals) from 0 up to (stops)
sum_poisson = function(lam, start=0, stop=5) {
  result = 0
  for (i in start:stop) {
    result = result + ((lam**(i))*(exp(-lam))/(factorial(i)))
  }
  return(result)
}

#the first timepoint data set
#lam_vals <- all_pats_wide$duplicate_count.0
lam_vals <- all_pats_wide$duplicate_count.0[1:500] #just taking a subset
#lam_vals <- c(30,0,2,0,2,30,1)

#the end timepoint dataset
#end_time_vals <- all_pats_wide$duplicate_count.7
end_time_vals <- all_pats_wide$duplicate_count.7[1:500] #just taking a subset
#end_time_vals <- c(10,0,2,7,3,3,0)

####Cannot have a poisson dist centred at 0. Therefore take significance boundary for poisson centred at 1 and take same limit. Make sig line straight at LHS of graph.
ar_bef_1 <- vector()
count2 <- 0
for (z in 1:max(end_time_vals)) {  #from 1 up to the maximum of the day 7 count values
  ar_bef_1 <- sapply(1, function(x) sum_poisson(x,start=0,stop=z))
  #print(paste("ar_bef_1",ar_bef_1))
  if (ar_bef_1 > (1-(p_value/2))) {
    count2 = count2+1
  }
}
count_2_thresh <- max(end_time_vals) - count2




