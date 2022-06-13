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

#replace duplicate count NA with duplicate count = 0.5 (so that all TCR's appear on log plot)
subset_pat_439679_alpha_prod_wide[is.na(subset_pat_439679_alpha_prod_wide)] <- 0.5

#find total number of TCR's on each day
TCRs_day_0 <- sum(as.numeric(subset_pat_439679_alpha_prod_wide$duplicate_count.0))
TCRs_day_7 <- sum(as.numeric(subset_pat_439679_alpha_prod_wide$duplicate_count.7))
TCRs_day_14 <- sum(as.numeric(subset_pat_439679_alpha_prod_wide$duplicate_count.14))

#######day 0 to day 7 significance

p_value <- 0.0001

###find the lower and upper bounds for 95% of the data for poisson centered on 1st time point

############# This code is summing poisson distributions of means (lam_vals) from 0 up to (stops)
sum_poisson = function(lam, start=0, stop=5) {
  result = 0
  for (i in start:stop) {
    result = result + ((lam**(i))*(exp(-lam))/(factorial(i)))
  }
  return(result)
}

#the first timepoint data set
lam_vals <- subset_pat_439679_alpha_prod_wide$duplicate_count.0
#lam_vals <- c(1,2,51) #dummy array of first timepoint values

#the timepoint 7 dataset
end_time_vals <- subset_pat_439679_alpha_prod_wide$duplicate_count.7
#end_time_vals <- c(1,2,52) #dummy array



lam_vals <- subset_pat_439679_alpha_prod_wide$duplicate_count.0[0:10]

####dealing with the errors in the significance colour for having no TCR abundance on day 7 or day 0
#first of all, use the p value above to calculate what the threshold in this case should be
ar_before_1 <- vector()
count <- 0
for (z in 1:max(lam_vals)) {
  ar_before_1 <- sapply(z, function(x) sum_poisson(x,start=0,stop=1))
  #print(paste("ar_before_1",ar_before_1))
  if (ar_before_1 < (p_value/2)) {
    count=count+1 
    #print(paste("count",count))
  }
}
count_1_thresh <- max(lam_vals) - count
print(paste("count_1_thresh",count_1_thresh,"max(lam_vals)",max(lam_vals),"count",count))

#the binary array to see if significant or not
sig_day_7_from_day_0 <- vector()

for (p in 1:length(lam_vals)) {
  greater_than_min_sig <- vector() #initialise empty array 
  low_counter <- 0 #initialize lower significance level counter
  lesser_than_max_sig <- vector() #initialise empty array 
  high_counter <- 0 #initialize higher significance level counter
  # print(lam_vals[p])
  #logic to find upper significance level
  if (lam_vals[p] < 50) {     #lambda is low so manually sum to find significance level
    stops = linspace(0, 3*lam_vals[p], n = (3*lam_vals[p])-(0)+1)
    summation = sapply(lam_vals[p], function(x) { sapply(stops, function(y) sum_poisson(x, start=0, stop=y))})
    for (k in 1:length(summation)) {  #for loop to check each value in the array of summation
      if (summation[k] > (p_value/2)) {   #logic for finding the lower significance level
        low_counter = low_counter + 1
        greater_than_min_sig[low_counter] = k
      }
      if (summation[k] < (1-(p_value/2))) {  #logic for finding the upper significance level
        high_counter = high_counter + 1
        lesser_than_max_sig[high_counter] = k
      }
    }
    low_sig_lim = min(greater_than_min_sig)
    high_sig_lim = max(lesser_than_max_sig)
    
  } else {  #lambda is big so approximate as normal dist. p=0.05, which is 95% sig, which is 2 std_devs = 2*sqrt(mean) 
    low_sig_lim <- ceiling(lam_vals[p] - (2*sqrt(lam_vals[p])))
    high_sig_lim <- floor(lam_vals[p] + (2*sqrt(lam_vals[p])))
  }
  # print(paste("low_sig_lim",low_sig_lim))
  # print(paste("high_sig_lim",high_sig_lim))
  
  #compare the timepoint 7 data set to the significance thresholds
  if (end_time_vals[p] < low_sig_lim) { #it is significant
    sig_day_7_from_day_0[p] = 1
  } else if (end_time_vals[p] > high_sig_lim) {
    sig_day_7_from_day_0[p] = 1
  } else {
    sig_day_7_from_day_0[p] = 0
  }
  
  #correct the significance if duplicate count is less than 1 to check with threshold
  if (lam_vals[p] < 1) { #these are just the zero readings in the first count
    if (end_time_vals[p] < count_1_thresh) {
      sig_day_7_from_day_0[p] = 0}
  }
  if (end_time_vals[p] < 1) { #these are just the zero readings in the first count
    if (lam_vals[p] < count_1_thresh) {
      sig_day_7_from_day_0[p] = 0
    }
  }
}

# print(paste("sig_day_7_from_day_0",sig_day_7_from_day_0))
# print(paste("first_timepoint_vals",lam_vals))
# print(paste("end_time_vals",end_time_vals))

#add significance column to data frame
subset_pat_439679_alpha_prod_wide$sig_day_7_from_day_0 <- sig_day_7_from_day_0

#add column for day 0 to day 7 significant?       #THIS HAS BEEN REPLACED BY HUGE CHUNK ABOVE
# subset_pat_439679_alpha_prod_wide <- subset_pat_439679_alpha_prod_wide %>% 
#   mutate(sig_day_7_from_day_0 = if_else(duplicate_count.7 < duplicate_count.0 - sqrt(duplicate_count.0) | duplicate_count.7 > duplicate_count.0 + sqrt(duplicate_count.0) , 1, 0))

#deal with the cases of 0's and 1's where we want not significant
subset_pat_439679_alpha_prod_wide <- subset_pat_439679_alpha_prod_wide %>% 
  mutate(sig_day_7_from_day_0_deal_with_zeros = if_else(duplicate_count.0 == 0 & duplicate_count.7 == 1, 0, sig_day_7_from_day_0))

#make graph of TCR change over time
#asp = 1 makes axes equal. logic in xlim and ylim to choose the maximum of the x and y data. col is conditional colouring.
#plot(log(subset_pat_439679_alpha_prod_wide$duplicate_count.0/sum(subset_pat_439679_alpha_prod_wide$duplicate_count.0)), log(subset_pat_439679_alpha_prod_wide$duplicate_count.7/sum(subset_pat_439679_alpha_prod_wide$duplicate_count.7)), main = "log TCR amount day 0 to day 7",
plot(log(subset_pat_439679_alpha_prod_wide$duplicate_count.0), log(subset_pat_439679_alpha_prod_wide$duplicate_count.7), main = "log TCR amount day 0 to day 7",
     xlab = "log(TCR per million day 0/sum(TCR per million day 0))", ylab = "log(TCR per million day 7/sum(TCR per million day 7))",asp = 1,
     xlim = c(log(min(subset_pat_439679_alpha_prod_wide$duplicate_count.0)), if_else(max(subset_pat_439679_alpha_prod_wide$duplicate_count.0) > max(subset_pat_439679_alpha_prod_wide$duplicate_count.7), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.0)), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.7)))),
     #xlim = c(log(min(subset_pat_439679_alpha_prod_wide$duplicate_count.0/sum(subset_pat_439679_alpha_prod_wide$duplicate_count.0))), if_else(max(subset_pat_439679_alpha_prod_wide$duplicate_count.0/sum(subset_pat_439679_alpha_prod_wide$duplicate_count.0)) > max(subset_pat_439679_alpha_prod_wide$duplicate_count.7/sum(subset_pat_439679_alpha_prod_wide$duplicate_count.7)), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.0/sum(subset_pat_439679_alpha_prod_wide$duplicate_count.0))), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.7/sum(subset_pat_439679_alpha_prod_wide$duplicate_count.7))))),
     #ylim = c(log(min(subset_pat_439679_alpha_prod_wide$duplicate_count.7/sum(subset_pat_439679_alpha_prod_wide$duplicate_count.7))), if_else(max(subset_pat_439679_alpha_prod_wide$duplicate_count.0/sum(subset_pat_439679_alpha_prod_wide$duplicate_count.0)) > max(subset_pat_439679_alpha_prod_wide$duplicate_count.7/sum(subset_pat_439679_alpha_prod_wide$duplicate_count.7)), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.0/sum(subset_pat_439679_alpha_prod_wide$duplicate_count.0))), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.7/sum(subset_pat_439679_alpha_prod_wide$duplicate_count.7))))), 
     ylim = c(log(min(subset_pat_439679_alpha_prod_wide$duplicate_count.7)), if_else(max(subset_pat_439679_alpha_prod_wide$duplicate_count.0) > max(subset_pat_439679_alpha_prod_wide$duplicate_count.7), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.0)), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.7)))), 
     col = ifelse(subset_pat_439679_alpha_prod_wide$sig_day_7_from_day_0_deal_with_zeros == 1,'red','blue'))
abline(a=0,b=1)



