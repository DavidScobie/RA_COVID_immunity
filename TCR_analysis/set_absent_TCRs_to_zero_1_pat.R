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

#replace duplicate count NA with duplicate count = 0 
subset_pat_439679_alpha_prod_wide[is.na(subset_pat_439679_alpha_prod_wide)] <- 0

#find total number of TCR's on each day
TCRs_day_0 <- sum(as.numeric(subset_pat_439679_alpha_prod_wide$duplicate_count.0))
TCRs_day_7 <- sum(as.numeric(subset_pat_439679_alpha_prod_wide$duplicate_count.7))
TCRs_day_14 <- sum(as.numeric(subset_pat_439679_alpha_prod_wide$duplicate_count.14))

#######day 0 to day 7 significance

p_value <- 0.0001

###find the lower and upper bounds for p value of the data for poisson centered on 1st time point

############# This function is summing poisson distributions of means (lam_vals) from 0 up to (stops)
sum_poisson = function(lam, start=0, stop=5) {
  result = 0
  for (i in start:stop) {
    result = result + ((lam**(i))*(exp(-lam))/(factorial(i)))
  }
  return(result)
}

#the first timepoint data set
lam_vals <- subset_pat_439679_alpha_prod_wide$duplicate_count.0
#lam_vals <- c(0,0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2) #dummy array of first timepoint values

#the timepoint 7 dataset
end_time_vals <- subset_pat_439679_alpha_prod_wide$duplicate_count.7
#end_time_vals <- c(4,5,6,7,8,9,10,6,7,8,9,10,6,7,8,9,10) #dummy array


#lam_vals <- subset_pat_439679_alpha_prod_wide$duplicate_count.0[0:10]

####dealing with the errors in the significance colour for having no TCR abundance on day 7 
ar_before_1 <- vector()
count <- 0
for (z in 1:max(lam_vals)) {
  ar_before_1 <- sapply(z, function(x) sum_poisson(x,start=0,stop=0)) #here we are finding height of 1st line in poisson centered at means from 1 to maximum day 0 value
  if (ar_before_1 < (p_value/2)) { #use the p value above to calculate what the threshold in this case should be
    count=count+1 #count increases if height of 1st line is below the threshold set by (p value/2). (So it is significant)
  }
}
count_1_thresh <- max(lam_vals) - count  #this finds the lowest value of the mean where the sum at x=0 is significant

####dealing with the errors in the significance colour for having no TCR abundance on day 0
ar_bef_1 <- vector()
count2 <- 0
for (z in 1:max(end_time_vals)) {
  ar_bef_1 <- sapply(1, function(x) sum_poisson(x,start=0,stop=z))
  #print(paste("ar_bef_1",ar_bef_1))
  if (ar_bef_1 > (1-(p_value/2))) {  
    count2 = count2+1
  }
}
count_2_thresh <- max(end_time_vals) - count2
#######

#the binary array to see if significant or not
sig_day_7_from_day_0 <- vector()

start_day_0_count = 0
start_day_7_count = 0
manual_x_threshold =50 #what value of x do we start approximating sig with (mean +- (num_std_devs*sqrt(lambda))) ? 

for (p in 1:length(lam_vals)) {
  greater_than_min_sig <- vector() #initialise empty array 
  low_counter <- 0 #initialize lower significance level counter
  lesser_than_max_sig <- vector() #initialise empty array 
  high_counter <- 0 #initialize higher significance level counter
  #logic to find upper significance level
  if (lam_vals[p] < manual_x_threshold) {     #lambda is low so manually sum to find significance level
    
    #this loop is in to save memory. We only need large num_lam_multip for small lambdas. Can always add steps in the loop to be more memory efficient
    if (lam_vals[p] < 10) {
      num_lam_multip = 12  #number of multiples of lambda to sum up to in the poisson distribution
    } else {
      num_lam_multip = 3
    }
    
    stops = linspace(0, num_lam_multip*lam_vals[p], n = (num_lam_multip*lam_vals[p])-(0)+1)  #need stops to go far beyond
    summation = sapply(lam_vals[p], function(x) { sapply(stops, function(y) sum_poisson(x, start=0, stop=y))})
    #print(paste("summation",summation))
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
    low_sig_lim = min(greater_than_min_sig) -1 #we need to allow this to be 0 if necessary. Need -1 to work with indices as k in 1:length(summation)
    high_sig_lim = max(lesser_than_max_sig) -1 #we need to allow this to be 0 if necessary. Need -1 to work with indices as k in 1:length(summation)
    
  } else {  #lambda is big so approximate as normal dist. p=0.05, which is 95% sig, which is 2 std_devs = 2*sqrt(mean) 
    num_std_devs_gauss <- qnorm(1-(p_value/2))
    low_sig_lim <- ceiling(lam_vals[p] - (num_std_devs_gauss*sqrt(lam_vals[p])))
    high_sig_lim <- floor(lam_vals[p] + (num_std_devs_gauss*sqrt(lam_vals[p])))
  }
  #print(paste("lesser_than_max_sig",lesser_than_max_sig))
  #print(paste("high_sig_lim",high_sig_lim))
  
  #compare the timepoint 7 data set to the significance thresholds
  if (end_time_vals[p] < low_sig_lim) { #it is significant at lower end
    #print(paste("SIGNIFICANT"))
    sig_day_7_from_day_0[p] = 1
  } else if (end_time_vals[p] > high_sig_lim) { #it is significant at higher end
    #print(paste("SIGNIFICANT"))
    sig_day_7_from_day_0[p] = 1
  } else { #it is not significant
    #print(paste("NOT SIGNIFICANT"))
    sig_day_7_from_day_0[p] = 0 
  }
  

  if (lam_vals[p] < 1) { #these are just the zero readings in the first count
    if (end_time_vals[p] <= count_2_thresh) {
      sig_day_7_from_day_0[p] = 0
      start_day_0_count = start_day_0_count + 1
    }
  }
  # 
  # if (end_time_vals[p] < 1) { #these are just the zero readings in the second count
  #   if (lam_vals[p] <= count_1_thresh) {
  #     sig_day_7_from_day_0[p] = 0
  #     start_day_7_count = start_day_7_count + 1
  #   }
  # }

}

# print(paste("sig_day_7_from_day_0",sig_day_7_from_day_0))
# print(paste("first_timepoint_vals",lam_vals))
# print(paste("end_time_vals",end_time_vals))

#change the entries in lam_vals and end_time_vals from 0 to 0.5 so that they appear in the log plot
lam_vals <- replace(lam_vals, lam_vals==0, 0.5)
end_time_vals <- replace(end_time_vals, end_time_vals==0, 0.5)

#print(paste("sig_day_7_from_day_0",sig_day_7_from_day_0))

#make graph of TCR change over time
#asp = 1 makes axes equal. logic in xlim and ylim to choose the maximum of the x and y data. col is conditional colouring.
plot(log(lam_vals), log(end_time_vals), main = "log TCR amount day 0 to day 7",
     xlab = "log(TCR per million day 0)", ylab = "log(TCR per million day 7)",asp = 1,
     xlim = c(log(min(lam_vals)), if_else(max(lam_vals) > max(end_time_vals), log(max(lam_vals)), log(max(end_time_vals)))),
     ylim = c(log(min(end_time_vals)), if_else(max(lam_vals) > max(end_time_vals), log(max(lam_vals)), log(max(end_time_vals)))),
     col = ifelse(sig_day_7_from_day_0 == 1,'red','blue'))
abline(a=0,b=1)

