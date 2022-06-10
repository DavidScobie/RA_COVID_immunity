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

###find the lower and upper bounds for 95% of the data for poisson centered on 1st time point

##first, add a column (top lim) to the data frame which is 5x the first timepoint (surely this is above the 95% range)
subset_pat_439679_alpha_prod_wide <- subset_pat_439679_alpha_prod_wide %>% 
  mutate(top_lim = if_else(duplicate_count.0 > 0.7 , 5*(duplicate_count.0), 0))  #Deal with the fake abundance of 0.5 if not present

##second, integrate poisson distribution centered at 1st time point between 0 to top lim to find (area_0_to_top_lim)
#(No point in doing this as the integral comes out as 1 in every case)

##third, starting at the mean, integrate the poisson dist between 0 to each point higher at increments of 1. 
#stop doing this when integral = 95% (0 to top lim area). The upper limit is the upper significance bound.

#############  This code is integrating poisson distributions of means (lam_vals) from 0 up to (upper_vals)
############# it seems that for means much above 20, the integration doesnt seem to work.
integrand<-function(x,lam)(((lam^x)*(exp(-lam)))/(factorial(x)))  #poisson distribution
  tmpfun <- function(lam,upper) {                                 #function to integrate
    integrate(integrand,lower=0,upper=upper,lam=lam)$value
  }

lam_vals <- c(1,4,2,5,3,3,2,3,2,20)
#lam_vals = subset_pat_439679_alpha_prod_wide$duplicate_count.0[0:10]
for (i in 1:length(lam_vals)) {
  print(paste("This iteration represents range value", i))
  print(paste("lam_vals[i]", lam_vals[i]))
  #upper_vals = linspace(lam_vals[i], 5*lam_vals[i], n = (5*lam_vals[i])-(lam_vals[i])+1)
  upper_vals = linspace(0, 5*lam_vals[i], n = (5*lam_vals[i])-(0)+1)
  print(paste("upper_vals", upper_vals))
  print(sapply(lam_vals[i], function(x) { sapply(upper_vals, function(y) tmpfun(x,y))}))
}
 ################## 

#


# this is just a quick way to do specific integrals if need be
integrand<-function(x,lam)(((lam^x)*(exp(-lam)))/(factorial(x)))        ######works
tmpfun <- function(lam) {
  integrate(integrand,lower=0,upper=4,lam=lam)$value
}
sapply(3,tmpfun)

############ summing a function between limits
foo1 = function(k, start = 3, stop = 5) {
  result = 0
  for (i in start:stop) {
    result = result + i / stop * cov(k, stats::lag(k, k = i))
  }
  return(result)
} 

k = rnorm(100)

stops = 5:20
res = sapply(stops, function(b) foo1(k, start = 3, stop = b))
print(res)

############# my attempt at summing a function between certain limits
sum_poisson = function(lam, start=0, stop=5) {
  result = 0
  for (i in start:stop) {
    result = result + ((lam**(i))*(exp(-lam))/(factorial(i)))
  }
  return(result)
}

lam_vals <- c(1,4,2,5,3,3,2,3,2,20)
for (p in 1:length(lam_vals)) {
  stops = linspace(0, 5*lam_vals[p], n = (5*lam_vals[p])-(0)+1)
  summation = sapply(lam_vals[p], function(x) { sapply(stops, function(y) sum_poisson(x, start=0, stop=y))})
  print(summation)
}

#fourth, repeat all this but start integral at (top lim) and work along left from mean of dist to find lower sig limit.


#add column for day 0 to day 7 significant?
subset_pat_439679_alpha_prod_wide <- subset_pat_439679_alpha_prod_wide %>% 
  mutate(sig_day_7_from_day_0 = if_else(duplicate_count.7 < duplicate_count.0 - sqrt(duplicate_count.0) | duplicate_count.7 > duplicate_count.0 + sqrt(duplicate_count.0) , 1, 0))

#deal with the cases of 0's and 1's where we want not significant
subset_pat_439679_alpha_prod_wide <- subset_pat_439679_alpha_prod_wide %>% 
  mutate(sig_day_7_from_day_0_deal_with_zeros = if_else(duplicate_count.0 == 0 & duplicate_count.7 == 1, 0, sig_day_7_from_day_0))

#make graph of TCR change over time
#asp = 1 makes axes equal. logic in xlim and ylim to choose the maximum of the x and y data. col is conditional colouring.
plot(log(subset_pat_439679_alpha_prod_wide$duplicate_count.0), log(subset_pat_439679_alpha_prod_wide$duplicate_count.7), main = "log TCR amount day 0 to day 7",
     xlab = "log(TCR per million day 0)", ylab = "log(TCR per million day 7)",asp = 1,
     xlim = c(log(min(subset_pat_439679_alpha_prod_wide$duplicate_count.0)), if_else(max(subset_pat_439679_alpha_prod_wide$duplicate_count.0) > max(subset_pat_439679_alpha_prod_wide$duplicate_count.7), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.0)), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.7)))),
     ylim = c(log(min(subset_pat_439679_alpha_prod_wide$duplicate_count.7)), if_else(max(subset_pat_439679_alpha_prod_wide$duplicate_count.0) > max(subset_pat_439679_alpha_prod_wide$duplicate_count.7), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.0)), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.7)))), 
     col = ifelse(subset_pat_439679_alpha_prod_wide$sig_day_7_from_day_0_deal_with_zeros == 1,'red','blue'))
abline(a=0,b=1)

#######day 7 to day 14 significance

#add column for day 7 to day 14 significant?
subset_pat_439679_alpha_prod_wide <- subset_pat_439679_alpha_prod_wide %>% 
  mutate(sig_day_14_from_day_7 = if_else(duplicate_count.14 < duplicate_count.7 - sqrt(duplicate_count.7) | duplicate_count.14 > duplicate_count.7 + sqrt(duplicate_count.7) , 1, 0))

#deal with the cases of 0's and 1's where we want not significant
subset_pat_439679_alpha_prod_wide <- subset_pat_439679_alpha_prod_wide %>% 
  mutate(sig_day_14_from_day_7_deal_with_zeros = if_else(duplicate_count.7 == 0 & duplicate_count.14 == 1, 0, sig_day_14_from_day_7))

#make graph of TCR change over time
plot(log(subset_pat_439679_alpha_prod_wide$duplicate_count.7), log(subset_pat_439679_alpha_prod_wide$duplicate_count.14), main = "log TCR amount day 7 to day 14",
     xlab = "log(TCR per million day 7)", ylab = "log(TCR per million day 14)",asp = 1,
     xlim = c(log(min(subset_pat_439679_alpha_prod_wide$duplicate_count.7)), if_else(max(subset_pat_439679_alpha_prod_wide$duplicate_count.7) > max(subset_pat_439679_alpha_prod_wide$duplicate_count.14), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.7)), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.14)))),
     ylim = c(log(min(subset_pat_439679_alpha_prod_wide$duplicate_count.14)), if_else(max(subset_pat_439679_alpha_prod_wide$duplicate_count.7) > max(subset_pat_439679_alpha_prod_wide$duplicate_count.14), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.7)), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.14)))), 
     col = ifelse(subset_pat_439679_alpha_prod_wide$sig_day_14_from_day_7_deal_with_zeros == 1,'red','blue'))
abline(a=0,b=1)

#######day 0 to day 14 significance

#add column for day 0 to day 14 significant?
subset_pat_439679_alpha_prod_wide <- subset_pat_439679_alpha_prod_wide %>% 
  mutate(sig_day_14_from_day_0 = if_else(duplicate_count.14 < duplicate_count.0 - sqrt(duplicate_count.0) | duplicate_count.14 > duplicate_count.0 + sqrt(duplicate_count.0) , 1, 0))

#deal with the cases of 0's and 1's where we want not significant
subset_pat_439679_alpha_prod_wide <- subset_pat_439679_alpha_prod_wide %>% 
  mutate(sig_day_14_from_day_0_deal_with_zeros = if_else(duplicate_count.0 == 0 & duplicate_count.14 == 1, 0, sig_day_14_from_day_0))

#make graph of TCR change over time
plot(log(subset_pat_439679_alpha_prod_wide$duplicate_count.0), log(subset_pat_439679_alpha_prod_wide$duplicate_count.14), main = "log TCR amount day 0 to day 14",
     xlab = "log(TCR per million day 0)", ylab = "log(TCR per million day 14)",asp = 1,
     xlim = c(log(min(subset_pat_439679_alpha_prod_wide$duplicate_count.0)), if_else(max(subset_pat_439679_alpha_prod_wide$duplicate_count.0) > max(subset_pat_439679_alpha_prod_wide$duplicate_count.14), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.0)), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.14)))),
     ylim = c(log(min(subset_pat_439679_alpha_prod_wide$duplicate_count.14)), if_else(max(subset_pat_439679_alpha_prod_wide$duplicate_count.0) > max(subset_pat_439679_alpha_prod_wide$duplicate_count.14), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.0)), log(max(subset_pat_439679_alpha_prod_wide$duplicate_count.14)))), 
     col = ifelse(subset_pat_439679_alpha_prod_wide$sig_day_14_from_day_0_deal_with_zeros == 1,'red','blue'))
abline(a=0,b=1)




