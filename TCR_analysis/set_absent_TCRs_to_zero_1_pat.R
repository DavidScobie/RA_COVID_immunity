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
data_path<-"C:/Research_Assistant/work/data/TCR_data/NS085/"
# pat_439679_day_0_alpha<-read.csv(paste0(data_path,"dcr_HVO_439679_pre_1_alpha.csv"))
# pat_439679_day_7_alpha<-read.csv(paste0(data_path,"dcr_HVO_439679_day7_1_alpha.csv"))
# pat_439679_day_14_alpha<-read.csv(paste0(data_path,"dcr_HVO_439679_day14_1_alpha.csv"))
pat_439679_day_0_alpha<-read.csv(paste0(data_path,"dcr_HVO_635729_pre_1_alpha.csv"))
pat_439679_day_7_alpha<-read.csv(paste0(data_path,"dcr_HVO_635729_day7_1_alpha.csv"))
pat_439679_day_14_alpha<-read.csv(paste0(data_path,"dcr_HVO_635729_day14_1_alpha.csv"))

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

# select variables junction_aa, duplicate_count.0, duplicate_count.7, duplicate_count.14
myvars <- c("junction_aa", "duplicate_count.0", "duplicate_count.7", "duplicate_count.14")
subset_pat_439679_alpha_prod_wide <- pat_439679_alpha_prod_wide[myvars]
#subset_pat_439679_alpha_prod_wide <- pat_439679_alpha_prod_wide[myvars][1:500,] #only want first 500 for speed

#replace duplicate count NA with duplicate count = 0
subset_pat_439679_alpha_prod_wide[is.na(subset_pat_439679_alpha_prod_wide)] <- 0

rownames(subset_pat_439679_alpha_prod_wide) <- NULL  #This resets the index of the rows of the dataframe

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
#lam_vals <- subset_pat_439679_alpha_prod_wide$duplicate_count.0[1:500] #just taking a subset
#lam_vals <- c(0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2)

#the timepoint 7 dataset
end_time_vals <- subset_pat_439679_alpha_prod_wide$duplicate_count.7
#end_time_vals <- subset_pat_439679_alpha_prod_wide$duplicate_count.7[1:500] #just taking a subset
#end_time_vals <- c(3,4,5,6,7,8,3,4,5,6,7,8,3,4,5,6,7,8,)

####Cannot have a poisson dist centred at 0. Therefore take significance boundary for poisson centred at 1 and take same limit. Make sig line straight at LHS of graph.
ar_bef_1 <- vector()
count2 <- 0
for (z in 1:max(end_time_vals)) {  #from 1 up to the maximum of the day 7 count values
  ar_bef_1 <- sapply(1, function(x) sum_poisson(x,start=0,stop=z))  #this is evluating poisson of mean 1 starting at 0 and ending at all possible values
  if (ar_bef_1 > (1-(p_value/2))) {   #this is for determining all end_time_vals which would theortically be above high_sig_lim
    count2 = count2+1
  }
}
count_2_thresh <- max(end_time_vals) - count2 #this finds the smallest possible end_time_vals value which is above high_sig_lim (6 would be a reasonable value)
#######

#the binary array to see if significant or not
sig_day_7_from_day_0 <- vector()

manual_x_threshold = 50 #what value of x do we start approximating sig with (mean +- (num_std_devs*sqrt(lambda))) ?
top_x_threshold = 180 #what value of x do we make the significance lines flat?
#what value of x do you want the straight line on the LHS beginning for the upper sig bound to go up to. Surely just 0 and 1? Hence 1.2 value
bottom_threshold = 1 #not right to have this as a decimal because it means we are considering a poisson dist centred not on an integer

low_sig_lim_list <- vector() #initialise list for plotting lower significance bound
high_sig_lim_list <- vector()
lam_vals_used_unique <- vector()

for (p in 1:length(lam_vals)) {
  greater_than_min_sig <- vector() #initialise empty array for greater than minimum significance level
  low_counter <- 0 #initialize lower significance level counter
  lesser_than_max_sig <- vector() #initialise empty array for less than max significance level
  high_counter <- 0 #initialize higher significance level counter

  greater_than_up_sig_bound <- vector() #error check vector for upper significance bound
  sig_error_check_counter <- 0 #counter for the error checker on upper significance bound

  #logic to find upper significance level
  if (lam_vals[p] <= bottom_threshold) { #the beginning cut off with straight significance lines

    num_lam_multip = 8  #how many multiples of lambda do we want the summation of poisson to go up to. This is important for memory reasons

    stops = linspace(0, num_lam_multip*bottom_threshold, n = (num_lam_multip*bottom_threshold)-(0)+1)  #need stops to go far beyond
    summation = sapply(bottom_threshold, function(x) { sapply(stops, function(y) sum_poisson(x, start=0, stop=y))})
    for (k in 1:length(summation)) {  #for loop to check each value in the array of summation
      
      if (is.nan(summation[k])) {
        print(paste0("summation",summation,"lam_vals[p]",lam_vals[p],"loop 1"))
      }
      
      if (summation[k] > (p_value/2)) {   #logic for finding the lower significance level
        low_counter = low_counter + 1
        greater_than_min_sig[low_counter] = k
      }
      if (summation[k] < (1-(p_value/2))) {  #logic for finding the upper significance level
        high_counter = high_counter + 1
        lesser_than_max_sig[high_counter] = k
      }

      if (summation[k] > (1-(p_value/2))) { #logic to check if any values in summation are above the upper sig bound
        sig_error_check_counter = sig_error_check_counter + 1  #we need the code to go through this loop at least once, or there is an error. We are not finding upper significance threshold properly
      }
    }
    #if no values in summation above upper sig bound throw error because using wrong indicie for significance
    if (sig_error_check_counter < 1) {   #we want to find the lam_vals[p] value that has broken it
      print(paste0("summation",summation,"lam_vals[p]",lam_vals[p],"loop 1"))
    }
    if( sig_error_check_counter < 1 ) stop('NOT CORRECTLY FINDING UPPER SIGNIFICANCE THRESHOLD. Raise num_lam_multip to fix this. Or increase p_value.')

    #low_sig_lim = min(greater_than_min_sig) -1 #this is the proper way, but this gives Inf as greater_than_min_sig is empty for small lam_vals. So have low_sig_lim=0
    low_sig_lim = 0  #as we are just dealing with the first few lambda values, the low_sig_lim will be zero.
    high_sig_lim = max(lesser_than_max_sig) -1 #we need to allow this to be 0 if necessary. Need -1 to work with indices as k in 1:length(summation)

    ###########################


  }

  else if (lam_vals[p] < manual_x_threshold) {     #lambda is low, so manually sum to find significance level

    #this loop is in to save memory. We only need large num_lam_multip for small lambdas (in order to find upper p value threshold). Can always add steps in the loop to be more memory efficient
    if (lam_vals[p] < 10) {
      num_lam_multip = 8  #number of multiples of lambda to sum up to in the poisson distribution
    } else {
      num_lam_multip = 3
    }

    stops = linspace(0, num_lam_multip*lam_vals[p], n = (num_lam_multip*lam_vals[p])-(0)+1)  #need stops to go far beyond
    summation = sapply(lam_vals[p], function(x) { sapply(stops, function(y) sum_poisson(x, start=0, stop=y))})
    for (k in 1:length(summation)) {  #for loop to check each value in the array of summation
      
      if (is.nan(summation[k])) {
        print(paste0("summation",summation,"lam_vals[p]",lam_vals[p],"loop 2"))
      }
        
      if (summation[k] > (p_value/2)) {   #logic for finding the lower significance level
        low_counter = low_counter + 1
        greater_than_min_sig[low_counter] = k
      }
      if (summation[k] < (1-(p_value/2))) {  #logic for finding the upper significance level
        high_counter = high_counter + 1
        lesser_than_max_sig[high_counter] = k
      }

      if (summation[k] > (1-(p_value/2))) { #logic to check if any values in summation are above the upper sig bound
        sig_error_check_counter = sig_error_check_counter + 1
      }
    }
    #if no values in summation above upper sig bound throw error because using wrong indicie for significance
    if (sig_error_check_counter < 1) {   #we want to find the lam_vals[p] value that has broken it
      print(paste0("summation",summation,"lam_vals[p]",lam_vals[p],"loop 2"))
    }
    if( sig_error_check_counter < 1 ) stop('NOT CORRECTLY FINDING UPPER SIGNIFICANCE THRESHOLD. Raise num_lam_multip to fix this. Or increase p_value.')

    low_sig_lim = min(greater_than_min_sig) -1 #we need to allow this to be 0 if necessary. Need -1 to work with indices as k in 1:length(summation)
    high_sig_lim = max(lesser_than_max_sig) -1 #we need to allow this to be 0 if necessary. Need -1 to work with indices as k in 1:length(summation)

  } else if (lam_vals[p] < top_x_threshold) {  #lambda is big so approximate as normal dist. use qnorm to find the number of standard deviations from the p value
    num_std_devs_gauss <- qnorm(1-(p_value/2))
    low_sig_lim <- ceiling(lam_vals[p] - (num_std_devs_gauss*sqrt(lam_vals[p])))
    high_sig_lim <- floor(lam_vals[p] + (num_std_devs_gauss*sqrt(lam_vals[p])))

  } else { #lambda is really big so have the significance lines go straight from this point onwards
    num_std_devs_gauss <- qnorm(1-(p_value/2))
    low_sig_lim <- ceiling(top_x_threshold - (num_std_devs_gauss*sqrt(top_x_threshold)))
    high_sig_lim <- ((10**100)*max(lam_vals))  #This is a huge number which makes the significance line go almost directly up
  }

  ######### making the low_sig_lim and high_sig_lim into lists ready to plot as dashed lines

  #for each lam_vals value (day 0 abundance value), append low_sig_lim to a list and high_sig_lim to a list
  if (p == 1) {  #the vector is of length zero in the beginning so we have to append in the first instance
    low_sig_lim_list <- append(low_sig_lim_list, low_sig_lim)
    high_sig_lim_list <- append(high_sig_lim_list, high_sig_lim)
  }
  #we only want to append low_sig_lim and high_sig_lim to lists if they are new values. Otherwise we will have too many points for plotting
  else if (setequal(lam_vals_used_unique[!(lam_vals_used_unique %in% lam_vals[p])], lam_vals_used_unique)) {   #here we are checking if lam_vals[p] appears in lam_vals_used_unique. If it does then we do not do the else if loop. If it does not then we append to list
    low_sig_lim_list <- append(low_sig_lim_list, low_sig_lim)
    high_sig_lim_list <- append(high_sig_lim_list, high_sig_lim)
  }

  lam_vals_used_unique <- unique(append(lam_vals_used_unique, lam_vals[p])) #we append the lam_vals value to this list each time to check whether low_sig_lim and high_sig_lim need to be appended onto their lists in future iterations

  ########## comparing the end_time_vals[p] value to the significance thresholds to determine significance and hence if sig_day_7_from_day_0 = 0 or 1

  #compare the end timepoint value to the significance thresholds
  if (end_time_vals[p] < low_sig_lim) { #it is significant at lower end
    sig_day_7_from_day_0[p] = 1
  } else if (end_time_vals[p] > high_sig_lim) { #it is significant at higher end
    sig_day_7_from_day_0[p] = 1
  } else { #it is not significant
    sig_day_7_from_day_0[p] = 0
  }

  ##########  for zero count in day 0 the significance checker comes out as always non-significant, so we have this to check the values properly

  if (lam_vals[p] < 1) { #these are just the zero readings in the first count
    if (end_time_vals[p] <= count_2_thresh) {
      sig_day_7_from_day_0[p] = 0
    }
  }

}

#we need to order the arrays of low_sig_lim_list, high_sig_lim_list and lam_vals_used_unique for plotting
ordered_low_sig_lim_list <- low_sig_lim_list[order(lam_vals_used_unique)]
ordered_high_sig_lim_list <- high_sig_lim_list[order(lam_vals_used_unique)]
ordered_lam_vals_used_unique <- sort(lam_vals_used_unique)

#replace the -Inf at the start of ordered_high_sig_lim_list with the value at indicie of 2. This is because high_sig_lim list is undefined when lam_vals[p]=0 as cannot have poisson centred at zero.
if (ordered_lam_vals_used_unique[1] == 0) {  #Only if first value in ordered_lam_vals_used_unique is 0. This is only time the high_sig_lim value is undefined
  ordered_high_sig_lim_list[1] = ordered_high_sig_lim_list[2] #we just call the high_sig_lim value the value when lam_vals[p]=2 as this is defined and we want same significance threshold
}


#change the entries in lam_vals and end_time_vals from 0 to 0.5 so that they appear in the log plot
lam_vals_used_unique <- replace(lam_vals_used_unique, lam_vals_used_unique==0, 0.5)
lam_vals <- replace(lam_vals, lam_vals==0, 0.5)
end_time_vals <- replace(end_time_vals, end_time_vals==0, 0.5)

#do the same for the lines I am plotting
ordered_lam_vals_used_unique <- replace(ordered_lam_vals_used_unique, ordered_lam_vals_used_unique==0, 0.5)
ordered_low_sig_lim_list <- replace(ordered_low_sig_lim_list, ordered_low_sig_lim_list==0, 0.5)
ordered_high_sig_lim_list <- replace(ordered_high_sig_lim_list, ordered_high_sig_lim_list==0, 0.5)

#find total number of TCR's on each day
total_TCRs_day_0 <- sum(as.numeric(subset_pat_439679_alpha_prod_wide$duplicate_count.0))
total_TCRs_day_7 <- sum(as.numeric(subset_pat_439679_alpha_prod_wide$duplicate_count.7))     #This is required for using units of TCR per million in the plot
total_TCRs_day_14 <- sum(as.numeric(subset_pat_439679_alpha_prod_wide$duplicate_count.14))

#need to scale the lam_vals and end_time_vals so that units are in TCR per million
lam_vals_u_p_m <- ((10**6)/(total_TCRs_day_0))*lam_vals
end_time_vals_u_p_m <- ((10**6)/(total_TCRs_day_7))*end_time_vals

#take all the important arrays and place them into a dataframe. First the points and significance binary
subset_pat_439679_alpha_prod_wide$log_day_0_for_plot_u_p_m <- log(lam_vals_u_p_m)
subset_pat_439679_alpha_prod_wide$log_day_7_for_plot_u_p_m <- log(end_time_vals_u_p_m)
subset_pat_439679_alpha_prod_wide$sig_day_7_from_day_0 <- sig_day_7_from_day_0

###introducing jitter, so that the points dont overlap, meaning we cant ascertain the number of points in each place on the plot
subset_pat_439679_alpha_prod_wide <- subset_pat_439679_alpha_prod_wide %>% rowwise() %>%
  mutate(log_day_0_for_plot_with_jitter_u_p_m = rnorm(1,mean=log_day_0_for_plot_u_p_m, sd=1/(log_day_0_for_plot_u_p_m + 10)))   #the SD is important in determining the amount of jitter

subset_pat_439679_alpha_prod_wide <- subset_pat_439679_alpha_prod_wide %>% rowwise() %>%
  mutate(log_day_7_for_plot_with_jitter_u_p_m = rnorm(1,mean=log_day_7_for_plot_u_p_m, sd=1/(log_day_7_for_plot_u_p_m + 10)))   #the SD is important in determining the amount of jitter

#need to scale the significance lines so that they too are in units of TCR per million (not actual count)
ordered_low_sig_lim_list_u_p_m <- ((10**6)/(total_TCRs_day_7))*ordered_low_sig_lim_list
ordered_high_sig_lim_list_u_p_m <- ((10**6)/(total_TCRs_day_7))*ordered_high_sig_lim_list
ordered_lam_vals_used_unique_u_p_m <- ((10**6)/(total_TCRs_day_0))*ordered_lam_vals_used_unique

#Put significance lines into a big data frame for plotting
log_ordered_lam_vals_used_unique_u_p_m <- log(ordered_lam_vals_used_unique_u_p_m)
log_ordered_low_sig_lim_list_u_p_m <- log(ordered_low_sig_lim_list_u_p_m)
log_ordered_high_sig_lim_list_u_p_m <- log(ordered_high_sig_lim_list_u_p_m)
subset_pat_439679_alpha_prod_wide_sig_lines <- data.frame(log_ordered_lam_vals_used_unique_u_p_m, log_ordered_low_sig_lim_list_u_p_m, log_ordered_high_sig_lim_list_u_p_m)

#we dont want to plot the start of the lower significance line, as it looks incorrect with the jitter. Therefore we chop off the start of the array which are all the points which have value zero
chopped_log_ordered_low_sig_lim_list_u_p_m <- log_ordered_low_sig_lim_list_u_p_m[max(which(log(ordered_low_sig_lim_list_u_p_m) == min(log(ordered_low_sig_lim_list_u_p_m)))):length(log_ordered_low_sig_lim_list_u_p_m)]
chopped_log_ordered_lam_vals_used_unique_u_p_m <- log_ordered_lam_vals_used_unique_u_p_m[max(which(log(ordered_low_sig_lim_list_u_p_m) == min(log(ordered_low_sig_lim_list_u_p_m)))):length(log_ordered_low_sig_lim_list_u_p_m)]
subset_pat_439679_alpha_prod_wide_chopped_sig_lines <- data.frame(chopped_log_ordered_low_sig_lim_list_u_p_m, chopped_log_ordered_lam_vals_used_unique_u_p_m)  #put the chopped values into a new dataframe for plotting

p1 <- ggplot(subset_pat_439679_alpha_prod_wide) #define the dataframe to plot
p2 <- p1 + geom_point(aes(log_day_0_for_plot_with_jitter_u_p_m, log_day_7_for_plot_with_jitter_u_p_m), shape = sig_day_7_from_day_0, colour = sig_day_7_from_day_0+1) + xlab("log(TCR per million day 0)") + ylab("log(TCR per million day 7)") #make the scatterplot and give axis labels
p3 <- p2 + geom_line(data=subset_pat_439679_alpha_prod_wide_chopped_sig_lines, aes(x=chopped_log_ordered_lam_vals_used_unique_u_p_m,y=chopped_log_ordered_low_sig_lim_list_u_p_m),linetype="dashed", color = 'blue', size = 2) #plot the chopped lower significance boundary line
p4 <- p3 + geom_line(data=subset_pat_439679_alpha_prod_wide_sig_lines, aes(x=log_ordered_lam_vals_used_unique_u_p_m,y=log_ordered_high_sig_lim_list_u_p_m),linetype="dashed", color = 'blue', size = 2) #plot the upper significance boundary line
p5 <- p4 + coord_cartesian(xlim(c(min(subset_pat_439679_alpha_prod_wide$log_day_0_for_plot_with_jitter_u_p_m), if_else(max(lam_vals) > max(end_time_vals), max(subset_pat_439679_alpha_prod_wide$log_day_0_for_plot_with_jitter_u_p_m), max(subset_pat_439679_alpha_prod_wide$log_day_7_for_plot_with_jitter_u_p_m))))) #define the xlim of the plot, coordinates cartesian ensures that we plot all of the line
p6 <- p5 + coord_cartesian(ylim=c(min(subset_pat_439679_alpha_prod_wide$log_day_7_for_plot_with_jitter_u_p_m), if_else(max(lam_vals) > max(end_time_vals), max(subset_pat_439679_alpha_prod_wide$log_day_0_for_plot_with_jitter_u_p_m), max(subset_pat_439679_alpha_prod_wide$log_day_7_for_plot_with_jitter_u_p_m)))) #define the ylim of the plot, coordinates cartesian ensures that we plot all of the line
p7 <- p6 + geom_abline() #plot the line y=x
p7 #show the plot

