#for taking .csv files from the cluster and plotting the results on this pc
options(stringsAsFactors=FALSE)
require(gdata)
require(dplyr)
require(ggplot2)
#library(XLConnect)
#library(xlsx)
library(data.table)
library(stringr)
library(tidyr)
library(pracma)
library(MASS)

#read the data frames in
data_path = "C:/Research_Assistant/work/data/TCR_data/significant_TCR_csvs/individual_pats/" 
sig_TCRs <- read.csv(paste0(data_path,"sig_TCRs_P_10exp(-7)_PCR_pos_alph_day0_day7_439679.csv"))
all_pats_wide_plotting<-read.csv(paste0(data_path,"all_pats_wide_plotting_P_10exp(-7)_PCR_pos_alph_day0_day7_439679.csv"))
all_pats_wide_chopped_sig_lines_plotting <- read.csv(paste0(data_path,"all_pats_wide_chopped_sig_lines_plotting_P_10exp(-7)_PCR_pos_alph_day0_day7_439679.csv"))
all_pats_wide_sig_lines_plotting <- read.csv(paste0(data_path,"all_pats_wide_sig_lines_plotting_P_10exp(-7)_PCR_pos_alph_day0_day7_439679.csv"))


p1 <- ggplot(all_pats_wide_plotting) #define the dataframe to plot
p2 <- p1 + geom_point(aes(log_day_0_for_plot_with_jitter_u_p_m, log_day_7_for_plot_with_jitter_u_p_m), shape = all_pats_wide_plotting$sig_day_7_from_day_0, colour = all_pats_wide_plotting$sig_day_7_from_day_0+1) + xlab("log(TCR per million day 0)") + ylab("log(TCR per million day 14)") #make the scatterplot and give axis labels
p3 <- p2 + geom_line(data=all_pats_wide_chopped_sig_lines_plotting, aes(x=chopped_log_ordered_lam_vals_used_unique_u_p_m,y=chopped_log_ordered_low_sig_lim_list_u_p_m),linetype="dashed", color = 'blue', size = 2) #plot the chopped lower significance boundary line
p4 <- p3 + geom_line(data=all_pats_wide_sig_lines_plotting, aes(x=log_ordered_lam_vals_used_unique_u_p_m,y=log_ordered_high_sig_lim_list_u_p_m),linetype="dashed", color = 'blue', size = 2) #plot the upper significance boundary line
p5 <- p4 + coord_cartesian(xlim(c(min(all_pats_wide_plotting$log_day_0_for_plot_with_jitter_u_p_m), if_else(max(all_pats_wide_plotting$log_day_0_for_plot_with_jitter_u_p_m) > max(all_pats_wide_plotting$log_day_7_for_plot_with_jitter_u_p_m), max(all_pats_wide_plotting$log_day_0_for_plot_with_jitter_u_p_m), max(all_pats_wide_plotting$log_day_7_for_plot_with_jitter_u_p_m))))) #define the xlim of the plot, coordinates cartesian ensures that we plot all of the line
p6 <- p5 + coord_cartesian(ylim=c(min(all_pats_wide_plotting$log_day_7_for_plot_with_jitter_u_p_m), if_else(max(all_pats_wide_plotting$log_day_0_for_plot_with_jitter_u_p_m) > max(all_pats_wide_plotting$log_day_7_for_plot_with_jitter_u_p_m), max(all_pats_wide_plotting$log_day_0_for_plot_with_jitter_u_p_m), max(all_pats_wide_plotting$log_day_7_for_plot_with_jitter_u_p_m)))) #define the ylim of the plot, coordinates cartesian ensures that we plot all of the line
p7 <- p6 + geom_abline() #plot the line y=x
p7 #show the plot