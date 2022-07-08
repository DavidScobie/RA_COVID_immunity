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
input_data_chop_all<-read.csv(paste0(data_path,"input_data_chop_all.csv"))

input_data_chop_all_wide <- reshape(input_data_chop_all, idvar=c('junction_aa','subject','chain','PCR_positive'),timevar='day',direction='wide')
