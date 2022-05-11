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

#the location where the data is
data_path<-"C:/Research_Assistant/work/data/TCR_data/NS085/"
#data_path<-"C:/Research_Assistant/work/data/TCR_data/small/"

# Read in a data set just to check it works
#file1<-read.csv(paste0(data_path,"dcr_HVO_439679_day7_1_alpha.csv"))
#print(file1)

#get a list of the files
# filenames = list.files(path=data_path, pattern=NULL, all.files=FALSE,full.names=FALSE)
# print(paste("filenames", filenames))

#get a list of all of the file paths
files <- list.files(path = data_path,pattern = ".csv")
#print(paste("files", files))

# #read in all of the files using read.csv function
# temp <- lapply(paste0(data_path,files),read.csv,sep=",")
# print(paste("temp", temp))
# 
# #add a new column to each of the dataframes which contains the name of the file
# temp2 <- for(i in seq_along(temp)) temp[[i]] = cbind(File=files[i],temp[[i]])
# print(paste("temp2", temp2))
# 
# #make one data table from a list of many
# data <- rbindlist( temp2 )
# 
# print(paste("data", data))

#alternative method
result <- rbindlist(sapply(paste0(data_path,files),read.csv,simplify = FALSE), idcol = 'filename')
# print(str_sub("adios", 3))

#only want the filename (not the filepath)
result$filename <- str_sub(result$filename, 48)


write.csv(result,"C:\\Research_Assistant\\work\\data\\TCR_data\\NS085\\result.csv", row.names = FALSE)
