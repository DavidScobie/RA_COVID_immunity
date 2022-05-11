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
library(fread)

data_path<-"C:/Research_Assistant/work/data/TCR_data/NS085/"

# Read in meta data for labels
file1<-read.csv(paste0(data_path,"dcr_HVO_439679_day7_1_alpha.csv"))

#print(file1)

#get a list of the files
filenames = list.files(path=data_path, pattern=NULL, all.files=FALSE,full.names=FALSE)
print(paste("filenames", filenames))

#get all of the file paths
files <- list.files(path = data_path,pattern = ".csv")
print(paste("files", files))

temp <- lapply(paste0(data_path,filenames),read.csv,sep=",")
temp <- for(i in seq_along(temp)) temp[[i]] = cbind(File=files[i],temp[[i]])
data <- rbindlist( temp )

print(paste("data", data))


