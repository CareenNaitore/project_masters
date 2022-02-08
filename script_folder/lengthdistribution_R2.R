########################################
#### Summarizing data mean, sd, se #####
######################################## 
#loading packages in r 
library(dplyr)
library(ggplot2)
library(tidyr)
library(plyr)
# reset the r program 
rm(list=ls())
#Set working directory.
setwd("/home/icipe/Documents/careenwork/objectiveoneactivity/sequencelengthdistribution/R_analysis/")
#loading data 
data <- read.delim("/home/icipe/Documents/careenwork/objectiveoneactivity/sequencelengthdistribution/R_analysis/lengthdistributionlifestages.csv", sep = "\t")
#reshape data 
data_long <- gather(data, condition, measurement, Adfem1:Wtmal3, factor_key=TRUE)
#view file
data_long
#attaching file to R 
attach(data_long)

cdata <- ddply(data_long, c("Length"), summarise, N = length(measurement), mean = mean(measurement), sd = sd(measurement), se = sd / sqrt(N))
# Create a box plott of Absolute reads vs Readlength, using the 'mean' column to determine the distribution of the length in each file and using the error bars to represent variations.
X<-ggplot(cdata, aes(x= Length, y= mean)) +
  scale_x_continuous("Readlength(in nucleotide)", labels = as.character(data$Length), breaks = data$Length) + geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(.9))+ theme_bw() 
X
# adding  labels
X + ggtitle("lifestages Read length distribution") +theme(plot.title = element_text(hjust = 0.5))+ xlab("Readlength(in nucleotide)") + ylab("Absolute number of Reads") 

