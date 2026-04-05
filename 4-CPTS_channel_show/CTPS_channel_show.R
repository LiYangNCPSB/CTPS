# setwd("E:/lung_GateID/lung_gateid_celseq_analysis_20200804/GateID-master")
source('Functions_gate_design.R', echo=TRUE)
library(ggplot2)
library(ggpubr)
library("sp")
#_________PART 1: Gate design_________#

##load gateID training dataset (rows= single cells & columns=ClusterID followed by index data in all available channels)
data <- as.data.frame(read.csv("../1-CTPS_exaple_data/celseq.train.csv", header = T, sep = ",", stringsAsFactors = F, row.names = 1))
colnames(data)[1] <- "ClusterID"

data[,8:25] <- log10(data[,8:25])
data[is.na(data)] <- 0

plot.data <- data

tgc <-t(combn(seq(2,dim(data)[2]), m = 2))

for (i in c(1:nrow(tgc))){
  p1 <- ggplot(plot.data, aes_string(colnames(plot.data)[tgc[i,1]],colnames(plot.data)[tgc[i,2]])) + 
    geom_point(aes_string(colour = "ClusterID"))
  plot.path <- paste("plot/",colnames(plot.data)[tgc[i,1]],"_",colnames(plot.data)[tgc[i,2]],".jpg", sep = "")
  ggsave(plot.path,p1)
}


