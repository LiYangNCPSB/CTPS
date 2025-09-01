# setwd("E:/lung_GateID/lung_gateid_celseq_analysis_20200804/GateID-master")
# source('Functions_gate_design.R', echo=TRUE)
library(ggplot2)
library(ggpubr)
library("sp")
#_________PART 1: Gate design_________#

##load gateID training dataset (rows= single cells & columns=ClusterID followed by index data in all available channels)
data <- as.data.frame(read.csv("../4-gateID-datacreate/celseq.train.csv", header = T, sep = ",", stringsAsFactors = F, row.names = 1))
colnames(data)[1] <- "ClusterID"

data[,8:25] <- log10(data[,8:25])
data[is.na(data)] <- 0

# gate.show.id <- c("PE.Cy5.A","FITC.BB515.A")
gate.show.id <- c("PE.CF594.A","BV510.A")

plot.data <- data[,c("ClusterID",gate.show.id)]

# plot.data <- plot.data[plot.data$ClusterID=="Club" | plot.data$ClusterID=="Ciliated",]

# plot.data <- plot.data[!plot.data$ClusterID == "other",]

# gate1 <- as.data.frame(matrix(0, ncol = 2 , nrow = 5))
# gate1[,1] <- c(1.5, 1.5, 2.8, 2.8, 1.5)
# gate1[,2] <- c(1.0, 2.0, 2.1, 1.0, 1.2)
# colnames(gate1)<-c(colnames(plot.data)[2],colnames(plot.data)[3])
# 
# gate1.res <- point.in.polygon(point.x = plot.data[,2], point.y = plot.data[,3], pol.x = gate1[,1], pol.y = gate1[,2])
# 
# plot.data$select = "no"
# plot.data[gate1.res==1,"select"] = "yes"

p1 <- ggplot(plot.data, aes_string(colnames(plot.data)[2],colnames(plot.data)[3])) + 
  geom_point(aes_string(colour = "ClusterID"))
p1
# 
# plot.data.in.gate1 <- plot.data[gate1.res==1,]
# 
# p2 <- ggplot(plot.data.in.gate1, aes_string(colnames(plot.data)[4],colnames(plot.data)[5])) + 
#   geom_point(aes_string(colour = "ClusterID", shape = "select"))
# 
# ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
