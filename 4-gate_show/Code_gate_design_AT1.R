# setwd("E:/lung_GateID/lung_gateid_celseq_analysis_20200804/GateID-master")
source('Functions_gate_design.R', echo=TRUE)
library(ggplot2)
library(ggpubr)
library("sp")
#_________PART 1: Gate design_________#

##load gateID training dataset (rows= single cells & columns=ClusterID followed by index data in all available channels)
data <- as.data.frame(read.csv("../4-gateID-datacreate/celseq.train.csv", header = T, sep = ",", stringsAsFactors = F, row.names = 1))
colnames(data)[1] <- "ClusterID"

data[,8:25] <- log10(data[,8:25])
data[is.na(data)] <- 0

# gate.show.id <- c("APC.A","BUV395.A","BUV496.A","BUV395.A","BUV496.A","BV421.A","APC.A","PE.A")

gate.show.id <- c("mCFP.A", "PE.Cy7.A", "BV711.A", "PE.Cy7.A","PE.A","BUV395.A")

plot.data <- data[,c("ClusterID",gate.show.id)]

# plot.data <- plot.data[plot.data$ClusterID=="Club" | plot.data$ClusterID=="Ciliated",]

# plot.data <- plot.data[!plot.data$ClusterID == "other",]

gate1 <- as.data.frame(matrix(0, ncol = 2 , nrow = 5))
gate1[,1] <- c(3.431728163, 3.085535548, 3.89697894, 4.485363693, 3.431728163)
gate1[,2] <- c(1.198229151, 2.005051897, 3.714291054, 3.323021989, 1.198229151)

# gate1[,1] <- c(3.5, 3.5, 5, 5, 3.5)
# gate1[,2] <- c(2.5, 4, 4, 2.5, 2.5)

colnames(gate1)<-c(colnames(plot.data)[2],colnames(plot.data)[3])
gate1.res <- point.in.polygon(point.x = plot.data[,2], point.y = plot.data[,3], pol.x = gate1[,1], pol.y = gate1[,2])

gate2 <- as.data.frame(matrix(0, ncol = 2 , nrow = 5))
# gate2[,1] <- c(2.763548373, 4.140295612, 0.349467613, 1.118899799, 2.763548373)
# gate2[,2] <- c(1.051133065,3.755027484,3.123088303,3.270514546,1.051133065)
colnames(gate2)<-c(colnames(plot.data)[4],colnames(plot.data)[5])

gate2.res <- point.in.polygon(point.x = plot.data[,4], point.y = plot.data[,5], pol.x = gate2[,1], pol.y = gate2[,2])


plot.data$select = "no"
plot.data[gate1.res==1 & gate2.res==1,"select"] = "yes"

plot.data.in.gate <- plot.data[which(plot.data$select == "yes"),]

p1 <- ggplot(plot.data, aes_string(colnames(plot.data)[2],colnames(plot.data)[3])) + 
  geom_point(aes_string(colour = "ClusterID", shape = "select")) +
  geom_polygon(data = gate1, alpha = 0.2)

p2 <- ggplot(plot.data, aes_string(colnames(plot.data)[4],colnames(plot.data)[5])) + 
  geom_point(aes_string(colour = "ClusterID", shape = "select"))+
  geom_polygon(data = gate2, alpha = 0.2)

p3 <- ggplot(plot.data, aes_string(colnames(plot.data)[6],colnames(plot.data)[7])) + 
  geom_point(aes_string(colour = "ClusterID", shape = "select"))

# p4 <- ggplot(plot.data.in.gate, aes_string(colnames(plot.data)[8],colnames(plot.data)[9])) + 
#   geom_point(aes_string(colour = "ClusterID", shape = "select"))

ggarrange(p1, p2, p3, ncol=2, nrow=2, common.legend = F, legend="bottom")

