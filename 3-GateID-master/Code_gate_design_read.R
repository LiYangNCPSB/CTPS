# setwd("E:/lung_GateID/lung_gateid_celseq_analysis_20200804/GateID-master")
source('Functions_gate_design.R', echo=TRUE)


#_________PART 1: Gate design_________#

##load gateID training dataset (rows= single cells & columns=ClusterID followed by index data in all available channels)
data <- as.data.frame(read.csv("../4-gateID-datacreate/celseq.train.csv", header = T, sep = ",", stringsAsFactors = F, row.names = 1))
datastart <- 2  #index data start at column 2 since colunm 1 is ClusterID
colnames(data)[1] <- "ClusterID"

data[data <= 0] = 0.1
data[,8:25] <- log10(data[,8:25])
data[is.na(data)] <- 0
# data[data[,1] == "Artery Cell" | data[,1] == "Vein",1] = "ArtVe"

channel.read <- read.csv("../5-GateID-featureselecte/gCap_feature.csv", row.names = 1)
channel.selected <- channel.read$x[1:6]
data <- data[,c("ClusterID",channel.selected)]

##define the desired cell type 
desired.cluster <- "gCap"
tgc <-t(combn(seq(datastart,dim(data)[2]), m = 2)) 
bool.desired <- data$ClusterID %in% desired.cluster

##define parameters for gate design
seed <- 23
thresh <- 0.3 #yield threshold
nverts <- 4 #define number of vertices for gates

ngc <- t(combn(seq(1,dim(tgc)[1]), m = 2))

##create dataframe for GateID results
gate_sol <- data.frame(matrix(0, ncol = 20 , nrow = nrow(ngc)))

gate_sol <- read.csv("gate_result/gate_gCap_log.csv", row.names = 1)

##picks the vertex coordinates (16 values) for the best gate combination

rowgate <- 1 #change the rowgate parameter to visualize different gates solutions
bestgate <- as.numeric(gate_sol[rowgate, c(5:ncol(gate_sol))])
numgates <- 2
ji <- seq(1,length(bestgate),length(bestgate)/numgates)
nd <- length(bestgate)/numgates
gate.ori <- list()
gate.ori <- lapply(ji, function(y) as.matrix(orderCW(bestgate[y:(y+(nd-1))],4)))
gate.ori <- lapply(gate.ori, function(y) as.matrix(rbind(y, y[1, ]))) #contains predicted gate coordinates
h.index <- as.numeric(gate_sol[rowgate, c(1,2)])
channel.comb <- t(combn(seq(datastart,dim(data)[2]), m = 2))
dims <- channel.comb[h.index, ] #channels in which the gates are designed (refers to the columns of the training which contain clusterID)

cal.gate.efficacy(trainingdata = data, bool.desired = bool.desired, gates = gate.ori, dims = dims, plot = T, verbose = T) #plot gate and print gate stats

# for (i in 1:nrow(ngc)){
#   rowgate <- i #change the rowgate parameter to visualize different gates solutions
#   bestgate <- as.numeric(gate_sol[rowgate, c(5:ncol(gate_sol))])
#   numgates <- 2
#   ji <- seq(1,length(bestgate),length(bestgate)/numgates)
#   nd <- length(bestgate)/numgates
#   gate.ori <- list()
#   gate.ori <- lapply(ji, function(y) as.matrix(orderCW(bestgate[y:(y+(nd-1))],4)))
#   gate.ori <- lapply(gate.ori, function(y) as.matrix(rbind(y, y[1, ]))) #contains predicted gate coordinates
#   h.index <- as.numeric(gate_sol[rowgate, c(1,2)])
#   channel.comb <- t(combn(seq(datastart,dim(data)[2]), m = 2))
#   dims <- channel.comb[h.index, ] #channels in which the gates are designed (refers to the columns of the training which contain clusterID)
# 
#   pdf(file = paste("cap_gate_figure/",rowgate,".pdf"))
#   cal.gate.efficacy(trainingdata = data, bool.desired = bool.desired, gates = gate.ori, dims = dims, plot = T, verbose = T) #plot gate and print gate stats
#   dev.off()
# }



