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

data[data[,1] == "Artery Cell" | data[,1] == "Vein",1] = "ArtVe"

channel.read <- read.csv("../5-GateID-featureselecte/ArtVe_feature.csv", row.names = 1)
channel.selected <- channel.read$x[1:6]
# channel.selected <- c("PE.Cy7.A", "BV510.A", "APC.A", "BV750.A", "PerCP.Cy5.5.A", "BUV395.A")
data <- data[,c("ClusterID",channel.selected)]

##define the desired cell type 
desired.cluster <- "ArtVe"
tgc <-t(combn(seq(datastart,dim(data)[2]), m = 2)) 
bool.desired <- data$ClusterID %in% desired.cluster

##define parameters for gate design
seed <- 23
thresh <- 0.3 #yield threshold
nverts <- 4 #define number of vertices for gates

ngc <- t(combn(seq(1,dim(tgc)[1]), m = 2))

##create dataframe for GateID results
gate_sol <- data.frame(matrix(0, ncol = 20 , nrow = nrow(ngc)))

##the for loop runs the gate design for each gate combination (runs in 5 minutes) 
##for gate design on a full dataset, run on an hpc and submit as array job
##gate_sol contains solutions

for(j in c(1:nrow(ngc))){
  print(j)
  gc <- ngc[j,]
  gs <- as.vector(t(tgc[gc, ]))
  dimensions <- gs
  dataset <- data
  numgates <- length(gs)/2
  d2.desired <- dataset[bool.desired,dimensions]
  d2.nondesired <- dataset[!bool.desired,dimensions]
  jumpindex <- seq(1,ncol(d2.desired), 2)
  d2.desired <- lapply(jumpindex, function(x) as.matrix(d2.desired[, c(x,(x+1))]))
  d2.nondesired <- lapply(jumpindex, function(x) as.matrix(d2.nondesired[, c(x,(x+1))]))
  l = lapply(1:length(d2.desired), function(x) t(apply(d2.desired[[x]], 2, function(y) quantile(y, probs = c(0.01, 0.99)))))
  quart <- as.matrix(l[[1]])
  for (i in 2:length(l)) {quart <- rbind(quart, l[[i]])}
  initials <- lapply(l, function(y) orderCW(as.matrix(expand.grid(y[1,], y[2,])), 4))
  initials <- unlist(lapply(initials, function(x) as.numeric(x)))
  numpoints.thresh <- round(dim(d2.desired[[1]])[1] * thresh)
  nd <- length(initials)/numgates
  ji <- seq(1,length(initials),length(initials)/numgates)
  lb <-  as.numeric(unlist(lapply(quart[,1], function(y) rep(y,4))))
  ub <-  as.numeric(unlist(lapply(quart[,2], function(y) rep(y,4))))
  lb2 <-  rep(0, length(lb))
  
  ##malsChains-based optimization
  initialpop <- initials
  cmares <- malschains(li, lower = lb, upper = ub, maxEvals = 20000, seed = seed, initialpop = initialpop, verbosity = 0,
                      control = malschains.control(optimum = 0, istep = 100, ls = "sw", effort = 0.55, alpha = 0.5, popsize = 50))
  impurity <- cmares$fitness
  initialpop <- cmares$sol
  cmares2  <- malschains(my, lower = lb2, upper = ub, maxEvals = 10000, initialpop = initialpop, seed = seed, verbosity = 0,
                       control = malschains.control(istep = 300, ls = "cmaes", effort = 0.55, alpha = 0.5, popsize = 50))
  yield.cma <- abs(cmares2$fitness)
  temp = sum(bool.desired)*yield.cma
  impurity.cma <- temp/(impurity+temp)
  spec = c(gc, yield.cma, impurity.cma, cmares2$sol)
  
  gate_sol[j,] <- spec
}

##order gates based on yield (column 3) and purity (column 4)
gate_sol <- gate_sol[rev(order(gate_sol$X4, gate_sol$X3)), ]

# write.csv(gate_sol,"gate_result/gate.ciliated.csv")

##picks the vertex coordinates (16 values) for the best gate combination
rowgate <-1 #change the rowgate parameter to visualize different gates solutions
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

write.csv(gate_sol,"gate_result/gate_ArtVe_log.csv")

for (i in 1:nrow(ngc)){
  rowgate <- i #change the rowgate parameter to visualize different gates solutions
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

  pdf(file = paste("ArtVe_gate_figure/",rowgate,".pdf"))
  cal.gate.efficacy(trainingdata = data, bool.desired = bool.desired, gates = gate.ori, dims = dims, plot = T, verbose = T) #plot gate and print gate stats
  dev.off()
}
