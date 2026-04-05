source('Functions_gate_design.R', echo=TRUE)

library(caret)
library(e1071)
library(pROC)

#parameter settings
target_cell_type = "Artery Cell"

##load gateID training dataset (rows= single cells & columns=ClusterID followed by index data in all available channels)
data <- as.data.frame(read.csv("../1-CTPS_exaple_data/celseq.train.csv", header = T, sep = ",", stringsAsFactors = F, row.names = 1))
colnames(data)[1] <- "ClusterID"

#Perform log transformation on the data from all channels other than FSC and SSC
data[,8:25] <- log10(data[,8:25])

data[is.na(data)] <- 0


data[data[,1] != target_cell_type ,1] = "other"
label <- as.factor(data[,1])
data.train <- as.matrix(data[,2:25])
data.train[is.infinite(data.train)] <- 0
data.train <- as.data.frame(data.train)

set.seed(123)
#default feature number is 6, can change ‘size’ to select more or less features.
lmProfile <- rfe(data.train, label,
                 sizes = 6, 
                 rfeControl = rfeControl(functions = rfFuncs))

write.csv(lmProfile$optVariables,"Selected_feature.csv")
