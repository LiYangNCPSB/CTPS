# setwd("E:/lung_GateID/lung_gateid_celseq_analysis_20200804/GateID-master")
source('Functions_gate_design.R', echo=TRUE)

library(caret)
library(e1071)
library(pROC)

#_________PART 1: Gate design_________#

##load gateID training dataset (rows= single cells & columns=ClusterID followed by index data in all available channels)
data <- as.data.frame(read.csv("../4-gateID-datacreate/celseq.train.csv", header = T, sep = ",", stringsAsFactors = F, row.names = 1))
colnames(data)[1] <- "ClusterID"

data[,8:25] <- log10(data[,8:25])
data[is.na(data)] <- 0

# gate.show.id <- c("FITC.BB515.A","BUV395.A")
# 
# plot.data <- data[,c("ClusterID",gate.show.id)]
# colnames(plot.data)[2]
# ggplot(plot.data, aes_string(colnames(plot.data)[2],colnames(plot.data)[3])) + geom_point(aes_string(colour = "ClusterID"))


data[data[,1] == "Artery Cell" | data[,1] == "Vein",1] = "ArtVe"
data[data[,1] != "ArtVe" ,1] = "other"
label <- as.factor(data[,1])
data.train <- as.matrix(data[,2:25])
data.train[is.infinite(data.train)] <- 0
data.train <- as.data.frame(data.train)

feature_auc <- list()
for (i in 1:24){
  model <- svm(x= data.train[i], y=label, kernel = )
  pre_svm <- predict(model,newdata = data.train[i])
  svm_roc <- roc(label,as.numeric(pre_svm))
  feature_auc[names(data.train)[i]] <- svm_roc$auc
}


set.seed(1)
lmProfile <- rfe(data.train, label,
                 sizes = 6,
                 rfeControl = rfeControl(functions = rfFuncs))

lmProfile$optVariables
lmProfile$results

write.csv(lmProfile$optVariables,"ArtVe_feature.csv")
