library(stringr)
library(data.table)
library(ROSE)
library(randomForest)
library(ROCR)
library(h2o)
library(readr)
# library(minfi)
library(GenomicRanges)
library(plyr)
library(tictoc)
library(h2o)
localH2O<-h2o.init(nthreads = 6,max_mem_size = "32G", startH2O = T, ip = 'localhost', port = 5432)
h2o.removeAll()

# all_data_1[is.na(all_data_1)] = 0

load("../02_Preparation/all_data_1_GM12878.Rdata")
region_GR_GM12878<-region_GR
all_data_1_GM12878<-all_data_1

all_data_1 <- all_data_1_GM12878

#Select training and testing data
train_data_1<- all_data_1[which(region_GR$train==1),] #70%
test_data_1 <- all_data_1[which(region_GR$train==0),] #30%


tfbs <- colnames(train_data_1)[str_detect(colnames(train_data_1),'tfbs')]
# tss <- colnames(train_data_1)[str_detect(colnames(train_data_1),'tss')]
# to_cent <- colnames(train_data_1)[str_detect(colnames(train_data_1),'to_cent')]
me <- c("me",colnames(train_data_1)[startsWith(colnames(train_data_1),'me_R_')|startsWith(colnames(train_data_1),'me_L_')])
# TF <- c("CTCF","POLR2A","EP300","RAD21","MAX","GABPA","SIN3A","ELF1","EGR1","SRF","PML","SMARCA5")
# ChIPseq <- setdiff(colnames(train_data_1),c(tfbs,tss,to_cent,me,"bdsGR","bds","id","chr","start","end","me"))


y="bds"
# x=c('chr',tfbs,tss,to_cent,me,ChIPseq)
x=c('chr', tfbs, me)


data_h2o<-as.h2o(train_data_1)
test_h2o<-as.h2o(test_data_1)

gbm <- h2o.gbm(x=x,y=y,training_frame = data_h2o, balance_classes = T, ntrees = 500, max_depth = 10,nfolds = 5,distribution="bernoulli")
vp <- h2o.varimp(gbm)
perf <- h2o.performance(model=gbm,newdata=test_h2o)
h2o.auc(gbm)

stats <- gbm@model$training_metrics@metrics$max_criteria_and_metric_scores
save(gbm,vp, perf, stats, file="m_compare_GM12878_training.Rdata")
model_path <- h2o.saveModel(object=gbm, path=getwd(), force=TRUE)
print(model_path)

pred.GBM<-h2o.performance(model=gbm,newdata=test_h2o)
jpeg(paste("GM12878.jpg",sep=''))
plot(pred.GBM,type="roc")
dev.off()

save(train_data_1,test_data_1, file="train_test_data_GM12878.Rdata")