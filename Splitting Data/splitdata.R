
######################################################################
# 区分训练集和测试集！

library(tidyverse)
library(stringr)
library(gridExtra)
library(future)
library(sva)
library(e1071)
library(pROC)
library(ROCit)
library(caret)
library(doParallel)
library(cancerclass)
rt=read.delim('expression.txt', header=T, sep="\t", check.names=F, row.names=1)
gene=read.table("postive.txt")
gene=gene$V1
rt = rt[gene,]
trainsample = read.delim("sample.txt",header = F,row.names = 1)
rt = as.data.frame(t(rt))
identical(row.names(rt),row.names(trainsample))
train = cbind(trainsample,rt)
colnames(train)[1] = "group"


# 划分数据集
set.seed(123)
combat = train
trainIndex <- createDataPartition(combat$group, p = .7,  ## 80% training set; 20% validation set
                                  list = FALSE, 
                                  times = 1)

## Dataset setting
training <- combat[trainIndex,]
test  <- combat[-trainIndex,]


save(training,file = "training.RData")
save(test,file = "test.RData")


# 存储train
out=cbind(id=row.names(training),training)
colnames(out)[1:3]
outsample = out[,c("id","group")]
write.table(outsample,file="trainsample.txt",sep="\t",row.names=F,quote=F)
out$group = NULL
out$id = NULL
out = as.data.frame(t(out))
out=cbind(id=row.names(out),out)
write.table(out,file="Trainexp.txt",sep="\t",row.names=F,quote=F)


# 存储test
out=cbind(id=row.names(test),test)
colnames(out)[1:3]
outsample = out[,c("id","group")]
write.table(outsample,file="testsample.txt",sep="\t",row.names=F,quote=F)
out$group = NULL
out$id = NULL
out = as.data.frame(t(out))
out=cbind(id=row.names(out),out)
write.table(out,file="Testexp.txt",sep="\t",row.names=F,quote=F)


