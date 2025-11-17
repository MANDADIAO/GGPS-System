######################
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
library(mlr3verse)

# GLO gen-set

gene=read.table("hubgene.txt")
gene=gene$V1

# 
rt=read.delim('Trainexp.txt', header=T, sep="\t", check.names=F, row.names=1)
rt = rt[gene,]
trainsample = read.delim("trainsample.txt",header = F,row.names = 1)
rt = as.data.frame(t(rt))
identical(row.names(rt),row.names(trainsample))
train = cbind(trainsample,rt)
colnames(train)[1] = "group"

# 
test=read.delim('Testexp.txt', header=T, sep="\t", check.names=F, row.names=1)
test = test[gene,]
testsample = read.delim("testsample.txt",header = F,row.names = 1)
test = as.data.frame(t(test))
identical(row.names(test),row.names(testsample))
test = cbind(testsample,test)
colnames(test)[1] = "group"


training = train
colnames(training)[1] = "group"
colnames(test)[1] = "group"
train = training
test = test

## label setting 
table(train$group)
train$group <- ifelse(train$group == 'P','1','0') %>% factor(.,levels = c('0','1'))
test$group <- ifelse(test$group == 'P','1','0') %>% factor(.,levels = c('0','1'))
train$group=factor(train$group)
test$group=factor(test$group)

positive_class<-"1"

task_train = as_task_classif(
  train,
  target = "group",
  positive=positive_class
)
task_test= as_task_classif(
  test,
  target = "group",
  positive = positive_class
)

# learners
learners = list(
  learner_logreg = lrn("classif.log_reg", predict_type = "prob",
                       predict_sets = c("train", "test")),
  learner_xgboost = lrn("classif.xgboost", predict_type = "prob", 
                        predict_sets = c("train", "test")),
  learner_lda = lrn("classif.lda", predict_type = "prob",
                    predict_sets = c("train", "test")),
  learner_svm = lrn("classif.svm", predict_type = "prob",
                    predict_sets = c("train", "test")),
  learner_nb = lrn("classif.naive_bayes", predict_type = "prob",
                   predict_sets = c("train", "test")),
  learner_knn = lrn("classif.kknn", scale = FALSE,
                    predict_type = "prob",predict_sets = c("train", "test")),
  learner_rpart = lrn("classif.rpart",
                      predict_type = "prob",predict_sets = c("train", "test")),
  learner_rf = lrn("classif.ranger", num.trees = 1000,
                   predict_type = "prob",predict_sets = c("train", "test"))
)


#### 
learners
final_model=learners$learner_rf
#RF
final_model$train(task_train)
#
prediction <- final_model$predict(task_test)

autoplot(prediction,type = "roc",auc=T)

prediction_tab=as.data.table(prediction)
pdf(file = "TestRoc.pdf",width = 4,height = 4)
pROC::plot.roc(test$group,prediction_tab$prob.1,print.auc=T)
dev.off()

save(final_model,file = "final_model.RData")

