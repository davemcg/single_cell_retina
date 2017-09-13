

#Impute categories for remaining ~30,000 cells not included in classification
library(here)
library(caret)
library(doMC)
registerDoMC(cores=4)
set.seed(1234)
load(here('data/retina_seurat_subSet.Rdata'))
trainIndex <- createDataPartition(retina@meta.data$res.1, p=0.5, times=1, list=T)
#macosko_trainIndex <- trainIndex <- createDataPartition(retina@meta.data$Macosko_Clusters, p=0.5, times=1, list=T)
retinaTrain <- t(retina@scale.data[,trainIndex$Resample1])
retinaTest <- t(retina@scale.data[,-trainIndex$Resample1])
outcomesTrain <- as.factor(retina@meta.data$res.1[trainIndex$Resample1])
outcomesTest <- as.factor(retina@meta.data$res.1[-trainIndex$Resample1])
myControl <- trainControl(method = "repeatedcv", repeats=5, number = 10)

# ~80 accuracy
rf_mod <- train(x=retinaTrain, y=outcomesTrain, method = 'rf', trControl = myControl)
# ~92% accuracy
svm_mod <- train(x=retinaTrain, y=outcomesTrain, method = 'svmLinear', trControl = myControl, tuneGrid = data.frame(.C = c(0.0001, 0.001, 0.01, 0.1)))
save(svm_mod, file=here('data/svm_imputation_mod.Rdata'))
#LogitBoost_mod <- train(x=retinaTrain, y=outcomesTrain, method = 'LogitBoost', trControl = myControl)
# ~20% accuracy
bayes_mod <- train(x=retinaTrain, y=outcomesTrain, method = 'bayesglm', trControl = myControl)
# 
gbm_mod <- train(x=retinaTrain, y=outcomesTrain, method = 'gbm', trControl = myControl)
#keras_net_mod <- train(x=retinaTrain, y=outcomesTrain, trControl = myControl, method = 'mlpKeras')
#naive_bayes_mod <-  train(x=retinaTrain, y=outcomesTrain, method = 'naive_bayes', trControl = myControl)


#svm works the best so far (still want to try keras neural net)
pred <- predict(svm_mod, retinaTest)
postResample(pred=pred, obs=outcomesTest)

# keras/Tensorflow attempt
registerDoMC(cores=1)
keras_net_mod <- train(x=retinaTrain, y=outcomesTrain, trControl = myControl, method = 'mlpKerasDropout')
save(keras_net_mod, file=here('data/keras_net_mod.Rdata'))

pred <- predict(svm_mod, retinaTest)
postResample(pred=pred, obs=outcomesTest)