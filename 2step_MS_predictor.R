## Implementation scripts for 2-step MS relapse predictor paper
## Author: Yuri Ahuja
## Last updated: Aug, 2020

library(foreach)
library(doParallel)
library(glmnet)
library(nlme)
library(lmmen)
library(lme4)
library(stats)
library(stringr)
source("library_v2.R")  ###
source("library_v3.r")  ###
source("Library.R")     ###AUC
library(ggplot2)
library(latex2exp)

expit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}


## Read in and order data

# 24 month window, 3 month subwindow
CC_all <- readRDS('CC_all_tf24_tr12withrelapse.rds')
CC_all <- CC_all[order(CC_all$Period),]; CC_all <- CC_all[order(CC_all$PatientNum),]
CC_all <- CC_all[!is.na(CC_all$DISEASE_DURA),]
CC_all$StartDate <- as.Date(CC_all$StartDate)

CC_val_separate <- readRDS('validation_set24_tr12.rds')
CC_val_separate <- CC_val_separate[order(CC_val_separate$Period),]
CC_val_separate <- CC_val_separate[order(CC_val_separate$PatientNum),]
CC_val_separate <- CC_val_separate[!is.na(CC_val_separate$DISEASE_DURA),]

# 12 month window, 3 month subwindow
CC_all2 <- readRDS('CC_all_tf12_tr12withrelapse.rds')
CC_all2 <- CC_all2[order(CC_all2$Period),]; CC_all2 <- CC_all2[order(CC_all2$PatientNum),]
CC_all2 <- CC_all2[!is.na(CC_all2$DISEASE_DURA),]
CC_all2$StartDate <- as.Date(CC_all2$StartDate)

CC_val2_separate <- readRDS('validation_set12_tr12.rds')
CC_val2_separate <- CC_val2_separate[order(CC_val2_separate$Period),]
CC_val2_separate <- CC_val2_separate[order(CC_val2_separate$PatientNum),]
CC_val2_separate <- CC_val2_separate[!is.na(CC_val2_separate$DISEASE_DURA),]


## Compute age, demographics

demographics <- read.csv('ms_demographics.csv')[,1:6]
CC_all$Age <- sapply(1:nrow(CC_all),function(i){
  id <- CC_all$PatientNum[i]
  match <- which(demographics$PATIENT_NUM == id)
  (as.numeric(CC_all$StartDate[i])-demographics$BIRTH_DATE[match]+25569)/365.25
})
CC_all2$Age <- sapply(1:nrow(CC_all2),function(i){
  id <- CC_all2$PatientNum[i]
  match <- which(demographics$PATIENT_NUM == id)
  (as.numeric(CC_all2$StartDate[i])-demographics$BIRTH_DATE[match]+25569)/365.25
})
CC_val_separate$Age <- sapply(1:nrow(CC_val_separate),function(i){
  id <- CC_val_separate$PatientNum[i]
  match <- which(demographics$PATIENT_NUM == id)
  (as.numeric(CC_val_separate$StartDate[i])-demographics$BIRTH_DATE[match]+25569)/365.25
})
CC_val2_separate$Age <- sapply(1:nrow(CC_val2_separate),function(i){
  id <- CC_val2_separate$PatientNum[i]
  match <- which(demographics$PATIENT_NUM == id)
  (as.numeric(CC_val2_separate$StartDate[i])-demographics$BIRTH_DATE[match]+25569)/365.25
})

intersection <- intersect(CC_all$PatientNum, CC_val_separate$PatientNum)
CC_all <- CC_all[!(CC_all$PatientNum %in% intersection),]

CC_val_separate_CLIMB <- CC_val_separate[CC_val_separate$PatientNum %in% intersection,]
CC_val_separate_nonCLIMB <- CC_val_separate[!(CC_val_separate$PatientNum %in% intersection),]

intersection <- intersect(CC_all2$PatientNum, CC_val2_separate$PatientNum)
CC_all2 <- CC_all2[!(CC_all2$PatientNum %in% intersection),]

CC_val2_separate_CLIMB <- CC_val2_separate[CC_val2_separate$PatientNum %in% intersection,]
CC_val2_separate_nonCLIMB <- CC_val2_separate[!(CC_val2_separate$PatientNum %in% intersection),]


## Final experiment with bootstrap

logfile <- "2step_experiment.txt"
writeLines(c(""), file(logfile,'w'))
clust <- makeCluster(20, outfile=logfile)
registerDoParallel(clust)

evaluate <- function(actual,prediction,T){
  result = ROC.Est.FUN(actual, prediction, yy0=0.5, fpr0=0.05)
  fscore = 2/(1/result$TPR + 1/result$PPV)
  
  prediction_calibrated = expit(cbind(1,prediction) %*% glm(actual~prediction,family='binomial')$coefficients)
  predround = as.integer(prediction_calibrated > mean(actual))
  
  #cutoff = quantile(prediction,1-prev)
  #predround = as.integer(prediction>cutoff)
  sens = sum(predround*actual) / sum(actual)
  spec = sum((1-predround)*(1-actual)) / sum(1-actual)
  ppv = sum(predround*actual) / sum(predround)
  npv = sum((1-predround)*(1-actual)) / sum(1-predround)
  fscorePrev = 2/(1/sens + 1/ppv)
  
  ord = order(T)
  actual = actual[ord]; prediction = prediction[ord]; T = T[ord]
  predround = unlist(sapply(min(T):max(T),function(t){
    keep <- which(T >= t-1 & T <= t+1)
    coefs <- glm(actual[keep]~prediction[keep],family='binomial')$coefficients
    prediction_calibrated <- expit(cbind(1,prediction[T==t]) %*% coefs)
    as.integer(prediction_calibrated > mean(actual[keep]))
  }))
  sensT = sum(predround*actual) / sum(actual)
  specT = sum((1-predround)*(1-actual)) / sum(1-actual)
  ppvT = sum(predround*actual) / sum(predround)
  npvT = sum((1-predround)*(1-actual)) / sum(1-predround)
  fscorePrevT = 2/(1/sensT + 1/ppvT)
  
  return(c(result$AUC,sens,spec,ppv,npv,fscorePrev,sensT,specT,ppvT,npvT,fscorePrevT))
}


bootstrap_results <- foreach(it=1:500, .combine=cbind, .packages=c('glmnet','nlme','lmmen','lme4','stats','stringr')) %dopar% {
  print(paste('On iteration',it))
  
  trainPats <- sample(unique(CC_all$PatientNum),length(unique(CC_all$PatientNum)), replace=TRUE)
  trainIndices <- unlist(sapply(trainPats,function(id){which(CC_all$PatientNum == id)}))
  CC_train <- CC_all[trainIndices,]
  CC_train2 <- CC_all2[trainIndices,]

  
  lasso_dems_24 <- cv.glmnet(as.matrix(CC_train[,c('Age','DISEASE_DURA','FEMALE','RACE')]),CC_train$CC,family='binomial',type.measure='auc')
  lasso_dems_pred_24_separate <- predict(lasso_dems_24,as.matrix(CC_val_separate[,c('Age','DISEASE_DURA','FEMALE','RACE')]))
  lasso_dems_24_result_separate <- evaluate(CC_val_separate$CC,lasso_dems_pred_24_separate,CC_val_separate$Period)

  lasso_dems_12 <- cv.glmnet(as.matrix(CC_train2[,c('Age','DISEASE_DURA','FEMALE','RACE')]),CC_train2$CC,family='binomial',type.measure='auc')
  lasso_dems_pred_12_separate <- predict(lasso_dems_12,as.matrix(CC_val2_separate[,c('Age','DISEASE_DURA','FEMALE','RACE')]))
  lasso_dems_12_result_separate <- evaluate(CC_val2_separate$CC,lasso_dems_pred_12_separate,CC_val2_separate$Period)

  
  lasso_dems_phecode_24 <- cv.glmnet(as.matrix(CC_train[,c('Age','DISEASE_DURA','FEMALE','RACE','PheCode.335_')]),CC_train$CC,family='binomial',type.measure='auc')
  lasso_dems_phecode_pred_24_separate <- predict(lasso_dems_phecode_24,as.matrix(CC_val_separate[,c('Age','DISEASE_DURA','FEMALE','RACE','PheCode.335_')]))
  lasso_dems_phecode_24_result_separate <- evaluate(CC_val_separate$CC,lasso_dems_phecode_pred_24_separate,CC_val_separate$Period)

  lasso_dems_phecode_12 <- cv.glmnet(as.matrix(CC_train2[,c('Age','DISEASE_DURA','FEMALE','RACE','PheCode.335_')]),CC_train2$CC,family='binomial',type.measure='auc')
  lasso_dems_phecode_pred_12_separate <- predict(lasso_dems_phecode_12,as.matrix(CC_val2_separate[,c('Age','DISEASE_DURA','FEMALE','RACE','PheCode.335_')]))
  lasso_dems_phecode_12_result_separate <- evaluate(CC_val2_separate$CC,lasso_dems_phecode_pred_12_separate,CC_val2_separate$Period)

  
  lasso_dems_Y_24 <- cv.glmnet(as.matrix(CC_train[,c('Age','DISEASE_DURA','FEMALE','RACE','PRIOR_RELAPSE_24MONS')]),
                            CC_train$CC,family='binomial',type.measure='auc')
  lasso_dems_Y_pred_24_separate <- predict(lasso_dems_Y_24,as.matrix(CC_val_separate[,c('Age','DISEASE_DURA','FEMALE','RACE','PRIOR_RELAPSE_24MONS')]))
  lasso_dems_Y_24_result_separate <- evaluate(CC_val_separate$CC,lasso_dems_Y_pred_24_separate,CC_val_separate$Period)

  lasso_dems_Y_12 <- cv.glmnet(as.matrix(CC_train2[,c('Age','DISEASE_DURA','FEMALE','RACE','PRIOR_RELAPSE_24MONS')]),
                               CC_train2$CC,family='binomial',type.measure='auc')
  lasso_dems_Y_pred_12_separate <- predict(lasso_dems_Y_12,as.matrix(CC_val2_separate[,c('Age','DISEASE_DURA','FEMALE','RACE','PRIOR_RELAPSE_24MONS')]))
  lasso_dems_Y_12_result_separate <- evaluate(CC_val2_separate$CC,lasso_dems_Y_pred_12_separate,CC_val2_separate$Period)


  omit <- which(colnames(CC_train) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS','PRIOR_RELAPSE_24MONS'))
  
  lasso_dems_EHR_24 <- cv.glmnet(as.matrix(CC_train[,-c(1:9,omit)]),CC_train$CC,family='binomial',type.measure='auc')
  lasso_dems_EHR_pred_24_separate <- predict(lasso_dems_EHR_24,as.matrix(CC_val_separate[,-c(1:9,omit)]))
  lasso_dems_EHR_24_result_separate <- evaluate(CC_val_separate$CC,lasso_dems_EHR_pred_24_separate,CC_val_separate$Period)

  
  omit <- which(colnames(CC_train2) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS','PRIOR_RELAPSE_24MONS'))
  
  lasso_dems_EHR_12 <- cv.glmnet(as.matrix(CC_train2[,-c(1:9,221,omit)]),CC_train2$CC,family='binomial',type.measure='auc')
  lasso_dems_EHR_pred_12_separate <- predict(lasso_dems_EHR_12,as.matrix(CC_val2_separate[,-c(1:9,221,omit)]))
  lasso_dems_EHR_12_result_separate <- evaluate(CC_val2_separate$CC,lasso_dems_EHR_pred_12_separate,CC_val2_separate$Period)


  omit <- which(colnames(CC_train) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS'))
  
  lasso_all_24 <- cv.glmnet(as.matrix(CC_train[,-c(1:9,omit)]),CC_train$CC,family='binomial',type.measure='auc')
  lasso_all_pred_24_separate <- predict(lasso_all_24,as.matrix(CC_val_separate[,-c(1:9,omit)]))
  lasso_all_24_result_separate <- evaluate(CC_val_separate$CC,lasso_all_pred_24_separate,CC_val_separate$Period)


  omit <- which(colnames(CC_train2) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS'))
  
  lasso_all_12 <- cv.glmnet(as.matrix(CC_train2[,-c(1:9,omit)]),CC_train2$CC,family='binomial',type.measure='auc')
  lasso_all_pred_12_separate <- predict(lasso_all_12,as.matrix(CC_val2_separate[,-c(1:9,omit)]))
  lasso_all_12_result_separate <- evaluate(CC_val2_separate$CC,lasso_all_pred_12_separate,CC_val2_separate$Period)

 
  omit <- which(colnames(CC_train) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS','PRIOR_RELAPSE_24MONS','RH_hat'))
  omit2 <- which(colnames(CC_train2) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS','PRIOR_RELAPSE_24MONS','RH_hat'))
  
  lasso_24_RHhat <- cv.glmnet(as.matrix(CC_train[,-c(1:9,omit)]),log(1+CC_train$PRIOR_RELAPSE_12MONS))
  CC_train$RH_hat <- predict(lasso_24_RHhat,as.matrix(CC_train[,-c(1:9,omit)]),type='response')
  CC_val_separate$RH_hat <- predict(lasso_24_RHhat,as.matrix(CC_val_separate[,-c(1:9,omit)]),type='response')

  lasso_12_RHhat <- cv.glmnet(as.matrix(CC_train2[,-c(1:9,omit2)]),log(1+CC_train2$PRIOR_RELAPSE_12MONS))
  CC_train2$RH_hat <- predict(lasso_12_RHhat,as.matrix(CC_train2[,-c(1:9,omit2)]),type='response')
  CC_val2_separate$RH_hat <- predict(lasso_12_RHhat,as.matrix(CC_val2_separate[,-c(1:9,omit2)]),type='response')

  lasso_RHhat_24 <- cv.glmnet(as.matrix(CC_train[,c('Age','DISEASE_DURA','FEMALE','RACE','RH_hat')]),
                               CC_train$CC,family='binomial',type.measure='auc')
  lasso_RHhat_pred_24_separate <- predict(lasso_RHhat_24,as.matrix(CC_val_separate[,c('Age','DISEASE_DURA','FEMALE','RACE','RH_hat')]))
  lasso_RHhat_24_result_separate_normal <- evaluate(CC_val_separate$CC,lasso_RHhat_pred_24_separate,CC_val_separate$Period)

  lasso_RHhat_12 <- cv.glmnet(as.matrix(CC_train2[,c('Age','DISEASE_DURA','FEMALE','RACE','RH_hat')]),
                              CC_train2$CC,family='binomial',type.measure='auc')
  lasso_RHhat_pred_12_separate <- predict(lasso_RHhat_12,as.matrix(CC_val2_separate[,c('Age','DISEASE_DURA','FEMALE','RACE','RH_hat')]))
  lasso_RHhat_12_result_separate_normal <- evaluate(CC_val2_separate$CC,lasso_RHhat_pred_12_separate,CC_val2_separate$Period)

  omit <- which(colnames(CC_train) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS','PRIOR_RELAPSE_24MONS'))
  omit2 <- which(colnames(CC_train2) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS','PRIOR_RELAPSE_24MONS'))
 
  lasso_EHR_RHhat_24 <- cv.glmnet(as.matrix(CC_train[,-c(1:9,omit)]),CC_train$CC,family='binomial',type.measure='auc')
  lasso_EHR_RHhat_pred_24_separate <- predict(lasso_EHR_RHhat_24,as.matrix(CC_val_separate[,-c(1:9,omit)]))
  lasso_EHR_RHhat_24_result_separate_normal <- evaluate(CC_val_separate$CC,lasso_EHR_RHhat_pred_24_separate,CC_val_separate$Period)

  lasso_EHR_RHhat_12 <- cv.glmnet(as.matrix(CC_train2[,-c(1:9,omit2)]),CC_train2$CC,family='binomial',type.measure='auc')
  lasso_EHR_RHhat_pred_12_separate <- predict(lasso_EHR_RHhat_12,as.matrix(CC_val2_separate[,-c(1:9,omit2)]))
  lasso_EHR_RHhat_12_result_separate_normal <- evaluate(CC_val2_separate$CC,lasso_EHR_RHhat_pred_12_separate,CC_val2_separate$Period)

  CC_val_separate <- CC_val_separate[,-which(colnames(CC_val_separate)=='RH_hat')]
  CC_val2_separate <- CC_val2_separate[,-which(colnames(CC_val2_separate)=='RH_hat')]
  
  c(lasso_dems_24_result_separate,lasso_dems_12_result_separate,lasso_dems_phecode_24_result_separate,lasso_dems_phecode_12_result_separate,
    lasso_dems_Y_24_result_separate,lasso_dems_Y_12_result_separate,lasso_dems_EHR_24_result_separate,lasso_dems_EHR_12_result_separate,
    lasso_all_24_result_separate,lasso_all_12_result_separate,lasso_RHhat_24_result_separate_normal,lasso_RHhat_12_result_separate_normal,
    lasso_EHR_RHhat_24_result_separate_normal,lasso_EHR_RHhat_12_result_separate_normal)
}
bootstrap_results <- array(bootstrap_results,dim=c(11,14,500))
rownames(bootstrap_results) <- c('AUC','SensPrev','SpecPrev','PPVPrev','NPVPrev','FscorePrev',
                                 'SensPrevT','SpecPrevT','PPVPrevT','NPVPrevT','FscorePrevT')
colnames(bootstrap_results) <- c('ASRD 24','ASRD 12','Phecode 24','Phecode 12','RH 24','RH 12','EHR 24','EHR 12',
                                 'RH+EHR 24','RH+EHR 12','RHhat 24','RHhat 12','RHhat+EHR 24','RHhat+EHR 12')
save(bootstrap_results, '2step_bootstrap_results.RData')


## Non-bootstrapped results

prev24 <- mean(CC_val_separate$CC)
prev12 <- mean(CC_val2_separate$CC)

lasso_dems_24 <- cv.glmnet(as.matrix(CC_all[,c('Age','DISEASE_DURA','FEMALE','RACE')]),CC_all$CC,family='binomial',type.measure='auc')
predictions_24_3_dems <- predict(lasso_dems_24,as.matrix(CC_val_separate[,c('Age','DISEASE_DURA','FEMALE','RACE')]),type='response')
evaluate(CC_val_separate$CC,predictions_24_3_dems,CC_val_separate$Period)

lasso_dems_12 <- cv.glmnet(as.matrix(CC_all2[,c('Age','DISEASE_DURA','FEMALE','RACE')]),CC_all2$CC,family='binomial',type.measure='auc')
predictions_12_3_dems <- predict(lasso_dems_12,as.matrix(CC_val2_separate[,c('Age','DISEASE_DURA','FEMALE','RACE')]),type='response')
dems12_result <- evaluate(CC_val2_separate$CC,predictions_12_3_dems,CC_val2_separate$Period)


lasso_phecode_24 <- cv.glmnet(as.matrix(CC_all[,c('Age','DISEASE_DURA','FEMALE','RACE')]),CC_all$CC,family='binomial',type.measure='auc')
predictions_24_3_phecode <- predict(lasso_phecode_24,as.matrix(CC_val_separate[,c('Age','DISEASE_DURA','FEMALE','RACE')]),type='response')
evaluate(CC_val_separate$CC,predictions_24_3_phecode,CC_val_separate$Period)

lasso_phecode_12 <- cv.glmnet(as.matrix(CC_all2[,c('Age','DISEASE_DURA','FEMALE','RACE')]),CC_all2$CC,family='binomial',type.measure='auc')
predictions_12_3_phecode <- predict(lasso_phecode_12,as.matrix(CC_val2_separate[,c('Age','DISEASE_DURA','FEMALE','RACE')]),type='response')
phecode12_result <- evaluate(CC_val2_separate$CC,predictions_12_3_phecode,CC_val2_separate$Period)


lasso_Y_24 <- cv.glmnet(as.matrix(CC_all[,c('Age','DISEASE_DURA','FEMALE','RACE','PRIOR_RELAPSE_12MONS')]),CC_all$CC,family='binomial',type.measure='auc')
predictions_24_3_Y <- predict(lasso_Y_24,as.matrix(CC_val_separate[,c('Age','DISEASE_DURA','FEMALE','RACE','PRIOR_RELAPSE_12MONS')]),type='response')
evaluate(CC_val_separate$CC,predictions_24_3_Y,CC_val_separate$Period)

lasso_Y_12 <- cv.glmnet(as.matrix(CC_all2[,c('Age','DISEASE_DURA','FEMALE','RACE','PRIOR_RELAPSE_12MONS')]),CC_all2$CC,family='binomial',type.measure='auc')
predictions_12_3_Y <- predict(lasso_Y_12,as.matrix(CC_val2_separate[,c('Age','DISEASE_DURA','FEMALE','RACE','PRIOR_RELAPSE_12MONS')]),type='response')
RH12_result <- evaluate(CC_val2_separate$CC,predictions_12_3_Y,CC_val2_separate$Period)


omit <- which(colnames(CC_all) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS','PRIOR_RELAPSE_24MONS'))

lasso_EHR_noY_24 <- cv.glmnet(as.matrix(CC_all[,-c(1:9,omit)]),CC_all$CC,family='binomial',type.measure='auc')
predictions_24_3_EHR_noY <- predict(lasso_EHR_noY_24,as.matrix(CC_val_separate[,-c(1:9,omit)]),type='response')
evaluate(CC_val_separate$CC,predictions_24_3_EHR_noY,CC_val_separate$Period)

omit <- which(colnames(CC_all2) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS','PRIOR_RELAPSE_24MONS'))

lasso_EHR_noY_12 <- cv.glmnet(as.matrix(CC_all2[,-c(1:9,omit)]),CC_all2$CC,family='binomial',type.measure='auc')
predictions_12_3_EHR_noY <- predict(lasso_EHR_noY_12,as.matrix(CC_val2_separate[,-c(1:9,omit)]),type='response')
EHR12_result <- evaluate(CC_val2_separate$CC,predictions_12_3_EHR_noY,CC_val2_separate$Period)

omit <- which(colnames(CC_all) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_24MONS'))

lasso_EHR_24 <- cv.glmnet(as.matrix(CC_all[,-c(1:9)]),CC_all$CC,family='binomial',type.measure='auc')
predictions_24_3_EHR <- predict(lasso_EHR_24,as.matrix(CC_val_separate[,-c(1:9)]),type='response')
evaluate(CC_val_separate$CC,predictions_24_3_EHR,CC_val_separate$Period)

omit <- which(colnames(CC_all2) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_24MONS'))

lasso_EHR_12 <- cv.glmnet(as.matrix(CC_all2[,-c(1:9)]),CC_all2$CC,family='binomial',type.measure='auc')
predictions_12_3_EHR <- predict(lasso_EHR_12,as.matrix(CC_val2_separate[,-c(1:9)]),type='response')
RHEHR12_result <- evaluate(CC_val2_separate$CC,predictions_12_3_EHR,CC_val2_separate$Period)


omit <- which(colnames(CC_all) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS','PRIOR_RELAPSE_24MONS','RH_hat'))

lasso_24_RHhat <- cv.glmnet(as.matrix(CC_all[,-c(1:9,omit)]),log(1+CC_all$PRIOR_RELAPSE_12MONS))
CC_all$RH_hat <- predict(lasso_24_RHhat,as.matrix(CC_all[,-c(1:9,omit)]),type='response')[,1]
CC_val_separate$RH_hat <- predict(lasso_24_RHhat,as.matrix(CC_val_separate[,-c(1:9,omit)]),type='response')[,1]
evaluate(CC_val_separate$PRIOR_RELAPSE_12MONS>0,CC_val_separate$RH_hat,CC_val_separate$Period)
cor(CC_val_separate$PRIOR_RELAPSE_12MONS,CC_val_separate$RH_hat)

omit <- which(colnames(CC_all2) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS','PRIOR_RELAPSE_24MONS','RH_hat'))

lasso_12_RHhat <- cv.glmnet(as.matrix(CC_all2[,-c(1:9,omit)]),log(1+CC_all2$PRIOR_RELAPSE_12MONS))
CC_all2$RH_hat <- predict(lasso_12_RHhat,as.matrix(CC_all2[,-c(1:9,omit)]),type='response')[,1]
CC_val2_separate$RH_hat <- predict(lasso_12_RHhat,as.matrix(CC_val2_separate[,-c(1:9,omit)]),type='response')[,1]
evaluate(CC_val2_separate$PRIOR_RELAPSE_12MONS>0,CC_val2_separate$RH_hat,CC_val2_separate$Period)
cor(CC_val2_separate$PRIOR_RELAPSE_12MONS,CC_val2_separate$RH_hat)


lasso_RHhat_24 <- cv.glmnet(as.matrix(CC_all[,c('Age','DISEASE_DURA','FEMALE','RACE','RH_hat')]),
                            CC_all$CC,family='binomial',type.measure='auc')
lasso_RHhat_pred_24_separate <- predict(lasso_RHhat_24,as.matrix(CC_val_separate[,c('Age','DISEASE_DURA','FEMALE','RACE','RH_hat')]),type='response')
lasso_RHhat_24_result_separate <- evaluate(CC_val_separate$CC,lasso_RHhat_pred_24_separate,CC_val_separate$Period)
print(lasso_RHhat_24_result_separate)

lasso_RHhat_12 <- cv.glmnet(as.matrix(CC_all2[,c('Age','DISEASE_DURA','FEMALE','RACE','RH_hat')]),
                            CC_all2$CC,family='binomial',type.measure='auc')
lasso_RHhat_pred_12_separate <- predict(lasso_RHhat_12,as.matrix(CC_val2_separate[,c('Age','DISEASE_DURA','FEMALE','RACE','RH_hat')]),type='response')
lasso_RHhat_12_result_separate <- evaluate(CC_val2_separate$CC,lasso_RHhat_pred_12_separate,CC_val2_separate$Period)
print(lasso_RHhat_12_result_separate)


omit <- which(colnames(CC_all2) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS','PRIOR_RELAPSE_24MONS'))

lasso_RHhat_EHR_12 <- cv.glmnet(as.matrix(CC_all2[,-c(1:9,omit)]),CC_all2$CC,family='binomial',type.measure='auc')
lasso_RHhat_EHR_pred_12_separate <- predict(lasso_RHhat_EHR_12,as.matrix(CC_val2_separate[,-c(1:9,omit)]),type='response')
lasso_RHhat_EHR_12_result_separate <- evaluate(CC_val2_separate$CC,lasso_RHhat_EHR_pred_12_separate,CC_val2_separate$Period)
print(lasso_RHhat_EHR_12_result_separate)

omit <- which(colnames(CC_all) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS','PRIOR_RELAPSE_24MONS'))

lasso_RHhat_EHR_24 <- cv.glmnet(as.matrix(CC_all[,-c(1:9,omit)]),CC_all$CC,family='binomial',type.measure='auc')
lasso_RHhat_EHR_pred_24_separate <- predict(lasso_RHhat_EHR_24,as.matrix(CC_val_separate[,-c(1:9,omit)]),type='response')
lasso_RHhat_EHR_24_result_separate <- evaluate(CC_val_separate$CC,lasso_RHhat_EHR_pred_24_separate,CC_val_separate$Period)
print(lasso_RHhat_EHR_24_result_separate)



omit <- which(colnames(CC_all) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS','PRIOR_RELAPSE_24MONS','RH_hat'))

lasso_24_RHhat <- cv.glmnet(as.matrix(CC_all[,-c(1:9,omit)]),CC_all$PRIOR_RELAPSE_12MONS>0,family='binomial')
CC_all$RH_hat <- predict(lasso_24_RHhat,as.matrix(CC_all[,-c(1:9,omit)]),type='response')[,1]
CC_val_separate$RH_hat <- predict(lasso_24_RHhat,as.matrix(CC_val_separate[,-c(1:9,omit)]),type='response')[,1]
evaluate(CC_val_separate$PRIOR_RELAPSE_12MONS>0,CC_val_separate$RH_hat,CC_val_separate$Period)
cor(CC_val_separate$PRIOR_RELAPSE_12MONS,CC_val_separate$RH_hat)

omit <- which(colnames(CC_all2) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS','PRIOR_RELAPSE_24MONS','RH_hat'))

lasso_12_RHhat <- cv.glmnet(as.matrix(CC_all2[,-c(1:9,omit)]),CC_all2$PRIOR_RELAPSE_12MONS>0,family='binomial')
CC_all2$RH_hat <- predict(lasso_12_RHhat,as.matrix(CC_all2[,-c(1:9,omit)]),type='response')[,1]
CC_val2_separate$RH_hat <- predict(lasso_12_RHhat,as.matrix(CC_val2_separate[,-c(1:9,omit)]),type='response')[,1]
evaluate(CC_val2_separate$PRIOR_RELAPSE_12MONS>0,CC_val2_separate$RH_hat,CC_val2_separate$Period)
cor(CC_val2_separate$PRIOR_RELAPSE_12MONS,CC_val2_separate$RH_hat)


lasso_RHhat_24 <- cv.glmnet(as.matrix(CC_all[,c('Age','DISEASE_DURA','FEMALE','RACE','RH_hat')]),
                            CC_all$CC,family='binomial',type.measure='auc')
lasso_RHhat_pred_24_separate <- predict(lasso_RHhat_24,as.matrix(CC_val_separate[,c('Age','DISEASE_DURA','FEMALE','RACE','RH_hat')]))
lasso_RHhat_24_result_separate <- evaluate(CC_val_separate$CC,lasso_RHhat_pred_24_separate)
print(lasso_RHhat_24_result_separate)

lasso_RHhat_12 <- cv.glmnet(as.matrix(CC_all2[,c('Age','DISEASE_DURA','FEMALE','RACE','RH_hat')]),
                            CC_all2$CC,family='binomial',type.measure='auc')
lasso_RHhat_pred_12_separate <- predict(lasso_RHhat_12,as.matrix(CC_val2_separate[,c('Age','DISEASE_DURA','FEMALE','RACE','RH_hat')]))
lasso_RHhat_12_result_separate <- evaluate(CC_val2_separate$CC,lasso_RHhat_pred_12_separate)
print(lasso_RHhat_12_result_separate)


omit <- which(colnames(CC_all) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS','PRIOR_RELAPSE_24MONS','RH_hat'))

lasso_24_RHhat <- cv.glmnet(as.matrix(CC_all[,-c(1:9,omit)]),CC_all$PRIOR_RELAPSE_12MONS,family='poisson')
CC_all$RH_hat <- predict(lasso_24_RHhat,as.matrix(CC_all[,-c(1:9,omit)]),type='response')[,1]
CC_val_separate$RH_hat <- predict(lasso_24_RHhat,as.matrix(CC_val_separate[,-c(1:9,omit)]),type='response')[,1]
evaluate(CC_val_separate$PRIOR_RELAPSE_12MONS>0,CC_val_separate$RH_hat,CC_val_separate$Period)
cor(CC_val_separate$PRIOR_RELAPSE_12MONS,CC_val_separate$RH_hat)

omit <- which(colnames(CC_all2) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS','PRIOR_RELAPSE_24MONS','RH_hat'))

lasso_12_RHhat <- cv.glmnet(as.matrix(CC_all2[,-c(1:9,omit)]),CC_all2$PRIOR_RELAPSE_12MONS,family='poisson')
CC_all2$RH_hat <- predict(lasso_12_RHhat,as.matrix(CC_all2[,-c(1:9,omit)]),type='response')[,1]
CC_val2_separate$RH_hat <- predict(lasso_12_RHhat,as.matrix(CC_val2_separate[,-c(1:9,omit)]),type='response')[,1]
evaluate(CC_val2_separate$PRIOR_RELAPSE_12MONS>0,CC_val2_separate$RH_hat,CC_val2_separate$Period)
cor(CC_val2_separate$PRIOR_RELAPSE_12MONS,CC_val2_separate$RH_hat)


lasso_RHhat_24 <- cv.glmnet(as.matrix(CC_all[,c('Age','DISEASE_DURA','FEMALE','RACE','RH_hat')]),
                            CC_all$CC,family='binomial',type.measure='auc')
lasso_RHhat_pred_24_separate <- predict(lasso_RHhat_24,as.matrix(CC_val_separate[,c('Age','DISEASE_DURA','FEMALE','RACE','RH_hat')]))
lasso_RHhat_24_result_separate <- evaluate(CC_val_separate$CC,lasso_RHhat_pred_24_separate)
print(lasso_RHhat_24_result_separate)

lasso_RHhat_12 <- cv.glmnet(as.matrix(CC_all2[,c('Age','DISEASE_DURA','FEMALE','RACE','RH_hat')]),
                            CC_all2$CC,family='binomial',type.measure='auc')
lasso_RHhat_pred_12_separate <- predict(lasso_RHhat_12,as.matrix(CC_val2_separate[,c('Age','DISEASE_DURA','FEMALE','RACE','RH_hat')]))
lasso_RHhat_12_result_separate <- evaluate(CC_val2_separate$CC,lasso_RHhat_pred_12_separate)
print(lasso_RHhat_12_result_separate)

