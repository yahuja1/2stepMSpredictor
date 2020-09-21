## Figure generating scripts for 2-step MS relapse predictor paper
## Author: Yuri Ahuja
## Last updated: Aug, 2020

library(ggplot2)
library(plotROC)
library(hrbrthemes)


# Figure 1

prevalences <- sapply(6:159,function(i){mean(CC_comb_cohort[,i]>0)})
aucs <- sapply(6:159,function(i){auc(CC_comb_cohort$Y,CC_comb_cohort[,i])})
order <- order(aucs,decreasing=TRUE)
X <- round(exp(CC_comb_cohort[,6:159]) - 1)[CC_comb_cohort$ID==138928,order[1:27][-c(13,25)]]
Y <- CC_comb_cohort$Y[CC_comb_cohort$ID==138928]
T <- CC_comb_cohort$T[CC_comb_cohort$ID==138928]
names <- sort(c('PheCode 335 (Multiple sclerosis)','CUI C0038317 (Steroid)','CUI C0028643 (Sensation loss)','CUI C0024485 (MRI)',
                'CUI C0026769 (Multiple sclerosis)','CUI C0701466 (Dexamethasone)','CUI C0025815 (Methylprednisolone)','CUI C0949766 (Physical therapy)',
                'CUI C0032743 (PET Scan)','CUI C1517205 (Flare)','CUI C0277556 (Recurrent disease)','CUI C1304680 (Migraine disorders)','CPTCode C96365 (Intravenous infusion)',
                'CUI C1299586 (Has difficulty doing)','CUI C0016911 (Gadolinium)','CUI C0030193 (Sense of pain)','CUI C0009450 (Communicable disease)',
                'CUI C2242996 (Tingling sensation)','CUI C0412675 (Nuclear MRI brain)','CPTCode C70553 (Imaging of Head & Neck)','CUI C0028738 (Nystagmus)','CUI C0029134 (Optic neuritis)',
                'CUI C0344232 (Blurred vision)','CPT Code C72142 (ED visit)','CUI C0012569 (Double vision)'), decreasing=TRUE)
f1.dataframe <- data.frame('Variable'=unlist(sapply(1:ncol(X),function(j){rep(names[j],sum(X[,j]))})),
                           'X'=unlist(sapply(1:ncol(X),function(j){rep(j,sum(X[,j]))})),
                           'T'=unlist(sapply(1:ncol(X),function(j){unlist(sapply(1:nrow(X),function(i){runif(X[i,j],T[i],T[i]+1)}))})))


rect1 <- data.frame (xmin=0, xmax=29, ymin=-Inf, ymax=Inf)
rect2 <- data.frame (xmin=29, xmax=31, ymin=-Inf, ymax=Inf)
rect3 <- data.frame (xmin=31, xmax=50, ymin=-Inf, ymax=Inf)
rect4 <- data.frame (xmin=50, xmax=53, ymin=-Inf, ymax=Inf)
rect5 <- data.frame (xmin=53, xmax=64, ymin=-Inf, ymax=Inf)


tiff('Figure1a.tiff', width = 1980, height = 880)
p1 <- ggplot(f1.dataframe, aes(x=T, y=factor(Variable,levels=names))) +
  geom_point(size=6, shape=23) +
  labs(x = "Months Elapsed", y = "(A)") +
  theme(axis.title.y=element_blank()) +
  theme(text = element_text(size=32,face="bold")) +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  geom_hline(yintercept=seq(1,25),size=0.75,linetype='dotted') +
  geom_rect(data=rect1, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='blue',alpha=0.1,inherit.aes=FALSE) +
  geom_rect(data=rect2, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='red',alpha=0.1,inherit.aes=FALSE) +
  geom_rect(data=rect3, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='blue',alpha=0.1,inherit.aes=FALSE) +
  geom_rect(data=rect4, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='red',alpha=0.1,inherit.aes=FALSE) +
  geom_rect(data=rect5, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='blue',alpha=0.1,inherit.aes=FALSE)
p1
dev.off()



# Figure 2

load('2step_bootstrap_results.RData')
bootstrap_results_mean <- apply(bootstrap_results,c(1,2),mean)
bootstrap_results_2.5ile <- apply(bootstrap_results,c(1,2),function(x){quantile(x,.025,na.rm=T)})
bootstrap_results_97.5ile <- apply(bootstrap_results,c(1,2),function(x){quantile(x,.975,na.rm=T)})

bootstrap_auc_dataframe2 <- data.frame('mean'=bootstrap_results_mean['AUC',],
                                       '2.5%ile'=bootstrap_results_2.5ile['AUC',],
                                       '97.5%ile'=bootstrap_results_97.5ile['AUC',],
                                       'period'=rep(c(24,12),7),
                                       'covariates'=rep(c('ASRD','ASRD + PheCode','ASRD + RH','ASRD + EHR','ASRD + RH + EHR',
                                                          'ASRD + RH^','ASRD + RH^ + EHR'),each=2),
                                       'alpha'=rep(c(1,1,1,0.75,0.75,1,0.75),each=2))

bootstrap_fscore_dataframe2 <- data.frame('mean'=bootstrap_results_mean[seq(11,length(bootstrap_results_mean),11)],
                                          '2.5%ile'=bootstrap_results_2.5ile[seq(11,length(bootstrap_results_mean),11)],
                                          '97.5%ile'=bootstrap_results_97.5ile[seq(11,length(bootstrap_results_mean),11)],
                                          'period'=rep(c(24,12),7),
                                          'covariates'=rep(c('ASRD','ASRD + PheCode','ASRD + RH','ASRD + EHR','ASRD + RH + EHR',
                                                             'ASRD + RH^','ASRD + RH^ + EHR'),each=2),
                                          'alpha'=rep(c(1,1,1,0.75,0.75,1,0.75),each=2))

bootstrap_auc_dataframe2$type <- c(rep('ASRD',2),rep('Baseline',2),rep('RH',2),rep('Baseline',2),rep('RH',2),rep('RH*',4))
bootstrap_fscore_dataframe2$type <- c(rep('ASRD',2),rep('Baseline',2),rep('RH',2),rep('Baseline',2),rep('RH',2),rep('RH*',4))


limits <- aes(ymax = bootstrap_auc_dataframe2$X2.5.ile[bootstrap_auc_dataframe2$period==12],
              ymin = bootstrap_auc_dataframe2$X97.5.ile[bootstrap_auc_dataframe2$period==12])
p1 = ggplot(bootstrap_auc_dataframe2[bootstrap_auc_dataframe2$period==12,], aes(x = factor(covariates, levels=c('ASRD','ASRD + PheCode','ASRD + EHR','ASRD + RH','ASRD + RH + EHR',
                                                                                                                'ASRD + RH^','ASRD + RH^ + EHR')),
                                                                                y = mean,
                                                                                fill = factor(type),
                                                                                alpha=alpha)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(limits, position=position_dodge(0.9), width=0.1) +
  scale_alpha_continuous(range=c(0.6,1)) +
  labs(y = "AUC", fill="") +
  theme(axis.title.y=element_blank()) +
  #theme(axis.title.x=element_blank()) +
  theme(legend.position="none") +
  theme(text = element_text(size=36,face="bold")) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Spectral") +
  coord_flip(ylim=c(0.66,0.72))
p1 

limits <- aes(ymax = bootstrap_fscore_dataframe2$X2.5.ile[bootstrap_auc_dataframe2$period==12],
              ymin = bootstrap_fscore_dataframe2$X97.5.ile[bootstrap_auc_dataframe2$period==12])
p2 = ggplot(bootstrap_fscore_dataframe2[bootstrap_fscore_dataframe2$period==12,], aes(x = factor(covariates, levels=c('ASRD','ASRD + PheCode','ASRD + EHR','ASRD + RH','ASRD + RH + EHR',
                                                                                                                      'ASRD + RH^','ASRD + RH^ + EHR')),
                                                                                      y = mean,
                                                                                      fill = factor(type),
                                                                                      alpha=alpha)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(limits, position=position_dodge(0.9), width=0.1) +
  scale_alpha_continuous(range=c(0.6,1)) +
  labs(y = "F score", fill="") +
  theme(axis.title.y=element_blank()) +
  #  theme(axis.title.x=element_blank()) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(text = element_text(size=36,face="bold")) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Spectral") +
  coord_flip(ylim=c(0.28,0.345))
p2 

library(gridExtra)
tiff("Figure2.tiff", width = 1980, height = 880)
grid.arrange(p1,p2,ncol=2, widths=c(1.25,1))
dev.off()



# Figure 1S

limits <- aes(ymax = bootstrap_auc_dataframe2$X2.5.ile[bootstrap_auc_dataframe2$period==24],
              ymin = bootstrap_auc_dataframe2$X97.5.ile[bootstrap_auc_dataframe2$period==24])
p1 = ggplot(bootstrap_auc_dataframe2[bootstrap_auc_dataframe2$period==24,], aes(x = factor(covariates, levels=c('ASRD','ASRD + PheCode','ASRD + EHR','ASRD + RH','ASRD + RH + EHR',
                                                                                                                'ASRD + RH^','ASRD + RH^ + EHR')),
                                                                                y = mean,
                                                                                fill = factor(type),
                                                                                alpha=alpha)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(limits, position=position_dodge(0.9), width=0.1) +
  scale_alpha_continuous(range=c(0.6,1)) +
  labs(y = "AUC", fill="") +
  theme(axis.title.y=element_blank()) +
  #theme(axis.title.x=element_blank()) +
  theme(legend.position="none") +
  theme(text = element_text(size=36,face="bold")) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Spectral") +
  coord_flip(ylim=c(0.66,0.725))
p1 

limits <- aes(ymax = bootstrap_fscore_dataframe2$X2.5.ile[bootstrap_auc_dataframe2$period==24],
              ymin = bootstrap_fscore_dataframe2$X97.5.ile[bootstrap_auc_dataframe2$period==24])
p2 = ggplot(bootstrap_fscore_dataframe2[bootstrap_auc_dataframe2$period==24,], aes(x = factor(covariates, levels=c('ASRD','ASRD + PheCode','ASRD + EHR','ASRD + RH','ASRD + RH + EHR',
                                                                                                                   'ASRD + RH^','ASRD + RH^ + EHR')),
                                                                                   y = mean,
                                                                                   fill = factor(type),
                                                                                   alpha=alpha)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(limits, position=position_dodge(0.9), width=0.1) +
  scale_alpha_continuous(range=c(0.6,1)) +
  labs(y = "F score", fill="") +
  theme(axis.title.y=element_blank()) +
  #  theme(axis.title.x=element_blank()) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(text = element_text(size=36,face="bold")) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Spectral") +
  coord_flip(ylim=c(0.39,0.465))
p2 


library(gridExtra)
tiff("Figure1S.tiff", width = 1980, height = 880)
grid.arrange(p1,p2,ncol=2, widths=c(1.25,1))
dev.off()



# Figure 3

dataframe_24_3 <- data.frame('Predicted_dems'=as.vector(predictions_24_3_dems), 'Predicted_phecode'=as.vector(predictions_24_3_phecode),
                             'Predicted_Y'=as.vector(predictions_24_3_Y), 'Predicted_EHR'=as.vector(predictions_24_3_EHR_noY),
                             'Predicted_Y_EHR'=as.vector(predictions_24_3_EHR),'Predicted_FHhat'=as.vector(lasso_RHhat_pred_24_separate),
                             'Actual'=CC_val_separate$CC)
dataframe_12_3 <- data.frame('Predicted_dems'=as.vector(predictions_12_3_dems), 'Predicted_phecode'=as.vector(predictions_12_3_phecode),
                             'Predicted_Y'=as.vector(predictions_12_3_Y), 'Predicted_EHR'=as.vector(predictions_12_3_EHR_noY),
                             'Predicted_Y_EHR'=as.vector(predictions_12_3_EHR),'Predicted_FHhat'=as.vector(prediction_12_3_RHhat),
                             'Actual'=CC_val2_separate$CC)
dataframe_combined <- cbind(dataframe_24_3,dataframe_12_3)

tiff('24_3_ROC_plot.tiff',width=880,height=880)
roc_24_3 <- ggplot(dataframe_24_3) + 
  geom_roc(aes(d=Actual, m=Predicted_dems, color='ASRD'), n.cuts=0) +
  geom_roc(aes(d=Actual, m=Predicted_phecode, color='ASRD + Phecode'), n.cuts=0) +
  geom_roc(aes(d=Actual, m=Predicted_Y, color='ASRD + RH'), n.cuts=0) +
  geom_roc(aes(d=Actual, m=Predicted_EHR, color='ASRD + EHR'), n.cuts=0) +
  geom_roc(aes(d=Actual, m=Predicted_Y_EHR, color='ASRD + RH + EHR'), n.cuts=0) +
  geom_roc(aes(d=Actual, m=Predicted_FHhat, color='ASRD + RH^'), n.cuts=0) +
  geom_abline(intercept = 0, slope = 1, color = "darkgrey", linetype = "dashed") +
  labs(x='False Positive Rate', y='True Positive Rate') +
  theme(legend.title=element_blank()) +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  theme(text = element_text(size=36)) +
  scale_color_brewer(palette = 'Set1')
roc_24_3
dev.off()

tiff('Figure3.tiff',width=1320,height=880)
roc_12_3 <- ggplot(dataframe_12_3) + 
  geom_roc(aes(d=Actual, m=Predicted_dems, color='ASRD'), n.cuts=0) +
#  geom_roc(aes(d=Actual, m=Predicted_phecode, color='ASRD + Phecode'), n.cuts=0) +
  geom_roc(aes(d=Actual, m=Predicted_Y, color='ASRD + RH'), n.cuts=0) +
#  geom_roc(aes(d=Actual, m=Predicted_EHR, color='ASRD + EHR'), n.cuts=0) +
#  geom_roc(aes(d=Actual, m=Predicted_Y_EHR, color='ASRD + RH + EHR'), n.cuts=0) +
  geom_roc(aes(d=Actual, m=Predicted_FHhat, color='ASRD + RH^'), n.cuts=0) +
  geom_abline(intercept = 0, slope = 1, color = "darkgrey", linetype = "dashed") +
  labs(x='False Positive Rate', y='True Positive Rate') +
  theme(legend.title=element_blank()) +
  theme(text = element_text(size=36)) +
  scale_color_brewer(palette = 'Set1')
roc_12_3
dev.off()

library(gridExtra)
tiff("Figure3.tiff", width = 1980, height = 880)
grid.arrange(roc_12_3,roc_24_3,ncol=2, widths=c(1,1.25))
dev.off()



# Figure 4

library(reshape2)
library(ggplot2)

X <- CC_all[,-c(1:9)]
lasso <- lasso_EHR_24
X_lasso <- X[,as.vector(coef(lasso))[-1]!=0]
corXlasso <- cor(X_lasso)

get_upper_tri <- function(dat){
  dat[lower.tri(dat)] <- NA
  return(dat)
}

switch <- function(mat,col1,col2){
  temp <- mat[,col1]
  mat[,col1] <- mat[,col2]
  mat[,col2] <- temp
  temp <- colnames(mat)[col1]
  colnames(mat)[col1] <- colnames(mat)[col2]
  colnames(mat)[col2] <- temp
  return(mat)
}

X_lasso_ordered <- switch(X_lasso,1,which(colnames(X_lasso)=='PheCode.335_'))
for (i in 2:(ncol(X_lasso_ordered)-1)){
  newcol <- (i-1) + which.max(sapply(i:ncol(X_lasso_ordered), function(j){cor(X_lasso_ordered[,(i-1)],X_lasso_ordered[,j])}))
  X_lasso_ordered <- switch(X_lasso_ordered,i,newcol)
}
corXlasso <- cor(X_lasso_ordered)
corXlassoNames <- sapply(colnames(corXlasso),function(col){
  truncated <- strsplit(col,' ')[[1]]
  truncated <- truncated[1:min(2,length(truncated))]
  paste(truncated, collapse=' ')})
names(corXlassoNames) <- c()
corXlassoNames <- sapply(colnames(corXlasso),function(col){substr(col,1,25)})

dataframe_lasso_heatplot <- data.frame(X=as.factor(rep(corXlassoNames,each=length(corXlassoNames))),
                                       Y=as.factor(rep(corXlassoNames,length(corXlassoNames))),
                                       Z=c(get_upper_tri(corXlasso)))

dataframe_lasso_heatplot$heat <- sapply(dataframe_lasso_heatplot$Z,function(x){
  if (is.na(x)){'0'}
  else if (x < -0.1){'< -0.1'}
  else if (x < -0.01){'-0.1 - 0.0'}
  else if (x < 0.01){'0'}
  else if (x < 0.1){'0.0 - 0.1'}
  else if (x < 0.2){'0.1 - 0.2'}
  else if (x < 0.3){'0.2 - 0.3'}
  else if (x < 0.4){'0.3 - 0.4'}
  else{'> 0.4'}
})

tiff('heatmap_lasso_24_3.tiff', width=1400, height=1200)
ggplot(dataframe_lasso_heatplot, aes(factor(X,levels=unique(X)),factor(Y,levels=unique(Y)),
                                     fill=factor(heat,levels=c('< -0.1','-0.1 - 0.0','0','0.0 - 0.1','0.1 - 0.2','0.2 - 0.3','0.3 - 0.4','> 0.4')))) + 
  geom_tile() + 
  labs(x='',y='',fill='') + 
  scale_fill_manual(values=c('indianred3','indianred1','white','dodgerblue1','dodgerblue2','dodgerblue3','dodgerblue4','darkblue')) +
  theme(axis.text.x=element_text(angle=90,hjust=1,size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  theme(legend.text=element_text(size=30)) +
  theme(legend.position = 'left') +
  scale_y_discrete(position='right')
dev.off()

omit <- which(colnames(CC_all2) %in% c('FOLLOWUP_DURA','PRIOR_RELAPSE_12MONS','PRIOR_RELAPSE_24MONS','RH_hat'))

X <- CC_all2[,-c(1:9,omit)]
lasso2 <- lasso_12_RHhat
X_lasso <- X[,as.vector(coef(lasso2))[-1]!=0]
corXlasso <- cor(X_lasso)

switch <- function(mat,col1,col2){
  temp <- mat[,col1]
  mat[,col1] <- mat[,col2]
  mat[,col2] <- temp
  temp <- colnames(mat)[col1]
  colnames(mat)[col1] <- colnames(mat)[col2]
  colnames(mat)[col2] <- temp
  return(mat)
}

X_lasso_ordered <- switch(X_lasso,1,which(colnames(X_lasso)=='CUI.C0277556'))
for (i in 2:(ncol(X_lasso_ordered)-1)){
  newcol <- (i-1) + which.max(sapply(i:ncol(X_lasso_ordered), function(j){cor(X_lasso_ordered[,(i-1)],X_lasso_ordered[,j])}))
  X_lasso_ordered <- switch(X_lasso_ordered,i,newcol)
}
corXlasso <- cor(X_lasso_ordered)
corXlassoNames <- sapply(colnames(corXlasso),function(col){
  truncated <- strsplit(col,'[.]')[[1]]
  truncated <- truncated[1:min(3,length(truncated))]
  paste(truncated, collapse=' ')})
names(corXlassoNames) <- c()
corXlassoNames <- sapply(colnames(corXlasso),function(col){substr(col,1,25)})

dataframe_lasso_heatplot <- data.frame(X=as.factor(rep(corXlassoNames,each=length(corXlassoNames))),
                                       Y=as.factor(rep(corXlassoNames,length(corXlassoNames))),
                                       Z=c(get_upper_tri(corXlasso)))

dataframe_lasso_heatplot$heat <- sapply(dataframe_lasso_heatplot$Z,function(x){
  if (is.na(x)){'0'}
  else if (x < -0.1){'< -0.1'}
  else if (x < -0.01){'-0.1 - 0.0'}
  else if (x < 0.01){'0'}
  else if (x < 0.1){'0.0 - 0.1'}
  else if (x < 0.2){'0.1 - 0.2'}
  else if (x < 0.3){'0.2 - 0.3'}
  else if (x < 0.4){'0.3 - 0.4'}
  else{'> 0.4'}
})

tiff('Figure4.tiff', width=1400, height=1200)
ggplot(dataframe_lasso_heatplot, aes(factor(X,levels=unique(X)),factor(Y,levels=unique(Y)),
                                     fill=factor(heat,levels=c('< -0.1','-0.1 - 0.0','0','0.0 - 0.1','0.1 - 0.2','0.2 - 0.3','0.3 - 0.4','> 0.4')))) + 
  geom_tile() + 
  labs(x='',y='',fill='') + 
  scale_fill_manual(values=c('indianred3','indianred1','white','dodgerblue1','dodgerblue2','dodgerblue3','dodgerblue4','darkblue')) +
  theme(axis.text.x=element_text(angle=90,hjust=1,size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  theme(legend.text=element_text(size=30)) +
  theme(legend.position = 'left') +
  scale_y_discrete(position='right')
dev.off()



# Figure 5
library(gridExtra)
library(grid)

trainPats <- sample(unique(CC_all$PatientNum),0.8*length(unique(CC_all$PatientNum)))
valPats <- setdiff(unique(CC_all$PatientNum),trainPats)
CC_comb <- CC_all[CC_all$PatientNum %in% trainPats,]
CC_val <- CC_all[CC_all$PatientNum %in% valPats,]

trainPats2 <- sample(unique(CC_all2$PatientNum),0.8*length(unique(CC_all2$PatientNum)))
valPats2 <- setdiff(unique(CC_all2$PatientNum),trainPats2)
CC_comb2 <- CC_all2[CC_all2$PatientNum %in% trainPats2,]
CC_val2 <- CC_all2[CC_all2$PatientNum %in% valPats2,]

omit <- which(colnames(CC_comb) %in% c('CC','Age','DISEASE_DURA','FOLLOWUP_DURA','PRIOR_RELAPSE_12'))
lasso_notime_24_3 <- cv.glmnet(as.matrix(CC_comb[,c('Age','FEMALE','RACE','DISEASE_DURA','RH_hat')]),CC_comb$CC)
predictions_notime_24_3 <- predict(lasso_notime_24_3,as.matrix(CC_val[,c('Age','FEMALE','RACE','DISEASE_DURA','RH_hat')]),type='response')

omit2 <- which(colnames(CC_comb2) %in% c('CC','Age','DISEASE_DURA','FOLLOWUP_DURA','PRIOR_RELAPSE_12'))
lasso_notime_12_3 <- cv.glmnet(as.matrix(CC_comb2[,c('Age','FEMALE','RACE','DISEASE_DURA','RH_hat')]),CC_comb2$CC)
predictions_notime_12_3 <- predict(lasso_notime_12_3,as.matrix(CC_val2[,c('Age','FEMALE','RACE','DISEASE_DURA','RH_hat')]),type='response')

DD <- seq(1, ceiling(max(CC_val$DISEASE_DURA)), 0.5)
RR_DD <- sapply(DD,function(d){
  matches <- which(abs(CC_val$DISEASE_DURA-d) <= 1)
  if (length(matches)>0){
    c(mean(CC_val$CC[matches]), sd(CC_val$CC[matches]), length(matches))
  }
  else{
    rep(NA,3)
  }
})
PR_DD <- sapply(DD,function(d){
  matches <- which(abs(CC_val$DISEASE_DURA-d) <= 1)
  if (length(matches)>0){
    c(mean(predictions_notime_24_3[matches]), sd(predictions_notime_24_3[matches]), length(matches))
  }
  else{
    rep(NA,3)
  }
})

df <- data.frame(DD=DD,RR_mean=RR_DD[1,],RR_lower=sapply(1:ncol(RR_DD),function(i){max(0,RR_DD[1,i]-1.96*RR_DD[2,i]/sqrt(RR_DD[3,i]))}),RR_upper=RR_DD[1,]+1.96*RR_DD[2,]/sqrt(RR_DD[3,]),count=RR_DD[3,],
                 PR_mean=PR_DD[1,],PR_lower=sapply(1:ncol(PR_Age),function(i){max(0,PR_Age[1,i]-PR_Age[2,i])}),PR_upper=PR_DD[1,]+PR_DD[2,],count=PR_DD[3,])
df <- df[df$count>=50 & df$DD <= 40,]
tiff('Disease_Duration_samePlot_24mo_lasso.tiff', width=1320, height=660)
p1 <- ggplot(df) +
  geom_line(aes(x = DD, y = RR_mean, color='Actual')) +
  geom_point(aes(x = DD, y = RR_mean, color='Actual'),size=2) +
  geom_ribbon(aes(x=DD,ymin=RR_lower,ymax=RR_upper,fill='Actual'),alpha=0.1) +
  geom_line(aes(x = DD, y = PR_mean,color='Predicted')) +
  geom_point(aes(x = DD, y = PR_mean,color='Predicted'),size=2) +
  geom_ribbon(aes(x=DD,ymin=PR_lower,ymax=PR_upper,fill='Predicted'),alpha=0.1) +
  theme(plot.margin=unit(c(27.5,5.5,5.5,5.5),'points')) +
  xlab('Disease Duration (Years)') +
  theme(legend.title=element_blank(), axis.title.y=element_blank()) +
  theme(text = element_text(size=36)) +
  scale_fill_brewer(palette = 'Set1') +
  scale_color_brewer(palette = 'Set1') +
  guides(fill=FALSE)
for (i in 1:length(df$count)){
  p1 <- p1 + annotation_custom(
    grob = textGrob(label = df$count[i], hjust = 0, gp = gpar(cex = 1), rot=45),
    ymin = 2.61,      # Vertical position of the textGrob
    ymax = 2.61,
    xmin = df$DD[i],         # Note: The grobs are positioned outside the plot area
    xmax = df$DD[i])
}
g1 <- ggplot_gtable(ggplot_build(p1))
g1$layout$clip[g1$layout$name == "panel"] <- "off"
grid.draw(g1)
dev.off()

Age <- seq(floor(min(CC_val_separate$Age)), floor(max(CC_val_separate$Age)))
RR_Age <- sapply(Age,function(a){
  matches <- which(abs(CC_val$Age-a) <= 2.5)
  if (length(matches)>0){
    c(mean(CC_val$CC[matches]), sd(CC_val$CC[matches]), length(matches))
  }
  else{
    rep(NA,3)
  }
})
PR_Age <- sapply(Age,function(a){
  matches <- which(abs(CC_val$Age-a) <= 2.5)
  if (length(matches)>0){
    c(mean(predictions_notime_24_3[matches]), sd(predictions_notime_24_3[matches]), length(matches))
  }
  else{
    rep(NA,3)
  }
})

df2 <- data.frame(Age=Age,RR_mean=RR_Age[1,],RR_lower=sapply(1:ncol(RR_Age),function(i){max(0,RR_Age[1,i]-1.96*RR_Age[2,i]/sqrt(RR_Age[3,i]))}),RR_upper=RR_Age[1,]+1.96*RR_Age[2,]/sqrt(RR_Age[3,]),count=RR_Age[3,],
                  PR_mean=PR_Age[1,],PR_lower=sapply(1:ncol(PR_Age),function(i){max(0,PR_Age[1,i]-PR_Age[2,i])}),PR_upper=PR_Age[1,]+PR_Age[2,],count=PR_Age[3,])
df2 <- df2[df2$count>=100,]
tiff('Age_samePlot_24mo_lasso.tiff', width=1320, height=660)
p1 <- ggplot(df2) +
  geom_line(aes(x = Age, y = RR_mean, color='Actual')) +
  geom_point(aes(x = Age, y = RR_mean, color='Actual'),size=2) +
  geom_ribbon(aes(x=Age,ymin=RR_lower,ymax=RR_upper,fill='Actual'),alpha=0.1) +
  geom_line(aes(x = Age, y = PR_mean,color='Predicted')) +
  geom_point(aes(x = Age, y = PR_mean,color='Predicted'),size=2) +
  geom_ribbon(aes(x=Age,ymin=PR_lower,ymax=PR_upper,fill='Predicted'),alpha=0.1) +
  theme(plot.margin=unit(c(27.5,5.5,5.5,5.5),'points')) +
  xlab('Age (Years)') +
  theme(legend.title=element_blank(), axis.title.y=element_blank()) +
  theme(text = element_text(size=36)) +
  scale_fill_brewer(palette = 'Set1') +
  scale_color_brewer(palette = 'Set1') +
  guides(fill=FALSE)
for (i in 1:length(df$count)){
  p1 <- p1 + annotation_custom(
    grob = textGrob(label = df$count[i], hjust = 0, gp = gpar(cex = 1), rot=45),
    ymin = 2.61,      # Vertical position of the textGrob
    ymax = 2.61,
    xmin = df$DD[i],         # Note: The grobs are positioned outside the plot area
    xmax = df$DD[i])
}
g1 <- ggplot_gtable(ggplot_build(p1))
g1$layout$clip[g1$layout$name == "panel"] <- "off"
grid.draw(g1)
dev.off()

RR_DD <- sapply(DD,function(d){
  matches <- which(abs(CC_val2$DISEASE_DURA-d) <= 1)
  if (length(matches)>0){
    c(mean(CC_val2$CC[matches]), sd(CC_val2$CC[matches]), length(matches))
  }
  else{
    rep(NA,3)
  }
})
PR_DD <- sapply(DD,function(d){
  matches <- which(abs(CC_val2$DISEASE_DURA-d) <= 1)
  if (length(matches)>0){
    c(mean(predictions_notime_12_3[matches]), sd(predictions_notime_12_3[matches]), length(matches))
  }
  else{
    rep(NA,3)
  }
})
df <- data.frame(DD=DD,RR_mean=RR_DD[1,],RR_lower=sapply(1:ncol(RR_DD),function(i){max(0,RR_DD[1,i]-1.96*RR_DD[2,i]/sqrt(RR_DD[3,i]))}),RR_upper=RR_DD[1,]+1.96*RR_DD[2,]/sqrt(RR_DD[3,]),count=RR_DD[3,],
                 PR_mean=PR_DD[1,],PR_lower=PR_DD[1,]-PR_DD[2,],PR_upper=PR_DD[1,]+PR_DD[2,],count=PR_DD[3,])
df <- df[df$count>=50 & df$DD <= 40,]
tiff('Disease_Duration_samePlot_12mo_lasso.tiff', width=1320, height=660)
p1 <- ggplot(df) +
  geom_line(aes(x = DD, y = RR_mean, color='red')) +
  geom_point(aes(x = DD, y = RR_mean, color='red'),size=2) +
  geom_ribbon(aes(x=DD,ymin=RR_lower,ymax=RR_upper,fill='red'),alpha=0.1) +
  geom_line(aes(x = DD, y = PR_mean,color='darkblue')) +
  geom_point(aes(x = DD, y = PR_mean,color='darkblue'),size=2) +
  geom_ribbon(aes(x=DD,ymin=PR_lower,ymax=PR_upper,fill='darkblue'),alpha=0.1) +
  theme(plot.margin=unit(c(27.5,5.5,5.5,5.5),'points')) +
  scale_color_manual(values=c('blue','red')) +
  scale_fill_manual(values=c('blue','red')) +
  xlab('Disease Duration (Years)') +
  theme(legend.title=element_blank()) +
  ylab("Actual Relapse Proportion") +
  scale_y_continuous(limits = c(0.0,0.20)) +
  theme(text = element_text(size=36)) +
  theme(legend.position="none") +
  theme(axis.title.y = element_text(color='red'), axis.text.y = element_text(color='red')) +
  guides(fill=FALSE)
for (i in 1:length(df$count)){
  p1 <- p1 + annotation_custom(
    grob = textGrob(label = df$count[i], hjust = 0, gp = gpar(cex = 1), rot=45),
    ymin = 2.61,      # Vertical position of the textGrob
    ymax = 2.61,
    xmin = df$DD[i],         # Note: The grobs are positioned outside the plot area
    xmax = df$DD[i])
}
g1 <- ggplot_gtable(ggplot_build(p1))
g1$layout$clip[g1$layout$name == "panel"] <- "off"
grid.draw(g1)
dev.off()

RR_Age <- sapply(Age,function(a){
  matches <- which(round(CC_val2$Age-a) <= 2.5)
  if (length(matches)>0){
    c(mean(CC_val2$CC[matches]), sd(CC_val2$CC[matches]), length(matches))
  }
  else{
    rep(NA,3)
  }
})
PR_Age <- sapply(Age,function(a){
  matches <- which(round(CC_val2$Age-a) <= 2.5)
  if (length(matches)>0){
    c(mean(predictions_notime_12_3[matches]), sd(predictions_notime_12_3[matches]), length(matches))
  }
  else{
    rep(NA,3)
  }
})
df2 <- data.frame(Age=Age,RR_mean=RR_Age[1,],RR_lower=sapply(1:ncol(RR_Age),function(i){max(0,RR_Age[1,i]-1.96*RR_Age[2,i]/sqrt(RR_Age[3,i]))}),RR_upper=RR_Age[1,]+1.96*RR_Age[2,]/sqrt(RR_Age[3,]),count=RR_Age[3,],
                  PR_mean=PR_Age[1,],PR_lower=sapply(1:ncol(PR_Age),function(i){max(0,PR_Age[1,i]-PR_Age[2,i])}),PR_upper=PR_Age[1,]+PR_Age[2,],count=PR_Age[3,])
df2 <- df2[df2$count>=100 & df2$Age>=30,]
tiff('Age_samePlot_12mo_lasso.tiff', width=1320, height=660)
p2 <- ggplot(df2) +
  geom_line(aes(x = Age, y = RR_mean, color='red')) +
  geom_point(aes(x = Age, y = RR_mean, color='red'),size=2) +
  geom_ribbon(aes(x=Age,ymin=RR_lower,ymax=RR_upper,fill='red'),alpha=0.1) +
  geom_line(aes(x = Age, y = PR_mean,color='darkblue')) +
  geom_point(aes(x = Age, y = PR_mean,color='darkblue'),size=2) +
  geom_ribbon(aes(x=Age,ymin=PR_lower,ymax=PR_upper,fill='darkblue'),alpha=0.1) +
  theme(plot.margin=unit(c(27.5,5.5,5.5,5.5),'points')) +
  scale_color_manual(values=c('blue','red')) +
  scale_fill_manual(values=c('blue','red')) +
  xlab('Age (Years)') +
  theme(legend.title=element_blank()) +
  ylab("Predicted Relapse Probability") +
  theme(text = element_text(size=36)) +
  theme(legend.position="none") +
  theme(axis.title.y = element_text(color='blue'), axis.text.y = element_text(color='blue')) +
  scale_y_continuous(limits = c(0.0,0.20), position='right') +
  guides(fill=FALSE)
for (i in 1:length(df$count)){
  p2 <- p2 + annotation_custom(
    grob = textGrob(label = df$count[i], hjust = 0, gp = gpar(cex = 1), rot=45),
    ymin = 2.61,      # Vertical position of the textGrob
    ymax = 2.61,
    xmin = df$DD[i],         # Note: The grobs are positioned outside the plot area
    xmax = df$DD[i])
}
g1 <- ggplot_gtable(ggplot_build(p2))
g1$layout$clip[g1$layout$name == "panel"] <- "off"
grid.draw(g1)
dev.off()

library(gridExtra)
tiff("Figure5.tiff", width = 1980, height = 880)
grid.arrange(p1,p2,ncol=2, widths=c(1,1))
dev.off()



# Table 2

cui_descriptors <- sapply(colnames(X_lasso)[19:77],function(name){
  print(name)
  name <- strsplit(name,'[.]')[[1]][2]
  matches <- which(dictionary[,4]==name)
  dictionary[min(matches),1]
})
write.csv(cui_descriptors,'Table2.csv')



# Table 1S

glm_24_DD_Actual <- glm(CC~DISEASE_DURA,data=CC_val_separate,family='quasibinomial')
summary(glm_24_DD_Actual)
cor(CC_val_separate$CC,CC_val_separate$DISEASE_DURA,method='spearman')

glm_24_DD_Predicted <- glm(predictions_notime_24_3~CC_val_separate$DISEASE_DURA,family='quasibinomial')
summary(glm_24_DD_Predicted)
cor(predictions_notime_24_3,CC_val_separate$DISEASE_DURA,method='spearman')

glm_12_DD_Actual <- glm(CC~DISEASE_DURA,data=CC_val2_separate,family='quasibinomial')
summary(glm_12_DD_Actual)
cor(CC_val2_separate$CC,CC_val2_separate$DISEASE_DURA,method='spearman')

glm_12_DD_Predicted <- glm(predictions_notime_12_3~CC_val2_separate$DISEASE_DURA,family='quasibinomial')
summary(glm_12_DD_Predicted)
cor(predictions_notime_12_3,CC_val2_separate$DISEASE_DURA,method='spearman')


glm_24_Age_Actual <- glm(CC~Age,data=CC_val_separate,family='quasibinomial')
summary(glm_24_Age_Actual)
cor(CC_val_separate$CC,CC_val_separate$Age,method='spearman')

glm_24_Age_Predicted <- glm(predictions_notime_24_3~CC_val_separate$Age,family='quasibinomial')
summary(glm_24_Age_Predicted)
cor(predictions_notime_24_3,CC_val_separate$Age,method='spearman')

glm_12_Age_Actual <- glm(CC~Age,data=CC_val2_separate,family='quasibinomial')
summary(glm_12_Age_Actual)
cor(CC_val2_separate$CC,CC_val2_separate$Age,method='spearman')

glm_12_Age_Predicted <- glm(predictions_notime_12_3~CC_val2_separate$Age,family='quasibinomial')
summary(glm_12_Age_Predicted)
cor(predictions_notime_12_3,CC_val2_separate$Age,method='spearman')



# Deprecated Figure

#plotActualPredictions <- function(id,dat,predictions){
#  matches <- which(dat$PatientNum==id)
#  dataframe_plot <- data.frame('T'=dat$DISEASE_DURA[matches], 'Actual'=dat$CC[matches],
#                               'Predicted'=predictions[matches])
#  p <- ggplot(dataframe_plot, aes(x=T,y=Actual)) + 
#    geom_point(size=2) + 
#    geom_line(aes(x=T,y=Predicted), color='red') + 
#    labs(x='Disease Duration (Years)', y='Relapse') + 
#    theme(text = element_text(size=36))
#  p
#}

#diffs <- t(sapply(unique(CC_comb_val$PatientNum),function(id){
#  matches <- which(CC_comb_val$PatientNum==id)
#  c(id,ifelse(length(unique(CC_comb_val$CC[matches]))==1,NA,auc(CC_comb_val$CC[matches],predictions_24_3[matches])),
#    max(abs(predictions_24_3[matches]-CC_comb_val$CC[matches])),
#    median(abs(predictions_24_3[matches]-CC_comb_val$CC[matches])))
#}))
#diffs <- diffs[!is.na(diffs[,2]),]
#print(diffs[order(diffs[,3]),][1:20,])


#ids <- c(94632,118903,119714,123698,138928,157163)
#ids <- c(118903,138928)
#for (id in ids){
#  matches <- which(CC_comb_val$PatientNum==id)
#  plot(CC_comb_val$CC[matches],type='p',col='black',xlab='Timepoint',ylab='Relapse Status')
#  lines(predictions_24_3[matches],type='l',col='red')
  
#  tiff(paste0('Actual_vs_Predicted_',id,'_24_3.tiff'),width=990,height=660)
#  plotActualPredictions(id,CC_comb_val,predictions_24_3)
#  dev.off()
#}

#ids <- c(94632,118903,119714,123698,138928,157163)
#ids <- c(118903,138928)
#for (id in ids){
#  matches <- which(CC_val2$PatientNum==id)
#  plot(CC_val2$CC[matches],type='p',col='black',xlab='Timepoint',ylab='Relapse Status')
#  lines(predictions_12_3[matches],type='l',col='red')
#  lines(predictions_12_3[matches],type='l',col='blue')
  
#  tiff(paste0('Actual_vs_Predicted_',id,'_12_3.tiff'),width=990,height=660)
#  plotActualPredictions(id,CC_comb_val,predictions_12_3)
#  dev.off()
#}

