Est.ALASSO.GLM.OFFSET = function(data,offset,Wi=NULL,rtn="EST",nopen.ind=NULL,regularize=yes.regularize, pwr=pwr){
  data = as.matrix(data); y = data[,1]; x = data[,-1,drop=F]; nn=length(y); if(is.null(Wi)){Wi=rep(1,nn)}; pp = ncol(x)
  if(regularize){
    ## adaptive lasso (aLASSO) with initial estimator obtained via ridge (bini); w.b creates adaptive weights for aLASSO ##
    lam.ridge = pp/nn; ##bini = plr(x,y,lambda=lam.ridge,weights=Wi)$coef; 
    bini = as.vector(coef(glmnet(x,y, offset=offset, weights=Wi,alpha=0,standardize=F,lambda=lam.ridge,family="binomial"))); ##print(bini)
    w.b = 1/abs(bini[-1]); x.t = x/VTM(w.b,nrow(x))
    
    ## glmpath provides solution path for a range of penalty parameters ##
    tmpfit = glmpath(x.t,y,offset=offset, nopenalty.subset=nopen.ind,family=binomial,weight=Wi,standardize=F,min.lambda=0,lambda2=lam.ridge)
    lam.all = c(seq(min(tmpfit$lambda),max(tmpfit$lambda),length=500))
    b.all = predict(tmpfit, s=lam.all, offset=offset, type="coefficients",mode="lambda")
    b.all = b.all/VTM(c(1,w.b),nrow(b.all)); m0 = length(lam.all)
    ## ================================================================================ ##
    ## calculates degree of freedom for all beta's (corresponding to different lam.all) ##
    ## ================================================================================ ##
    df.all = apply(b.all[,-1,drop=F]!=0,1,sum); x = as.matrix(x)
    ## =============================================================================================================== ##
    ## calculates modified BIC, log(n) is modified as min{sum(y)^0.1, log(n)} to avoid over shrinkage in finite sample ##
    ## =============================================================================================================== ##
    BIC.lam = -2*apply(predict(tmpfit,offset=offset, newx=x.t,newy=y,s=lam.all,type="loglik",mode="lambda"),2,sum)+min(sum(y)^pwr,log(sum(y)))*df.all 
    #BIC.lam = -2*apply(predict(tmpfit,offset=offset, newx=x.t,newy=y,s=lam.all,type="loglik",mode="lambda"),2,sum)+log(sum(y))*df.all 
    
    m.opt = (1:m0)[BIC.lam==min(BIC.lam)]; bhat = b.all[m.opt,]; lamhat = lam.all[m.opt]
  }else{
    ## ========================================================================= ##
    ## if regularize = F, then use standard logistic regression w/o penalization ##
    ## ========================================================================= ##
    bhat=bini=glm(y~x,family=binomial,weight=Wi)$coef;lamhat = 0; lam.all=BIC.lam=b.all=NULL
  }
  out = c("b"=bhat, "bini"=bini,"lamhat"=lamhat,"lam.all"=lam.all,"BIC.lam"=BIC.lam,"b.all"=b.all)
  if(rtn=="EST"){return(out)}else{return(list(out,"b.all"=b.all,"lam.all"=lam.all,"fit"=tmpfit,"BIC.lam"=BIC.lam))}
}


ROC.FUN.ALASSO.OFFSET.632boot  <- function(data, offset, ind=NULL, wgti=NULL, Wi=NULL, yy0 = 0.5, nopen.ind=NULL,
                                    FPR0=seq(.01,.99,by=.01),rtn="ALL",yes.CV=F,yes.seed=F,rep=10,regularize=T){
  ## ========================================================================= ##
  ## ==== ind: for bootstrap methods; data[ind,] returns bootstrap samples === ##
  ## ==== 1st column is binary disease status;                            ==== ##
  ## ==== 2nd column and on are individual markes                         ==== ##
  ## ==== yy0: prespecified cut-off;  FPR0: prespecified Specificity level ==== ##
  ## ========================================================================= ##
  nn <- nrow(data); if(is.null(wgti)){wgti=rep(1,nn)}; 
  if(!is.null(ind)){data <- data[ind,]; wgti = wgti[ind]; offset=offset[ind]}    
  n.set <- 1; pp = ncol(data);  yyi.vmat = matrix(NA, nrow=nn, ncol=rep)
  ## === Apparent Accuracy === ##
  beta.offset = try(Est.ALASSO.GLM.OFFSET(data,offset, Wi=wgti, rtn="EST",regularize=regularize,nopen.ind=nopen.ind),silent=T)[1:(pp)]
  betahat=c(beta.offset[1], beta0[-1]+beta.offset[-1])
  if(length(betahat)>1){
    yyi = g.logit(cbind(1,as.matrix(data[,-1]))%*%betahat)    
    q.nonzero = sum(abs(betahat[2:(pp)])!=0)
    if(q.nonzero==0){rochat = NULL}else{rochat = ROC.Est.FUN(data[,1],yyi,yy0,FPR0,wgti=wgti)}
    if(yes.seed){set.seed(315)}
    if(yes.CV&q.nonzero>0){
      ## === 0.632 bootstrap === ##
      roc.cv = NULL; 
      for(i in 1:rep){
        tmpind=sample(1:nn,replace=T); ind.v = setdiff(1:nn,unique(tmpind))
        dat.t = data[tmpind,]; wgti.t = wgti[tmpind]; offset.t=offset[tmpind]; dat.v = data[ind.v,];  wgti.v = wgti[ind.v]; offset.v=offset[ind.v]
        beta.t.offset = try(Est.ALASSO.GLM.OFFSET(data,offset, Wi=wgti, rtn="EST",regularize=regularize,nopen.ind=nopen.ind),silent=T)[1:(pp)]
        beta.t=c(beta.t.offset[1], beta0[-1]+beta.t.offset[-1])
        if(length(beta.t)>1){
          beta.t = beta.t; yyi.v = g.logit(cbind(1,as.matrix(dat.v[,-1]))%*%beta.t); yyi.vmat[ind.v,i]=yyi.v
          roc.k = try(ROC.Est.FUN(dat.v[,1],yyi.v,yy0,FPR0,wgti.v),silent=T)
          if(length(roc.k)>1){roc.cv = cbind(roc.cv, roc.k)}  } }
      roc.cv = apply(roc.cv,1,mean)*0.632+rochat*0.368  
    }else{roc.cv = NULL}
  }else{roc.cv=beta.hat=NULL}
  if(rtn=="EST"){
    if(yes.CV){out=c(roc.cv,betahat)}else{out=c(rochat,betahat)}
    return(out)
  }else{return(list("rochat"=rochat, "roc.cv"=roc.cv, beta.offset=beta.offset, "beta"=betahat, "Si"=yyi, "yyi.v"=yyi.vmat))}
}

RA.sampling.prob.FUN=function(dat){
msrandom = which(dat$Enriched450Review_samplingcohort=="MS")
bdrandom = which(dat$Enriched450Review_samplingcohort=="BD")
rarandom = which(dat$Enriched450Review_samplingcohort=="RA")
ucrandom = which(dat$Enriched450Review_samplingcohort=="UC")
cdrandom = which(dat$Enriched450Review_samplingcohort=="CD")

mspool = which(dat$MS_MeetsFilterCriteria=="Y")
sampled = msrandom
bdpool = setdiff(which(dat$BD_MeetsFilterCriteria=="Y"),sampled)
sampled = c(sampled,bdrandom)
rapool = setdiff(which(dat$RA_MeetsFilterCriteria=="Y"),sampled)
sampled = c(sampled,rarandom)
ucpool = setdiff(which(dat$UC_MeetsFilterCriteria=="Y"),sampled)
sampled = c(sampled,ucrandom)
cdpool = setdiff(which(dat$CD_MeetsFilterCriteria=="Y"),sampled)
sampled = c(sampled,cdrandom)

p.ms = 100/length(mspool)
p.bd = 100/length(bdpool)
p.ra = 100/length(rapool)
p.uc = 75/length(ucpool)
p.cd = 75/length(cdpool)
p.rand = 100/nrow(dat)
p.t2dm = 100/sum(dat$T2DM_MeetsFilterCriteria=='Y')
dat = data.frame(dat,"pi"=1-(1-p.rand)*(1-p.ms*(dat$MS_MeetsFilterCriteria=="Y"))*(1-p.bd*(dat$BD_MeetsFilterCriteria=="Y"))*(1-p.ra*(dat$RA_MeetsFilterCriteria=="Y"))*(1-p.uc*(dat$UC_MeetsFilterCriteria=="Y"))*(1-p.cd*(dat$CD_MeetsFilterCriteria=="Y"))*(1-p.t2dm*(dat$T2DM_MeetsFilterCriteria=='Y')))
return(dat)
}

DAT.FILTER.PREP.FUN=function(dat, dis, patient_num,
                      nm.label,nm.common,nm.filter,nm.cov,nm.notrans, standardize=F){
  rownames(dat)=dat$patient_num
  dat[,nm.filter] = 1*(dat[,nm.filter]=="Y")
  nm.trans=setdiff(nm.cov, nm.notrans)

  p0 = length(nm.cov)
  dat.filter = dat[dat[,"patient_num"]%in%as.character(patient_num),c(nm.label, nm.cov)]
  
  N.f = nrow(dat.filter)
  
  dat.filter[, nm.trans] = log(dat.filter[, nm.trans] + 1)
  if(standardize==T){
  sd.filter = apply(dat.filter[, nm.trans], 2, sd)
  dat.filter[, nm.trans] = dat.filter[, nm.trans]/VTM(sd.filter, nrow(dat.filter))} # standardize
  return(dat.filter=dat.filter)
}


SUR.SEL.FUN=function(dat.filter, nm.label, nm.cov, nm.sur, rep0, rerunSelection){
  if(rerunSelection ==TRUE){
    p0 = length(nm.cov)
    N.f = nrow(dat.filter)
    K.sur = length(nm.sur)
    prev.filter = sum(is.element(dat.filter[,nm.label],c("Y", "P")))/sum(dat.filter[,nm.label]!="NULL")
    cut.sur = rep(0,K.sur)
    names(cut.sur) = nm.sur
    p0 = length(nm.cov)
    cut.sur = apply(dat.filter[,nm.sur],2,quantile,1-prev.filter)
    bet.mat.screen = matrix(100, nrow = K.sur * rep0, ncol = p0)
    colnames(bet.mat.screen) = nm.cov
    row.names(bet.mat.screen) = paste(rep(nm.sur, each=rep0), rep(1:rep0, K.sur), sep="_")
    set.seed(100)
    for (tmpnm in c(nm.sur)) {
      print(tmpnm)
      for (rr in 1:rep0) {
        print(rr)
        tmpdat = dat.filter[sample(1:N.f, min(N.f/2,1000)), c(tmpnm, setdiff(nm.cov, tmpnm))]
        tmpdat[,tmpnm] = 1 * (tmpdat[,tmpnm] > cut.sur[tmpnm])
        tmpbet = Est.ALASSO.GLM(data = as.matrix(tmpdat), Wi = NULL, rtn = "EST", regularize = T)
        bet.mat.screen[paste(tmpnm,rr,sep="_"), setdiff(nm.cov,tmpnm)] = tmpbet[2:length(nm.cov)]
      }
    }
  nm.cov.keep = nm.cov[apply(bet.mat.screen != 0, 2, mean) >= 1/2]
  save(bet.mat.screen,file=paste("bet.mat.screen ",dis,".RData",sep=""))
  }
  if(rerunSelection == FALSE){
  load(file=paste("bet.mat.screen ",dis,".RData",sep=""))
  nm.cov.keep = nm.cov[apply(bet.mat.screen != 0, 2, mean) >= 1/2]
  }
  return(list(bet.mat.screen=bet.mat.screen, nm.cov.keep=nm.cov.keep))
}

DAT.TRAIN.PREP.FUN=function(dat.filter, dat.pi, patient_num, nm.label, nm.cov.keep,
                            yes.case){
    dat = dat.filter[rownames(dat.filter)%in%as.character(patient_num),]
    dat.pi=dat.pi[rownames(dat.pi)%in%rownames(dat),]
    dat.sort=dat[order(rownames(dat)),]
    dat.pi.sort=dat.pi[order(rownames(dat.pi)),]
    dat=cbind(dat.sort, dat.pi.sort[,"pi"])
    colnames(dat)=c(colnames(dat.sort), "pi")
    dat[,nm.label] = 1*is.element(dat[,nm.label],yes.case)
    dat=dat[,c(nm.label, nm.cov.keep, "pi")]
    return(dat)
 }

convert2=function(fit) {
  rochat.auc = fit$rochat[1]
  rochat.values = matrix(fit$rochat[-1],ncol=6)
  colnames(rochat.values) = c("cutoff","pos.rate","FPR","TPR","PPV","NPV")
  roc.cv.auc = fit$roc.cv[1]
  roc.cv.values = matrix(fit$roc.cv[-1],ncol=6)
  colnames(roc.cv.values) = c("cutoff","pos.rate","FPR","TPR","PPV","NPV")
  betahat = fit$beta
  return(list(ROC.hat.auc=rochat.auc, ROC.hat.values=rochat.values, ROC.cv.auc=roc.cv.auc,ROC.cv.values=roc.cv.values, beta=betahat, Y.hat=fit$Si))
}

    
ROC.FUN.ALASSO.mycv  <- function(ind=NULL,data, wgti=NULL, Wi=NULL, yy0 = 0.5, nopen.ind=NULL,K=2,
                                 FPR0=seq(.01,.99,by=.01),rtn="ALL",yes.CV=F,yes.seed=F,rep=10,regularize=T,lam.ridge=0)
{
  ## ========================================================================= ##
  ## ==== ind: for bootstrap methods; data[ind,] returns bootstrap samples === ##
  ## ==== 1st column is binary disease status;                            ==== ##
  ## ==== 2nd column and on are individual markes                         ==== ##
  ## ==== yy0: prespecified cut-off;  FPR0: prespecified Specificity level ==== ##
  ## ========================================================================= ##
  
  nn <- nrow(data); if(is.null(wgti)){wgti=rep(1,nn)}; 
  if(!is.null(ind)){data <- data[ind,]; wgti = wgti[ind]}    
  n.set <- 1; pp = ncol(data); yyi.vmat = matrix(NA, nrow=nn, ncol=rep)
  
  ## ========================= ##
  ## === Apparent Accuracy === ##
  ## ========================= ##
  
  betahat = try(Est.ALASSO.GLM(data,Wi=wgti, rtn="EST",regularize=regularize,nopen.ind=nopen.ind),silent=T)
  if(length(betahat)>1)
  {
    yyi = g.logit(cbind(1,as.matrix(data[,-1]))%*%betahat[1:pp]); 
    q.nonzero = sum(abs(betahat[2:(pp)])!=0)
    if(q.nonzero==0){rochat = NULL}else{rochat = ROC.Est.FUN(data[,1],yyi,yy0,FPR0,wgti=wgti)}
    if(yes.seed){set.seed(315)}
    if(yes.CV & q.nonzero >0)
    {
      ## =============================== ##
      ## === K-fold cross validation === ##
      ## =============================== ##
      nv = floor(nn/K); roc.cv = NULL
      for(i in 1:rep){
        tmpind=sample(1:nn); data = data[tmpind,]; wgti = wgti[tmpind]
        for(k in 1:K){
          ind.v = 1:nv + (k-1)*nv; ind.t = setdiff(1:nn,ind.v)
          dat.t = data[ind.t,]; dat.v = data[ind.v,]            
          ## ============================================================================== ##
          ## ==== Calculating (1) Coef for Combining Markers (beta) with Training Data ==== ##
          ## ====             (2) Combined Marker Value (yyi.new) with Validation Data ==== ##
          ## ============================================================================== ##
          beta.t = try(Est.ALASSO.GLM(dat.t,Wi=wgti[ind.t],rtn="EST",regularize=regularize,nopen.ind=nopen.ind),silent=T)
          if(length(beta.t)>1){
            beta.t = beta.t[1:pp]
            yyi.v = g.logit(cbind(1,as.matrix(dat.v[,-1]))%*%beta.t);  yyi.vmat[ind.v,i]=yyi.v
            yyi.t = g.logit(cbind(1,as.matrix(dat.t[,-1]))%*%beta.t)
            roc.k = try(ROC.Est.FUN(dat.v[,1],yyi.v,yy0,FPR0,wgti[ind.v]),silent=T)
            if(length(roc.k)>1){roc.cv = cbind(roc.cv, roc.k)}
          }
        }
      }
      print(ncol(roc.cv)); auc.cv = roc.cv[1,]; roc.cv = apply(roc.cv,1,mean,na.rm=T)  
    }else{
      roc.cv = NULL
    }
  }else{
    roc.cv=beta.hat=NULL
  }
  if(rtn=="EST"){
    if(yes.CV){out=c(roc.cv,betahat)}else{out=c(rochat,betahat)}
    return(out)
  }else{return(list("rochat"=rochat, "roc.cv"=roc.cv, "beta"=betahat, "Si"=yyi, "auc.cv"=auc.cv, "yyi.v"=yyi.vmat))}
}
   

ROC.FUN.ALASSO.revise.mycv  <- function(ind=NULL,data, wgti=NULL, Wi=NULL, yy0 = 0.5, nopen.ind=NULL,K=2,
                                 FPR0=seq(.01,.99,by=.01),rtn="ALL",yes.CV=F,yes.seed=F,rep=10,regularize=T,lam.ridge=0)
{
  ## ========================================================================= ##
  ## ==== ind: for bootstrap methods; data[ind,] returns bootstrap samples === ##
  ## ==== 1st column is binary disease status;                            ==== ##
  ## ==== 2nd column and on are individual markes                         ==== ##
  ## ==== yy0: prespecified cut-off;  FPR0: prespecified Specificity level ==== ##
  ## ========================================================================= ##
  
  nn <- nrow(data); if(is.null(wgti)){wgti=rep(1,nn)}; 
  if(!is.null(ind)){data <- data[ind,]; wgti = wgti[ind]}    
  n.set <- 1; pp = ncol(data); yyi.vmat = matrix(NA, nrow=nn, ncol=rep)
  
  ## ========================= ##
  ## === Apparent Accuracy === ##
  ## ========================= ##
  
  betahat = try(Est.ALASSO.revise.GLM(data,Wi=wgti, rtn="EST",regularize=regularize,nopen.ind=nopen.ind),silent=T)
  if(length(betahat)>1)
  {
    yyi = g.logit(cbind(1,as.matrix(data[,-1]))%*%betahat[1:pp]); 
    q.nonzero = sum(abs(betahat[2:(pp)])!=0)
    if(q.nonzero==0){rochat = NULL}else{rochat = ROC.Est.FUN(data[,1],yyi,yy0,FPR0,wgti=wgti)}
    if(yes.seed){set.seed(315)}
    if(yes.CV & q.nonzero >0)
    {
      ## =============================== ##
      ## === K-fold cross validation === ##
      ## =============================== ##
      nv = floor(nn/K); roc.cv = NULL
      for(i in 1:rep){
        tmpind=sample(1:nn); data = data[tmpind,]; wgti = wgti[tmpind]
        for(k in 1:K){
          ind.v = 1:nv + (k-1)*nv; ind.t = setdiff(1:nn,ind.v)
          dat.t = data[ind.t,]; dat.v = data[ind.v,]            
          ## ============================================================================== ##
          ## ==== Calculating (1) Coef for Combining Markers (beta) with Training Data ==== ##
          ## ====             (2) Combined Marker Value (yyi.new) with Validation Data ==== ##
          ## ============================================================================== ##
          beta.t = try(Est.ALASSO.revise.GLM(dat.t,Wi=wgti[ind.t],rtn="EST",regularize=regularize,nopen.ind=nopen.ind),silent=T)
          if(length(beta.t)>1){
            beta.t = beta.t[1:pp]
            yyi.v = g.logit(cbind(1,as.matrix(dat.v[,-1]))%*%beta.t);  yyi.vmat[ind.v,i]=yyi.v
            yyi.t = g.logit(cbind(1,as.matrix(dat.t[,-1]))%*%beta.t)
            roc.k = try(ROC.Est.FUN(dat.v[,1],yyi.v,yy0,FPR0,wgti[ind.v]),silent=T)
            if(length(roc.k)>1){roc.cv = cbind(roc.cv, roc.k)}
          }
        }
      }
      print(ncol(roc.cv)); auc.cv = roc.cv[1,]; roc.cv = apply(roc.cv,1,mean,na.rm=T)  
    }else{
      roc.cv = NULL
    }
  }else{
    roc.cv=beta.hat=NULL
  }
  if(rtn=="EST"){
    if(yes.CV){out=c(roc.cv,betahat)}else{out=c(rochat,betahat)}
    return(out)
  }else{return(list("rochat"=rochat, "roc.cv"=roc.cv, "beta"=betahat, "Si"=yyi, "auc.cv"=auc.cv, "yyi.v"=yyi.vmat))}
}


Est.ALASSO.POR=function(y,x){
fit.cv <- try(ordinalNetCV(x, y, family="cumulative", tuneMethod = "bic", link="logit", alpha=1,
                     parallelTerms=TRUE, nonparallelTerms=FALSE), silent=T)

fit <- ordinalNet(x, y, family="cumulative", link="logit", alpha=1,
                   parallelTerms=TRUE, nonparallelTerms=FALSE, lambdaVals=summary(fit.cv)$lambda[1])
coef(fit)
}

convert2score=function(data, betahat, K){
  b0.hat=betahat[1:(K-1)]
  b.hat=betahat[-c(1:(K-1))]
  myscore <- c(as.matrix(data[,-1])%*% b.hat) + matrix(b0.hat, nrow=dim(data)[1], ncol=K-1, byrow=TRUE)
  invlogit <- function(x) 1 / (1+exp(-x))
  cumprob <- t(apply(myscore, 1, invlogit))
  prob <- cbind(cumprob, 1) - cbind(0, cumprob)
  myscore2=prob[,1]*0+prob[,2]*0.5+prob[,3]*1
  return(myscore2)
}


ROC.FUN.ALASSO.ordinal=function(ind=NULL,data, wgti=NULL, Wi=NULL, yy0 = 0.5, nopen.ind=NULL,nfold=5,K,
         FPR0=seq(.01,.99,by=.01),rtn="ALL",yes.CV=T,yes.seed=F,rep=10,regularize=T,lam.ridge=0)
{
  ## ========================================================================= ##
  ## ==== ind: for bootstrap methods; data[ind,] returns bootstrap samples === ##
  ## ==== 1st column is binary disease status;                            ==== ##
  ## ==== 2nd column and on are individual markes                         ==== ##
  ## ==== yy0: prespecified cut-off;  FPR0: prespecified Specificity level ==== ##
  ## ========================================================================= ##
  
  nn <- nrow(data); if(is.null(wgti)){wgti=rep(1,nn)}; 
  if(!is.null(ind)){data <- data[ind,]; wgti = wgti[ind]}    
  n.set <- 1; pp = ncol(data); yyi.vmat = matrix(NA, nrow=nn, ncol=rep)
  
  ## ========================= ##
  ## === Apparent Accuracy === ##
  ## ========================= ##
  y=as.factor(data[,1]);x=as.matrix(data[,-1])
  betahat=Est.ALASSO.POR(y,x)
  if(length(betahat)>1)
  {
    yyi = convert2score(data, betahat, K)
    q.nonzero = sum(abs(betahat[-(1:(K-1))])!=0)
    if(q.nonzero==0){rochat = NULL}else{rochat = ROC.Est.FUN(data[,1],yyi,yy0,FPR0,wgti=wgti)}
    if(yes.seed){set.seed(315)}
    if(yes.CV & q.nonzero >0)
    {
      ## =============================== ##
      ## === K-fold cross validation === ##
      ## =============================== ##
      nv = floor(nn/nfold); roc.cv = NULL
      for(i in 1:rep){
        tmpind=sample(1:nn); data = data[tmpind,]; wgti = wgti[tmpind]
        for(k in 1:nfold){
          ind.v = 1:nv + (k-1)*nv; ind.t = setdiff(1:nn,ind.v)
          dat.t = data[ind.t,]; dat.v = data[ind.v,]            
          ## ============================================================================== ##
          ## ==== Calculating (1) Coef for Combining Markers (beta) with Training Data ==== ##
          ## ====             (2) Combined Marker Value (yyi.new) with Validation Data ==== ##
          ## ============================================================================== ##
          beta.t = try(Est.ALASSO.POR(as.factor(dat.t[,1]), as.matrix(dat.t[,-1])),silent=T)
          if(length(beta.t)>1){
            yyi.v = convert2score(dat.v, beta.t, K);  yyi.vmat[ind.v,i]=yyi.v
            yyi.t = convert2score(dat.t, beta.t, K);  yyi.vmat[ind.v,i]=yyi.v
            roc.k = try(ROC.Est.FUN(dat.v[,1],yyi.v,yy0,FPR0,wgti[ind.v]),silent=T)
            if(length(roc.k)>1){roc.cv = cbind(roc.cv, roc.k)}
          }
        }
      }
      print(ncol(roc.cv)); auc.cv = roc.cv[1,]; roc.cv = apply(roc.cv,1,mean,na.rm=T)  
      
    }
  }
return(list("rochat"=rochat, "roc.cv"=roc.cv,  "beta"=betahat))
}

ROC.FUN.ALASSO.mycv2=function(ind=NULL,data, nm.label.train, nm.label.validation, wgti.train=NULL, wgti.validation=NULL, Wi=NULL, yy0 = 0.5, nopen.ind=NULL,K=2,
         FPR0=seq(.01,.99,by=.01),rtn="ALL",yes.CV=T,yes.seed=T,rep=10,regularize=T,lam.ridge=0)
{
  ## ========================================================================= ##
  ## ==== ind: for bootstrap methods; data[ind,] returns bootstrap samples === ##
  ## ==== 1st column is binary disease status;                            ==== ##
  ## ==== 2nd column and on are individual markes                         ==== ##
  ## ==== yy0: prespecified cut-off;  FPR0: prespecified Specificity level ==== ##
  ## ========================================================================= ##
  
  nn <- nrow(data); if(is.null(wgti.train)){wgti.train=rep(1,nn)}; if(is.null(wgti.validation)){wgti.validation=rep(1,nn)}; 
  if(!is.null(ind)){data <- data[ind,]; wgti.train = wgti.train[ind]; wgti.validation = wgti.validation[ind]}    
  n.set <- 1; pp = ncol(data)-1; yyi.vmat = matrix(NA, nrow=nn, ncol=rep)
  
  ## ========================= ##
  ## === Apparent Accuracy === ##
  ## ========================= ##
  id.rm.train=which(is.na(data[,nm.label.train])==1)
  betahat = try(Est.ALASSO.GLM(data[-id.rm.train,setdiff(colnames(data), nm.label.validation)],Wi=wgti.train[-id.rm.train], rtn="EST",regularize=regularize,nopen.ind=nopen.ind),silent=T)
  if(length(betahat)>1)
  {
    yyi = g.logit(cbind(1,as.matrix(data[,setdiff(colnames(data), c(nm.label, nm.label2))]))%*%betahat[1:pp]); 
    q.nonzero = sum(abs(betahat[2:(pp)])!=0)
    if(q.nonzero==0){rochat = NULL}else{rochat = ROC.Est.FUN(data[,nm.label.validation],yyi,yy0,FPR0,wgti=wgti.validation)}
    if(yes.seed){set.seed(315)}
    if(yes.CV & q.nonzero >0)
    {
      ## =============================== ##
      ## === K-fold cross validation === ##
      ## =============================== ##
      nv = floor(nn/K); roc.cv = NULL
      for(i in 1:rep){
        tmpind=sample(1:nn); data = data[tmpind,]; wgti.train = wgti.train[tmpind]; wgti.validation = wgti.validation[tmpind]
        for(k in 1:K){
          ind.v = 1:nv + (k-1)*nv; ind.t = setdiff(1:nn,ind.v)
          dat.t = data[ind.t,setdiff(colnames(data), nm.label2)]; dat.v = data[ind.v,setdiff(colnames(data), nm.label)]   
          id.rm.train.t=which(is.na(dat.t[,nm.label])==1)
          ## ============================================================================== ##
          ## ==== Calculating (1) Coef for Combining Markers (beta) with Training Data ==== ##
          ## ====             (2) Combined Marker Value (yyi.new) with Validation Data ==== ##
          ## ============================================================================== ##
          beta.t = try(Est.ALASSO.GLM(dat.t[-id.rm.train.t, ],Wi=wgti.train[ind.t[-id.rm.train.t]],rtn="EST",regularize=regularize,nopen.ind=nopen.ind),silent=T)
          if(length(beta.t)>1){
            beta.t = beta.t[1:pp]
            yyi.v = g.logit(cbind(1,as.matrix(dat.v[,-1]))%*%beta.t);  yyi.vmat[ind.v,i]=yyi.v
            yyi.t = g.logit(cbind(1,as.matrix(dat.t[,-1]))%*%beta.t)
            roc.k = try(ROC.Est.FUN(dat.v[,nm.label2],yyi.v,yy0,FPR0,wgti.validation[ind.v]),silent=T)
            if(length(roc.k)>1){roc.cv = cbind(roc.cv, roc.k)}
          }
        }
      }
      #print(ncol(roc.cv)); 
      auc.cv = roc.cv[1,]; roc.cv = apply(roc.cv,1,mean,na.rm=T)  
    }else{
      roc.cv = NULL
    }
  }else{
    roc.cv=beta.hat=NULL
  }
  if(rtn=="EST"){
    if(yes.CV){out=c(roc.cv,betahat)}else{out=c(rochat,betahat)}
    return(out)
  }else{return(list("rochat"=rochat, "roc.cv"=roc.cv, "beta"=betahat, "Si"=yyi, "auc.cv"=auc.cv, "yyi.v"=yyi.vmat))}
}


Est.ALASSO.GLM.GAUSSIAN=function(data,Wi=NULL,rtn="EST",nopen.ind=NULL,regularize=yes.regularize,pwr=0.1){
  data = as.matrix(data); y = data[,1]; x = data[,-1,drop=F]; nn=length(y); if(is.null(Wi)){Wi=rep(1,nn)}; pp = ncol(x)
  if(regularize){
    ## adaptive lasso (aLASSO) with initial estimator obtained via ridge (bini); w.b creates adaptive weights for aLASSO ##
    lam.ridge = pp/nn; ##bini = plr(x,y,lambda=lam.ridge,weights=Wi)$coef; 
    bini = as.vector(coef(glmnet(x,y,weights=Wi,alpha=0,standardize=F,lambda=lam.ridge,family="gaussian"))); ##print(bini)
    w.b = 1/abs(bini[-1]); x.t = x/VTM(w.b,nrow(x))
    
    ## glmpath provides solution path for a range of penalty parameters ##
    tmpfit = glmpath(x.t,y,nopenalty.subset=nopen.ind,family=gaussian,weight=Wi,standardize=F,min.lambda=0,lambda2=lam.ridge)
    lam.all = c(seq(min(tmpfit$lambda),max(tmpfit$lambda),length=500))
    b.all = predict(tmpfit, s=lam.all, type="coefficients",mode="lambda")
    b.all = b.all/VTM(c(1,w.b),nrow(b.all)); m0 = length(lam.all)
    ## ================================================================================ ##
    ## calculates degree of freedom for all beta's (corresponding to different lam.all) ##
    ## ================================================================================ ##
    df.all = apply(b.all[,-1,drop=F]!=0,1,sum); x = as.matrix(x)
    ## =============================================================================================================== ##
    ## calculates modified BIC, log(n) is modified as min{sum(y)^0.1, log(n)} to avoid over shrinkage in finite sample ##
    ## =============================================================================================================== ##
    BIC.lam = -2*apply(predict(tmpfit,newx=x.t,newy=y,s=lam.all,type="loglik",mode="lambda"),2,sum)+min(sum(y)^pwr,log(sum(y)))*df.all 
    m.opt = (1:m0)[BIC.lam==min(BIC.lam)]; bhat = b.all[m.opt,]; lamhat = lam.all[m.opt]
  }else{
    ## ========================================================================= ##
    ## if regularize = F, then use standard logistic regression w/o penalization ##
    ## ========================================================================= ##
    bhat=bini=glm(y~x,family=gaussian,weight=Wi)$coef;lamhat = 0; lam.all=BIC.lam=b.all=NULL
  }
  out = c("b"=bhat, "bini"=bini,"lamhat"=lamhat,"lam.all"=lam.all,"BIC.lam"=BIC.lam,"b.all"=b.all)
  if(rtn=="EST"){return(out)}else{return(list(out,"b.all"=b.all,"lam.all"=lam.all,"fit"=tmpfit,"BIC.lam"=BIC.lam))}
}

Est.ALASSO.revise.GLM = function(data,Wi=NULL,rtn="EST",nopen.ind=NULL,regularize=yes.regularize,pwr=0.1){
  data = as.matrix(data); y = data[,1]; x = data[,-1,drop=F]; nn=length(y); if(is.null(Wi)){Wi=rep(1,nn)}; pp = ncol(x)
  if(regularize){
    ## adaptive lasso (aLASSO) with initial estimator obtained via ridge (bini); w.b creates adaptive weights for aLASSO ##
    lam.ridge = log(pp)/nn; ##bini = plr(x,y,lambda=lam.ridge,weights=Wi)$coef; 
    bini = as.vector(coef(glmnet(x,y,weights=Wi,alpha=0,standardize=F,lambda=lam.ridge,family="binomial"))); ##print(bini)
    w.b = 1/abs(bini[-1]); x.t = x/VTM(w.b,nrow(x))
    
    ## glmpath provides solution path for a range of penalty parameters ##
    tmpfit = glmpath(x.t,y,nopenalty.subset=nopen.ind,family=binomial,weight=Wi,standardize=F,min.lambda=0,lambda2=lam.ridge)
    lam.all = c(seq(min(tmpfit$lambda),max(tmpfit$lambda),length=500))
    b.all = predict(tmpfit, s=lam.all, type="coefficients",mode="lambda")
    b.all = b.all/VTM(c(1,w.b),nrow(b.all)); m0 = length(lam.all)
    ## ================================================================================ ##
    ## calculates degree of freedom for all beta's (corresponding to different lam.all) ##
    ## ================================================================================ ##
    df.all = apply(b.all[,-1,drop=F]!=0,1,sum); x = as.matrix(x)
    ## =============================================================================================================== ##
    ## calculates modified BIC, log(n) is modified as min{sum(y)^0.1, log(n)} to avoid over shrinkage in finite sample ##
    ## =============================================================================================================== ##
    BIC.lam = -2*apply(predict(tmpfit,newx=x.t,newy=y,s=lam.all,type="loglik",mode="lambda"),2,sum)+min(sum(y)^pwr,log(sum(y)))*df.all 
    m.opt = (1:m0)[BIC.lam==min(BIC.lam)]; bhat = b.all[m.opt,]; lamhat = lam.all[m.opt]
  }else{
    ## ========================================================================= ##
    ## if regularize = F, then use standard logistic regression w/o penalization ##
    ## ========================================================================= ##
    bhat=bini=glm(y~x,family=binomial,weight=Wi)$coef;lamhat = 0; lam.all=BIC.lam=b.all=NULL
  }
  out = c("b"=bhat, "bini"=bini,"lamhat"=lamhat,"lam.all"=lam.all,"BIC.lam"=BIC.lam,"b.all"=b.all)
  if(rtn=="EST"){return(out)}else{return(list(out,"b.all"=b.all,"lam.all"=lam.all,"fit"=tmpfit,"BIC.lam"=BIC.lam))}
}

    