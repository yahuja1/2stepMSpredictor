# install.packages(c("MASS","glmpath"), repos = "http://cran.us.r-project.org")
library(MASS); library(glmpath); library(glmnet)
g.logit = function(xx){exp(xx)/(exp(xx)+1)}
logit = function(xx){log(xx/(1-xx))}
dg.logit = function(xx){exp(xx)/(exp(xx)+1)^2}
Intercept.GLM = function(yi,bxi){glm(yi~bxi,family="binomial")}

ROC.FUN.Sup.Par = function(St,Yt,fpr=jump.u){
  mhati = g.logit(St); ss = unique(sort(St)); mu1 = mean(Yt); mu0 = 1-mu1; nt = length(Yt)
  S1.ss = sum.I(ss, "<=", St, mhati)/sum(mhati)
  S0.ss = sum.I(ss, "<=", St, 1-mhati)/sum(1-mhati); auc = sum(S1.ss[-1]*(S0.ss[-nt]-S0.ss[-1]))
  PPV.ss = S1.ss*mu1/(S1.ss*mu1+S0.ss*mu0); NPV.ss = (1-S0.ss)*mu0/((1-S0.ss)*mu0+(1-S1.ss)*mu1)
  out=cbind("cut"=g.logit(ss),"FPR"=S0.ss,"TPR"=S1.ss,"PPV"=PPV.ss,"NPV"=NPV.ss); tmpnm=colnames(out)
  out=sapply(1:ncol(out),function(kk){approx(out[,"FPR"],out[,kk],fpr)$y}); colnames(out)=tmpnm
  list(auc,out)
}


ROC.FUN.SSL.NP = function(Sv,St,Yt,fpr=jump.u,yes.CV=F,ss0=NULL, Xt=NULL,Xv=NULL,rep=10,regularize=yes.regularize){
  nv = length(Sv); nt = length(Yt); Pt = sum.I(St, ">=", Sv)/nv
  Pv = sum.I(Sv,">=",Sv)/nv; bw = sd(Pt)/nt^0.3; ss = c(ss0,sort(Pv))
  mhat.v = predict(locfit(Yt~lp(Pt, deg=0, h=bw),ev=Pv))
  mu1 = mean(Yt); mu0 = 1-mu1; 
  if(!yes.CV){
    S1.ss = sum.I(ss, "<=", Pv, mhat.v)/sum(mhat.v)
    S0.ss = sum.I(ss, "<=", Pv, 1-mhat.v)/sum(1-mhat.v); 
  }else{
    K.fold = 2; nt.t = floor(nt/K.fold); S1.ss = S0.ss = matrix(NA,nrow=length(ss),ncol=K.fold*rep); dat.t = cbind(Yt,Xt)
    for(rr in 1:rep){
      dat.t = dat.t[sample(1:nt),]
      for(kk in 1:K.fold){
        indt.v = 1:nt.t + (kk-1)*nt.t; indt.t = setdiff(1:nt, indt.v)
        bhat.t = Est.ALASSO.GLM(dat.t[indt.t,],rtn="EST",regularize=regularize)[1+0:ncol(Xv)];
        bhat.v = Est.ALASSO.GLM(dat.t[indt.v,],rtn="EST",regularize=regularize)[1+0:ncol(Xv)];
        mvi = g.logit(cbind(1,Xv)%*%bhat.t); Pvi = sum.I(cbind(1,Xv)%*%bhat.v,">=",Sv)/nv
        S1.ss[,kk+K.fold*(rr-1)] = sum.I(ss,"<=",Pvi,mvi)/sum(mvi)
        S0.ss[,kk+K.fold*(rr-1)] = sum.I(ss,"<=",Pvi,1-mvi)/sum(1-mvi)
      }}
    S1.ss = apply(S1.ss,1,mean); S0.ss = apply(S0.ss,1,mean)
  }  
  auc = sum(S1.ss[length(ss0)+2:nv]*(S0.ss[length(ss0)+1:(nv-1)]-S0.ss[length(ss0)+2:nv]))
  PPV.ss = S1.ss*mu1/(S1.ss*mu1+S0.ss*mu0); NPV.ss = (1-S0.ss)*mu0/((1-S0.ss)*mu0+(1-S1.ss)*mu1)
  out=cbind("cut"=ss,"FPR"=S0.ss,"TPR"=S1.ss,"PPV"=PPV.ss,"NPV"=NPV.ss); tmpnm=colnames(out)
  out2=sapply(1:ncol(out),function(kk){approx(out[,"FPR"],out[,kk],fpr,rule=2)$y}); colnames(out2)=tmpnm
  list(auc,rbind(out[-(length(ss0)+1:nv),],out2))
}


ROC.FUN.SSL.Par = function(Xv,dat.t,fpr=jump.u,yes.CV=F,rep=10,wgt=NULL,rtn="list",bptb=NULL,bptb2=NULL,regularize=yes.regularize){
  if(is.null(wgt)){wgt=rep(1,nrow(dat.t))}
  if(is.null(bptb)){bptb = Est.ALASSO.GLM(dat.t,rtn="EST",regularize=regularize,Wi=wgt)[1+0:ncol(Xv)]}
  Yt = dat.t[,1]; if(is.null(bptb2)){bptb2 = bptb}; 
  Sv = c(cbind(1,Xv)%*%bptb); mhat.v = g.logit(c(cbind(1,Xv)%*%bptb2)); 
  ss = unique(sort(Sv)); mu1 = mean(mhat.v); mu0 = 1-mu1; nv = length(Sv); nt = length(Yt)
  if(!yes.CV){
    S1.ss = sum.I(ss, "<=", Sv, mhat.v)/sum(mhat.v)
    S0.ss = sum.I(ss, "<=", Sv, 1-mhat.v)/sum(1-mhat.v); 
  }else{
    K.fold = 2; nt.t = floor(nt/K.fold); S1.ss = S0.ss = matrix(NA,nrow=length(ss),ncol=K.fold*rep)
    for(rr in 1:rep){
      dat.t = dat.t[sample(1:nt),]
      for(kk in 1:K.fold){
        indt.v = 1:nt.t + (kk-1)*nt.t; indt.t = setdiff(1:nt, indt.v)
        bhat.t = Est.ALASSO.GLM(dat.t[indt.t,],rtn="EST",regularize=regularize,BIC.power=0.1)[1+0:ncol(Xv)];
        bhat.v = Est.ALASSO.GLM(dat.t[indt.v,],rtn="EST",regularize=regularize,BIC.power=0.1)[1+0:ncol(Xv)];
        mvi = g.logit(cbind(1,Xv)%*%bhat.t); Svi = cbind(1,Xv)%*%bhat.v
        S1.ss[,kk+K.fold*(rr-1)] = sum.I(ss,"<=",Svi,mvi)/sum(mvi)
        S0.ss[,kk+K.fold*(rr-1)] = sum.I(ss,"<=",Svi,1-mvi)/sum(1-mvi)
    }}
    S1.ss = apply(S1.ss,1,mean); S0.ss = apply(S0.ss,1,mean)
  }  
  auc = sum(S1.ss[-1]*(S0.ss[-nv]-S0.ss[-1]))
  PPV.ss = S1.ss*mu1/(S1.ss*mu1+S0.ss*mu0); NPV.ss = (1-S0.ss)*mu0/((1-S0.ss)*mu0+(1-S1.ss)*mu1)
  out=cbind("cut"=g.logit(ss),"FPR"=S0.ss,"TPR"=S1.ss,"PPV"=PPV.ss,"NPV"=NPV.ss); tmpnm=colnames(out)
  out=sapply(1:ncol(out),function(kk){approx(out[,"FPR"],out[,kk],fpr,rule=2)$y}); colnames(out)=tmpnm
  if(rtn=="vec"){return(c(auc,out))}else{return(list(auc,out))}
}


## ================================================================================================================== ##
## estimate beta w/ adaptive LASSO regularization (if regularize=T) or standard logistic regression (if regularize=F) ##
## data: 1st column y; remaining x; nopen.ind indexes which subset of x should not be penalized                       ##
## Wi: weights for resampling if interested in obtaining standard errors                                              ## 
## ================================================================================================================== ##
Est.ALASSO.GLM = function(data,Wi=NULL,rtn="EST",nopen.ind=NULL,regularize=yes.regularize,pwr=0.1){
  data = as.matrix(data); y = data[,1]; x = data[,-1,drop=F]; nn=length(y); if(is.null(Wi)){Wi=rep(1,nn)}; pp = ncol(x)
  if(regularize){
    ## adaptive lasso (aLASSO) with initial estimator obtained via ridge (bini); w.b creates adaptive weights for aLASSO ##
    lam.ridge = pp/nn; ##bini = plr(x,y,lambda=lam.ridge,weights=Wi)$coef; 
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


Est.ALASSO.GLM2 = function(data,Wi=NULL,rtn="EST",nopen.ind=NULL,regularize=yes.regularize, pwr=0.1){
  data = as.matrix(data); y = data[,1]; x = data[,-1,drop=F]; nn=length(y); if(is.null(Wi)){Wi=rep(1,nn)}; pp = ncol(x)
  if(regularize){
    ## adaptive lasso (aLASSO) with initial estimator obtained via ridge (bini); w.b creates adaptive weights for aLASSO ##
    lam.ridge = pp/nn; ##bini = plr(x,y,lambda=lam.ridge,weights=Wi)$coef; 
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


SIM.FUN = function(nn,rtn="data.t")
  {
	xx = mvrnorm(nn,mu=rep(0,p.x),Sigma=Sig0.X); 
	icd.B = rbinom(nn,size=1,prob=g.logit(xx[,1]*3+2))
	xx[,1] = (xx[,1]*(xx[,1]>0)+(rexp(nn,rate=0.1)+5)*rbinom(nn,size=1,prob=0.1))*icd.B
	prob.x = g.logit(-alp0+c(xx%*%beta0)-3*(1-icd.B)+0.2*(xx[,1]>15))
	yy = rbinom(nn,prob=prob.x,size=1); dat = cbind(yy,xx)
	if(rtn=="data.t"){return(dat)}else{
		zz = rbinom(nn, size=2, prob=g.logit(log(maf)+gam.z*yy))
		return(list(dat,cbind("D"=yy,"P.x"=prob.x,"G"=zz)))}
  }

Est.ALASSO.GLM.new = function(data,Wi=NULL,fit.type=fit.typ,nopen.ind=NULL,BIC.factor=0.1,offset=NULL,regularize=T){
  if(fit.type=="exact"){
    return(Est.ALASSO.GLM(data,Wi=Wi,nopen.ind=nopen.ind,BIC.factor=BIC.factor,offset=offset,fam0="binomial",regularize=regularize))
  }else{ 
    return(Est.ALASSO.GLM.Approx(data,Wi=Wi,nopen.ind=nopen.ind,BIC.factor=BIC.factor,offset=offset))
  }
}

## ================================================================================================================== ##
## estimate beta w/ adaptive LASSO regularization (if regularize=T) or standard logistic regression (if regularize=F) ##
## data: 1st column y; remaining x; nopen.ind indexes which subset of x should not be penalized                       ##
## Wi: weights for resampling if interested in obtaining standard errors                                              ## 
## ================================================================================================================== ##
Est.ALASSO.GLM.Approx = function(data,Wi=NULL,rtn="EST",nopen.ind=NULL,BIC.factor=0.1,offset=NULL){
  data = as.matrix(data); y = data[,1]; x = data[,-1,drop=F]; nn=length(y); pp = ncol(x)
  if(is.null(Wi)){Wi=rep(1,nn)}; if(is.null(offset)){offset=rep(0,nn)}
  btilde = glm(y~x,family=binomial,weight=Wi)
  Ahat = solve(summary(btilde)$cov.unscaled)
  btilde = btilde$coef
  Ahat.half = svd(Ahat); Ahat.half = Ahat.half$u%*%diag(sqrt(Ahat.half$d))%*%t(Ahat.half$u)
  Xtilde = Ahat.half; Ytilde = Ahat.half%*%btilde
  w.b = 1/abs(btilde); Xtilde.t = Xtilde/VTM(w.b,nrow(Xtilde))
  tmpfit = lars(Xtilde.t,Ytilde,type="lasso",normalize=F,intercept=F)
  lam.all = c(seq(min(tmpfit$lambda),max(tmpfit$lambda),length=500))
  b.all = predict(tmpfit, s=lam.all, type="coefficients",mode="lambda")$coef
  b.all = b.all/VTM(w.b,nrow(b.all)); m0 = length(lam.all)
  df.all = apply(b.all[,-1,drop=F]!=0,1,sum)+1; 
  BIC.lam = -2*logitlik.fun(t(b.all),dat=data)+min(nn^BIC.factor,log(nn))*df.all 
  m.opt = (1:m0)[BIC.lam==min(BIC.lam)]; bhat = b.all[m.opt,]; lamhat = lam.all[m.opt]
  return(bhat)
}


ROC.FUN.ALASSO.632boot  <- function(data, ind=NULL, wgti=NULL, Wi=NULL, yy0 = 0.5, nopen.ind=NULL,
                            FPR0=seq(.01,.99,by=.01),rtn="ALL",yes.CV=F,yes.seed=F,rep=10,regularize=T){
  ## ========================================================================= ##
  ## ==== ind: for bootstrap methods; data[ind,] returns bootstrap samples === ##
  ## ==== 1st column is binary disease status;                            ==== ##
  ## ==== 2nd column and on are individual markes                         ==== ##
  ## ==== yy0: prespecified cut-off;  FPR0: prespecified Specificity level ==== ##
  ## ========================================================================= ##
  nn <- nrow(data); if(is.null(wgti)){wgti=rep(1,nn)}; 
  if(!is.null(ind)){data <- data[ind,]; wgti = wgti[ind]}    
  n.set <- 1; pp = ncol(data);  yyi.vmat = matrix(NA, nrow=nn, ncol=rep)
  ## === Apparent Accuracy === ##
  betahat = try(Est.ALASSO.GLM(data,Wi=wgti, rtn="EST",regularize=regularize,nopen.ind=nopen.ind),silent=T)
  if(length(betahat)>1){
    yyi = g.logit(cbind(1,as.matrix(data[,-1]))%*%betahat[1:pp])    
    q.nonzero = sum(abs(betahat[2:(pp)])!=0)
    if(q.nonzero==0){rochat = NULL}else{rochat = ROC.Est.FUN(data[,1],yyi,yy0,FPR0,wgti=wgti)}
    if(yes.seed){set.seed(315)}
    if(yes.CV&q.nonzero>0){
      ## === 0.632 bootstrap === ##
      roc.cv = NULL; 
      for(i in 1:rep){
        tmpind=sample(1:nn,replace=T); ind.v = setdiff(1:nn,unique(tmpind))
        dat.t = data[tmpind,]; wgti.t = wgti[tmpind]; dat.v = data[ind.v,];  wgti.v = wgti[ind.v]
        beta.t = try(Est.ALASSO.GLM(dat.t,Wi=wgti.t,rtn="EST",regularize=regularize,nopen.ind=nopen.ind),silent=T)
        if(length(beta.t)>1){
          beta.t = beta.t[1:pp]; yyi.v = g.logit(cbind(1,as.matrix(dat.v[,-1]))%*%beta.t); yyi.vmat[ind.v,i]=yyi.v
          roc.k = try(ROC.Est.FUN(dat.v[,1],yyi.v,yy0,FPR0,wgti.v),silent=T)
          if(length(roc.k)>1){roc.cv = cbind(roc.cv, roc.k)}  } }
       roc.cv = apply(roc.cv,1,mean)*0.632+rochat*0.368  
     }else{roc.cv = NULL}
   }else{roc.cv=beta.hat=NULL}
   if(rtn=="EST"){
     if(yes.CV){out=c(roc.cv,betahat)}else{out=c(rochat,betahat)}
     return(out)
   }else{return(list("rochat"=rochat, "roc.cv"=roc.cv, "beta"=betahat, "Si"=yyi, "yyi.v"=yyi.vmat))}
}

ROC.FUN.ALASSO.cv  <- function(ind=NULL,data, wgti=NULL, Wi=NULL, yy0 = 0.5, nopen.ind=NULL,K=2,
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
  
  betahat = try(Est.ALASSO.GLM.new(data,Wi=wgti, rtn="EST",regularize=regularize,nopen.ind=nopen.ind),silent=T)
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
          beta.t = try(Est.ALASSO.GLM.new(dat.t,Wi=wgti[ind.t],rtn="EST",regularize=regularize,nopen.ind=nopen.ind),silent=T)
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


ROC.FUN.ALASSO.cv2  <- function(ind=NULL,data, wgti=NULL, Wi=NULL, yy0 = 0.5, nopen.ind=NULL,K=2,
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
    n.set <- 1; pp = ncol(data)

    ## ========================= ##
    ## === Apparent Accuracy === ##
    ## ========================= ##

  	betahat = try(Est.ALASSO.GLM(data,Wi=wgti, rtn="EST",regularize=regularize,nopen.ind=nopen.ind,lam.ridge=lam.ridge),silent=T)
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
            nv = floor(nn/K); bias.cv = NULL
            for(i in 1:rep)
              {
                tmpind=sample(1:nn); data = data[tmpind,]; wgti = wgti[tmpind]
                for(k in 1:K)
                 {
                    ind.v = 1:nv + (k-1)*nv; ind.t = setdiff(1:nn,ind.v)
                    dat.t = data[ind.t,]; dat.v = data[ind.v,]            
                    ## ============================================================================== ##
                    ## ==== Calculating (1) Coef for Combining Markers (beta) with Training Data ==== ##
                    ## ====             (2) Combined Marker Value (yyi.new) with Validation Data ==== ##
                    ## ============================================================================== ##
                    beta.t = try(Est.ALASSO.GLM(dat.t,Wi=wgti[ind.t],rtn="EST",regularize=regularize,nopen.ind=nopen.ind,lam.ridge=lam.ridge),silent=T)
                    if(length(beta.t)>1)
                      {
                        beta.t = beta.t[1:pp]
                        yyi.v = g.logit(cbind(1,as.matrix(dat.v[,-1]))%*%beta.t)
                        yyi.t = g.logit(cbind(1,as.matrix(dat.t[,-1]))%*%beta.t)
                        bias.k = try(ROC.Est.FUN(dat.t[,1],yyi.t,yy0,FPR0,wgti[ind.t])-ROC.Est.FUN(dat.v[,1],yyi.v,yy0,FPR0,wgti[ind.v]),silent=T)
                        if(length(bias.k)>1){bias.cv = cbind(bias.cv, bias.k)}
                      }
                }
              }
            print(ncol(bias.cv)); bias.cv = apply(bias.cv,1,mean,trim=0.05,na.rm=T)/2; roc.cv = rochat - bias.cv  
        }else{
            roc.cv = NULL
            }
     }else{
      roc.cv=beta.hat=NULL
     }
   if(rtn=="EST"){
     if(yes.CV){out=c(roc.cv,betahat)}else{out=c(rochat,betahat)}
     return(out)
   }else{return(list(rochat, roc.cv, betahat, yyi))}
  }

ROC.Est.FUN <- function(Di,yyi,yy0,fpr0=NULL,wgti=NULL,yes.smooth=F)
  {
    out.yy <- out.pp <- out.AUC <- out.TPR <- out.FPR <- out.PPV <- out.NPV <- NULL
    if(is.null(wgti)){wgti=rep(1,length(Di))}; yyi = as.matrix(yyi); pp=ncol(as.matrix(yyi));  
    mu0 = sum(wgti*(1-Di))/sum(wgti); mu1 = 1-mu0  
    for(k in 1:pp)
      {
       yy = yy0;
       if(!is.null(fpr0)){
         tpr.all = S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth);
         fpr.all = S.FUN(yyi[,k],Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth);
         TPR = approx(c(0,fpr.all,1),c(0,tpr.all,1),fpr0,method="linear",rule=2)$y;
         TPR = c(S.FUN(yy0,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth), TPR);
          yy = c(yy,Sinv.FUN(fpr0,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth))        
         FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
       }else{
         TPR = S.FUN(yy,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth);
         FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
       }
       out.yy = cbind(out.yy, yy)
       out.pp = cbind(out.pp, S.FUN(yy,Yi=yyi[,k],wgti,yes.smooth=yes.smooth))
       out.TPR = cbind(out.TPR,  TPR);  out.FPR  <- cbind(out.FPR,  FPR)
       PPV <- 1/(1+FPR*mu0/(TPR*mu1)); NPV <- 1/(1+(1-TPR)*mu1/((1-FPR)*mu0))
       out.PPV <- cbind(out.PPV, PPV); out.NPV <- cbind(out.NPV, NPV)
       #AUC <- sum((sum.I(yyi[,k],"<=",Yi=yyi[,k],Vi=Di*wgti)+sum.I(yyi[,k],"<",Yi=yyi[,k],Vi=Di*wgti))*(1-Di)*wgti/2
       #             )/(sum((1-Di)*wgti)*sum(Di*wgti))
       AUC = sum(S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth)*(1-Di)*wgti)/sum((1-Di)*wgti)
       out.AUC <- c(out.AUC, AUC)
     }
    out = list("AUC"=out.AUC,"FPR"=out.FPR,"TPR"=out.TPR,"PPV"=out.PPV,"NPV"=out.NPV)
    out
  }

AUC.FUN = function(data)
  {
	dd = data[,1]; xx = data[,2]; n0 = sum(1-dd); n1 = sum(dd) 
	x0 = xx[dd==0]; x1 = xx[dd==1]
  	sum((sum.I(x0, "<=", x1)+sum.I(x0,"<",x1))/2)/(n0*n1)
  }

chat.FUN.632boot <- function(data,fpr0=0.05,betahat,regularize=T,nopen.ind=NULL)
  {
    yy = data[,1]; xx = as.matrix(data[,-1]); p0 = ncol(data); data = as.matrix(data); nn = length(yy)
    set.seed(1202); tmpb.all = c.t.all = c.all = NULL
    for(rep in 1:20)
      { 
      	ind.t = sample(1:nn,rep=T); ind.v = setdiff(1:nn,ind.t)
        datat = data[ind.t,]; datav = data[ind.v,]
        beta = try(Est.ALASSO.GLM(datat,nopen.ind=nopen.ind,regularize=regularize,rtn="EST")[1:p0],silent=T)
        if(length(beta)>1){
            phatv = (cbind(1,as.matrix(datav[,-1]))%*%beta)
            c.fpr   = as.numeric(try(quantile(phatv[datav[,1]==0],1-fpr0),silent=T))
            c.all = cbind(c.all,c.fpr); ##c.t.all = c(c.t.all, c.fpr.t);  
         }
      }
    phat0 = c(cbind(1,as.matrix(data[,-1]))%*%betahat[1:p0])
    c.fpr.0   = quantile(phat0[data[,1]==0],1-fpr0)
    c.fpr.bc  = 0.368*c.fpr.0 + 0.632*apply(as.matrix(c.all),1,mean,na.rm=T)     
    c(g.logit(c.fpr.0),g.logit(c.fpr.bc))
  }

chat.FUN <- function(data,fpr0=0.05,betahat,regularize=T,nopen.ind=NULL)
  {
    yy = data[,1]; xx = as.matrix(data[,-1]); p0 = ncol(data)
    phat.0 = (cbind(1,as.matrix(data[,-1]))%*%betahat[1:p0])
    c.fpr.0   = quantile(phat.0[data[,1]==0],1-fpr0)
    set.seed(1202); nt = floor(nrow(data)/2); tmpb.all = c.t.all = c.all = NULL
    for(rep in 1:10)
      { 
        data = data[sample(1:nrow(data)),]
        for(k in 1:2)
          {
            tmpind = 1:nt+(k-1)*nt
            datat = data[-tmpind,]; datav = data[tmpind,]
        	beta = try(Est.ALASSO.GLM(datat,nopen.ind=nopen.ind,regularize=regularize,rtn="EST")[1:p0],silent=T)
            if(length(beta)>1){
                phatv = (cbind(1,as.matrix(datav[,-1]))%*%beta)
                phatt = (cbind(1,as.matrix(datat[,-1]))%*%beta)
                c.fpr   = quantile(phatv[datav[,1]==0],1-fpr0)
                c.fpr.t = quantile(phatt[datat[,1]==0],1-fpr0) 
                c.all = rbind(c.all,c.fpr); c.t.all = rbind(c.t.all, c.fpr.t)
                tmpb.all = cbind(tmpb.all,beta)}
          }
      }
    c.fpr.bc  = c.fpr.0 + apply(c.all-c.t.all,2,mean,trim=0.05)
    g.logit(c(c.fpr.0,c.fpr.bc))
  }

Predict.FUN <- function(data,newdata,fpr0=0.05,betahat)
  {
    yy = data[,1]; xx = as.matrix(data[,-1]); p0 = ncol(data)
    set.seed(1202); nt = floor(nrow(data)/3); tmpb.all = c.t.all = c.all = NULL
    for(rep in 1:10)
      { 
        data = data[sample(1:nrow(data)),]
        for(k in 1:3)
          {
            tmpind = 1:nt+(k-1)*nt
            datat = data[-tmpind,]; datav = data[tmpind,]
            beta = try(Est.ALASSO.GLM(datat,nopen.ind=c(1,2))[[1]][1:p0])
            if(length(beta)>1){
                phatv = g.logit(cbind(1,as.matrix(datav[,-1]))%*%beta)
                phatt = g.logit(cbind(1,as.matrix(datat[,-1]))%*%beta)
                c.fpr   = quantile(phatv[datav[,1]==0],1-fpr0)
                c.fpr.t = quantile(phatt[datat[,1]==0],1-fpr0) 
                c.all = c(c.all,c.fpr); c.t.all = c(c.t.all, c.fpr.t)
                tmpb.all = cbind(tmpb.all,beta)}
          }
      }
    #pnew.all = g.logit(cbind(1,as.matrix(newdata[,-1]))%*%tmpb.all)
    pnew.0 = g.logit(cbind(1,as.matrix(newdata[,-1]))%*%betahat[1:p0])
    phat.0 = g.logit(cbind(1,as.matrix(data[,-1]))%*%betahat[1:p0])
    c.fpr.0   = quantile(phat.0[data[,1]==0],1-fpr0)
    c.fpr.bc  = c.fpr.0 + mean(c.all-c.t.all)
    data.frame("patient_num"=newdata$patient,"Prob.RA"=round(pnew.0,5),
    "Cut-off.0"=round(c.fpr.0,5), "Cut-off.bc"=round(c.fpr.bc,5))
  }

ProbTestPos <- function(tpr,fpr,prev)
  {
    tpr*prev + fpr*(1-prev)
  }
      
    
L2Norm <- function(data,coef,intercept=F)
  {
    yy = data[,1]; xx = data[,-1]; if(intercept){xx=cbind(1,xx)}
    apply((yy - xx%*%t(coef))^2,2,mean)
  }

                        
PPV.FUN <- function(fpr,SE, mu0){ 1/(1+fpr/SE*(1-mu0)/mu0)}
NPV.Project <- function(npv.e0y0,p.e1,p.y0.e0,npv.e1=1)
  {
  	## P(D=0 | E=1 or E=0&Y=0) = {P(D=0|E=1)P(E=1)+P(D=0|E=0&Y=0)P(E=0&Y=0)}/{P(E=1)+P(E=0&Y=0)
  	npv.all = (npv.e1 * p.e1 + npv.e0y0*p.y0.e0*(1-p.e1))/(p.e1+p.y0.e0*(1-p.e1))
  	npv.all
  }

S.FUN <- function(yy,Yi,Di,yes.smooth=F)
  {
  	if(yes.smooth){
		Y1i = Yi[Di==1]; n1 = sum(Di); bw = bw.nrd(Y1i)/n1^0.6
		c(t(rep(1/n1,n1))%*%pnorm((Y1i-VTM(yy,n1))/bw))
  	}else{
		return((sum.I(yy,"<",Yi,Vi=Di)+sum.I(yy,"<=",Yi,Vi=Di))/sum(Di)/2)
  	}
    ##sum.I(yy,"<=",Yi,Vi=Di)/sum(Di)
  }

Sinv.FUN <- function(uu,Yi,Di,yes.smooth=F)
  {
    yy0<-unique(sort(Yi,decreasing=TRUE)); ss0 <- S.FUN(yy0,Yi,Di,yes.smooth=yes.smooth) 
    return(approx(ss0[!duplicated(ss0)],yy0[!duplicated(ss0)],uu,method="linear",f=0,rule=2)$y)
  }



CV.FUN <- function(data)
  {
    y = data[,1]; x = data[,-1]; nn=length(y)
    w = 1/abs(lm(y~x)$coef[-1]); x.w = x/VTM(w,nrow(x))
    s0.all = (1:100)/100
    fit.all = l1ce(y~x.w, standardize=F, bound=s0.all)
    K = 5; L2cv = NULL
    for(k in 1:K)
      {
        indv = 1:floor(nn/K) + (k-1)*floor(nn/K)
        indt = setdiff(1:nn,indv)
        fitk = l1ce(y~x.w,subset=indt, standardize=F, bound=s0.all)
        L2cv = rbind(L2cv, L2Norm(cbind(y,x.w)[indv,],coef(fitk),intercept=T))
      }
    L2cv = apply(L2cv,2,mean)
    s0 = min(s0.all[L2cv==min(L2cv)])
    bhat  = l1ce(y~x.w, standardize=F, bound=s0)
    list("b"=coef(bhat)/w, "s0"=s0, "lam0" =bhat$Lagrangian,"b0"=bhat$bound)    
  }

VTM<-function(vc, dm){
    matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

sum.I <- function(yy,FUN,Yi,Vi=NULL)
## sum_i I(yy FUN Yi)Vi
# Vi weight
  {
    if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
    # for each distinct ordered failure time t[j], number of Xi < t[j]
    pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')    
    if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos # number of Xi>= t[j]
    if (!is.null(Vi)) {
	   ## if FUN contains '=', tmpind is the order of decending
        if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
        ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
        Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
        return(rbind(0,Vi)[pos+1,])
    } else return(pos)
  }
convert <- function(fit) {
	rochat.auc = fit$rochat[1]
	rochat.values = matrix(fit$rochat[-1],ncol=6)
	colnames(rochat.values) = c("cutoff","pos.rate","FPR","TPR","PPV","NPV")
	roc.cv.auc = fit$roc.cv[1]
	roc.cv.values = matrix(fit$roc.cv[-1],ncol=6)
	colnames(roc.cv.values) = c("cutoff","pos.rate","FPR","TPR","PPV","NPV")
	betahat = fit$beta[1:(grep("\\<bini",names(fit$beta))[1]-1)]
	names(betahat) = gsub("\\<b\\.","",names(betahat))
	return(list(ROC.hat.auc=rochat.auc, ROC.hat.values=rochat.values, ROC.cv.auc=roc.cv.auc,ROC.cv.values=roc.cv.values, beta=betahat, Y.hat=fit$Si))
}

AUC = function(D,score,wgt=NULL){
	if(is.null(wgt)) {wgt=rep(1,length(D))}
	auc = sum(S.FUN(score,Yi=score,D*wgt,yes.smooth=F)*(1-D)*wgt)/sum((1-D)*wgt)
	return(auc)
}

ROC = function(D,score){
	roc = ROC.Est.FUN(D,score,0.5,seq(.01,.99,by=.01))
	roc = matrix(roc[-1],ncol=6)
	colnames(roc) = c("cutoff","est.pos.rate","FPR","TPR","PPV","NPV")
	return(roc)
}