options(stringsAsFactors=F)
myread.csv = function(fn,...){read.csv(fn,stringsAsFactors =F,...)}
Find.String = function(vec,string)
  {
  	junk = strsplit(vec,string,fixed=T); 
  	lapply(junk,mypaste) != vec
  }

Find.String2 = function(vec,string)
  {
  	tmpf = function(string,vec){junk = strsplit(vec,string,fixed=T); (1:length(vec))[lapply(junk,mypaste) != vec]}
	unlist(lapply(string,tmpf,vec=vec))
  }


MOD <- function(aa, a0){aa - a0*floor(aa/a0)}

RUN.File.FUN <- function(newseed,oldseed,newfn,oldfn,newp=NULL,oldp="5",newrhostep=NULL,oldrhostep="1")
  {
    ## replace myseed=oldseed in fn0 by myseed=newseed 
    ## if newseed=315; fn0 should be "XX.R" and the new fn should be "XX315.R"##
    k0 = nchar(oldseed)
    junk=scan(file=oldfn,what="character",sep="\n")
    tmpind=substring(junk,1,6)=="myseed"
    substr(junk[tmpind],8,7+k0)<-as.character(newseed)
    if(!is.null(newp)){
        tmpind=substring(junk,1,4)=="beta"; k0=nchar(junk[tmpind])
        substr(junk[tmpind],k0-1,k0-1)<-as.character(newp)}   
    if(!is.null(newrhostep)){
        tmpind=substring(junk,1,3)=="rho"
        junk[tmpind]=gsub(paste("by=",oldrhostep,sep=""),paste("by=",newrhostep,sep=""),junk[tmpind])}
    newfn = paste(newfn,newseed,".R",sep="")
    write(junk,file=newfn)
  }

RUN.File.FUN2 <- function(newseed,oldseed,newfn=NULL,newfn.full=NULL,oldfn,str0="myseed")
  {
    ## replace myseed=oldseed in fn0 by myseed=newseed 
    ## if newseed=315; fn0 should be "XX.R" and the new fn should be "XX315.R"##
    k0 = nchar(oldseed); k0new = nchar(newseed)
    junk=scan(file=oldfn,what="character",sep="\n")
    k1 = nchar(str0); tmpind=substring(junk,1,k1)==str0
    junk[tmpind]=mypaste(substr(junk[tmpind],1,k1+1),as.character(newseed))
    if(is.null(newfn.full)){
      if(!is.null(newfn)){newfn = paste(newfn,newseed,".R",sep="")}else{newfn=oldfn}
    }else{
      newfn = newfn.full}
    write(junk,file=newfn)
  }


Trace <- function(mat){sum(diag(mat))}

Matrix.TableLatex2 <- function(mat,fn="result.txt",...)
  {
	nm.Row = row.names(mat); nm.col = names(mat); 
    mat = as.matrix(mat); nc = ncol(mat); mode(mat)="character"; nl = max(nchar(mat))
    tmpind = substring(mat,nl,nl) == "" ; 
    while(sum(tmpind)>0)
      {
      	mat[tmpind]=paste(mat[tmpind],"0",sep="")
    	tmpind = substring(mat,nl,nl) == "" 
      }
    tmpind = substring(mat,1,1)   == "0"; mat[tmpind] = substring(mat,2,nl)[tmpind]
    mat = matrix(mat,ncol=nc)
    mat = rbind(nm.col,mat)
    mat = paste(c("",nm.Row),apply(mat,1,paste,collapse="&"),sep="&")
    mat = paste(mat,"\\\\")  
    mat[length(mat)] = paste(mat[length(mat)],"\\hline")
    write.table(mat,fn,row.names=F,col.names=F,quote=F,...)
  }

Matrix.TableLatex <- function(mat, fn="result.txt",...)
  {
    nm.Row = row.names(mat); nm.col = names(mat); nc = ncol(mat); nr = nrow(mat)
    for(i in 1:nr){
      for(j in 1:nc){
        if(j < nc){mat[i,j] = paste(mat[i,j],"&")}
        if(j==1){mat[i,j]=paste(nm.Row[i],"&",mat[i,j])}
        if(j==nc){
            mat[i,j]=paste(mat[i,j],paste("\\", "\\",sep=""))
            if(i==nr){
              mat[i,j]=paste(mat[i,j],"\\hline")}}
        }
    }
    nm.col[nc] = paste(nm.col[nc],paste("\\", "\\",sep=""))
    mat = rbind(paste("&",nm.col),mat)
    write.table(mat,fn,row.names=F,col.names=F,quote=F,...)
   }

N.L.E <- function(yy, YY)  ## sum I(YY <= yy[i])
  {
    rank(c(yy+1e-8,YY))[1:length(yy)] - rank(yy)
  }

N.G.E <- function(yy, YY)  ## sum I(YY >= yy[i])
  {
    length(YY)-(rank(c(yy-1e-8,YY))[1:length(yy)] - rank(yy))
  }
lines.step <- function(x, y, left.cont=T, ...)
  {
    nx <- length(x)
    if(left.cont == T){
        x <- c(rep(x[-nx], rep(2, nx-1)), x[nx])
        y <- c(y[1], rep(y[-1], rep(2, nx-1)))
    }else{
        x <- c(x[1],rep(x[-1], rep(2,nx-1)))
        y <- c(rep(y[-nx],rep(2,nx-1)),y[nx]) }
    lines(x, y, ...)
  }

shade <- function(x, yl, yu, ...)
  {
    nl <- length(x)
    for(i in 1:(nl-1))
      {
        polygon(c(x[i], x[i], x[i+1], x[i+1]),
                c(yl[i], yu[i], yu[i+1], yl[i+1]), ...)
      }
  }

mypaste <- function(...)
  {
    paste(..., sep="")
  }
seq.gen <- function(cluster.size)
  {
    K <- max(cluster.size)
    tmp <- matrix(1:K-1, ncol=length(cluster.size),nrow=K)
    tmp[row(tmp) <= VTM(cluster.size, K)]
  }

mypaste <- function(...)
  {
    paste(..., sep="")
  }

## *********************************** ##
##    Generate step function plots     ##
## *********************************** ##

plot.step <- function(x, y, ...)
  {
    nx <- length(x)
    x <- c(rep(x[-nx], rep(2, nx-1)), x[nx])
    y <- c(y[1], rep(y[-1], rep(2, nx-1)))
    plot(x, y, ...)
  }

##################### Splus function : cumsum2 #####################
## cumsum2 gets cumsum of each column, that is, apply(X, 2, cumsum)

cumsum2 <- function(mydat)     #cumsum by row, col remains the same
  {
    if(is.null(dim(mydat))) return(cumsum(mydat))
    else{
      out <- matrix(cumsum(mydat), nrow=nrow(mydat))
      out <- out - VTM(c(0, out[nrow(mydat), -ncol(mydat)]), nrow(mydat))
      return(out)
    }
  }

#################### End of function : cumsum2 #####################
sum.I <- function(yy,FUN,Yi,Vi=NULL)
  {
    if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
    pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')
    if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos
    if (!is.null(Vi)) {
        if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
        ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
        Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
        return(rbind(0,Vi)[pos+1,])
    } else return(pos)
  }

sum.I.old <- function(yy,FUN,Yi,Vi=NULL)   ## sum I(yy FUN Yi)Vi
  {  
    if(FUN=="<"|FUN=="<=") { yy <- -yy; Yi <- -Yi}
    if(substring(FUN,2,2)=="=") yy <- yy + 1e-8 else yy <- yy - 1e-8
    pos <- rank(c(yy,Yi))[1:length(yy)] - rank(yy)
    if(is.null(Vi)){return(pos)}else{
      Vi <- cumsum2(as.matrix(Vi)[order(Yi),,drop=F])
      out <- matrix(0, nrow=length(yy), ncol=dim(as.matrix(Vi))[2])
      out[pos!=0,] <- Vi[pos,]
      if(is.null(dim(Vi))) out <- c(out)
      return(out) ## n.y x p
    }
  } 

################ Splus function : NR #################

## Newton Raphson algorithm

NR<-function(para, fn, fh, epn=1e-7, ...)
    {
      xold<-para; dif<-1;
      library(Matrix)
      while(dif>epn)
        {
      xnew<-xold-c(solve.Matrix(fh(xold, ...))%*%fn(xold,...))
      dif<-sqrt(sum((xnew-xold)^2)); 
      xold<-xnew ;  
        }
      xnew
    }

################## End of function NR ###################


        
### Splus Function : simplex : minimiazation ###

simplex0 <- function(b, fun, del = 1, ftol = 9.9999999999999995e-07,
                    itmax = 1000, ...)
        {
        # minimization using the Nelder-Mead simplex method, based on num-rec amoeba
        # output gives the minimizing parameter value (b), the minimum value
        # of fun, and comp giving the number of function evaluations and an 
        # error code (1=did not converge)
        # b=initial values, fun=function to be minimized, called as fn(b,...)
        # del=the step used from b in each coordinate direction to set up the
        # initial simplex, ftol=convergence criterion--max relative difference
        # between highest and lowest point in the simplex
        # setup
    iter <- 0
    np <- length(b)
    p <- matrix(b, np, np); browser()
    dim(p) <- c(np, np)
    p <- t(p)
    diag(p) <- diag(p) + del
    p <- rbind(b, p)
    y <- rep(0, np + 1)
    for(i in 1:(np + 1))
        y[i] <- fun(p[i, ,drop=F ], ...)    #
    psum <- apply(p, 2, sum)
    while(iter <= itmax) {
        o <- order(y)   
    # don't need a full sort, only smallest and two largest
        p <- p[o, ,drop=F ] # so could be done more efficiently
        y <- y[o]
        ilo <- 1
        ihi <- np + 1
        inhi <- np
        rtol <- (2 * abs(y[ihi] - y[ilo]))/(abs(y[ihi]) + abs(y[ilo]))
        if(rtol < ftol)
            return(list(b = p[ilo,  ], value = y[ilo], comp = c(
                iter = iter, error = 0)))
        if(iter >= itmax)
            return(list(b = p[ilo,  ], value = y[ilo], comp = c(
                iter = iter, error = 1)))
        iter <- iter + 2    
    # new point chosen by reflecting the worst current through the plane
              # of the others
        z <- smptry0(p, y, psum, fun, ihi, -1, ...)
        if(z[[1]] <= y[ilo]) {
        # new point is best--try going further
            z <- smptry0(z[[4]], z[[2]], z[[3]], fun, ihi, 2, ...)
            y <- z[[2]]
            psum <- z[[3]]
            p <- z[[4]]
        }
        else if(z[[1]] >= y[inhi]) {
            ysave <- z[[2]][ihi]    
    #new point is still worst, try smaller step
            z <- smptry0(z[[4]], z[[2]], z[[3]], fun, ihi, 0.5, ...)
            y <- z[[2]]
            psum <- z[[3]]
            p <- z[[4]]
            if(z[[1]] >= ysave) {
        # still bad, shrink simplex
                for(i in (1:(np + 1))[ - ilo]) {
                  psum <- (p[i,  ] + p[ilo,  ])/2
                  p[i,  ] <- psum
                  y[i] <- fun(psum, ...)
                }
                iter <- iter + np
                psum <- apply(p, 2, sum)
            }
        }
        else {
            y <- z[[2]]
            psum <- z[[3]]
            p <- z[[4]]
            iter <- iter - 1
        }
    }
}

########## End of Splus Function : simplex ###############

### Splus Function : smptry0 : a subfuntion for simplex ###
smptry0 <- function(p, y, psum, fun, ihi, fac, ...)
        {
          ndim <- ncol(p); 
          fac1 <- (1 - fac)/ndim
          fac2 <- fac1 - fac
          ptry <- psum * fac1 - p[ihi, ] * fac2
          ytry <- fun(ptry, ...)
          if(ytry < y[ihi]) {
            y[ihi] <- ytry
            psum <- psum - p[ihi,  ] + ptry
            p[ihi,  ] <- ptry
          }
          list(ytry, y, psum, p)
        }



SUMMARY.PAR <- function(output,truepar,yes.na.rm=T)
  {
    pn <- dim(output)[2]/2; mm <- dim(output)[1]
    mean <- apply(output[,1:pn,drop=F],2,mean,na.rm=yes.na.rm)
    emp.se <- sqrt(diag(var(output[,1:pn,drop=F],na.rm=yes.na.rm)))
    est.se <- apply(output[,-(1:pn),drop=F],2,mean,na.rm=yes.na.rm)
    covp <-apply(abs(output[,1:pn,drop=F]-VTM(truepar,mm))/
    output[,-(1:pn)]<=1.96, 2, mean,na.rm=yes.na.rm)
    out <- cbind("Bias"=mean-truepar,"EMP.SE"=emp.se, "EST.SE"=est.se,
                 "CovProb"=covp)

    round(out, 5)
  }




### generate 0,1,2,0,1,2,0,1,2,0,1,2,3...
### if the cluster sizes are 3, 3, 3, 4


EMP.ROC <- function(fpr, yd1, yd0)
  {
    yd0 <- sort(yd0); 
    nd0 <- length(yd0); nd1 <- length(yd1)
    ##rep(1/nd1, nd1) %*% ( 1*(yd1 > VTM(yd0, nd1)) + (yd1 == VTM(yd0, nd1))/2)
    c(rep(1/nd1, nd1) %*% (yd1 >= VTM(yd0[ceiling(nd0*(1-fpr))], nd1)))
  }
ROC.Binorm <- function(u, theta = theta0.true)
  {
    a <- theta[1]
    b <- theta[2]
    pnorm(a + b * qnorm(u))
  }


          

############  Splus Function VTM  #############
# change an n by 1 vector x to an m by n matrix
# with each row equal to x 

# Last changed: 25/05/00

VTM<-function(vc, dm){
     matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
    }



UNI<-function(vec)
        {
          new.ind<-duplicated(vec)
          repeat{
            vec<-vec+0.000001*new.ind
            new.ind<-duplicated(vec)
            if(sum(new.ind)==0) break
          }
      vec
        }


SUM.bycluster <- function(mymat, mycluster)
  {
    ord.cluster <- order(mycluster)
    mycluster <- sort(mycluster)
    clustersize <- table(mycluster)
    ind.tmp <- cumsum(clustersize)

    if(is.null(dim(mymat)))
      {
        mymat <- mymat[ord.cluster]
        aout <- cumsum(mymat)[ind.tmp]
        aout <- aout - c(0, aout[-length(ind.tmp)])
      }
    else
      {
        mymat <- mymat[ord.cluster,]
        aout <- as.matrix(cumsum2(mymat)[ind.tmp,])
        aout <- aout - rbind(0, aout[-length(ind.tmp),])
      }
    aout
  }




SUMMARY.PAR <- function(output,truepar, novar=T)
  {
    pn <- length(truepar); mm <- dim(output)[1]; print(mm)
    mean <- apply(output[,1:pn,drop=F],2,mean)
    emp.se <- sqrt(diag(var(output[,1:pn,drop=F])))
    est.se <- apply(output[,-(1:pn),drop=F],2,mean)
    covp <- apply(abs(output[,1:pn,drop=F]-VTM(truepar,mm))/
                  output[,-(1:pn)]<=1.96, 2, mean)
    out <- cbind("Bias"=mean-truepar,"EMP.SE"=emp.se, "EST.SE"=est.se,
                 "CovProb"=covp)
    out
  }
                   
      

simplex <- function(b,del=1,ftol=1e-6,itmax=1000,y0,delta0,x0,g0, fn) {
  iter <- 0
  np <- length(b)
  p0 <- rep(b,np)
  dim(p0) <- c(np,np)
  p0 <- t(p0)
  diag(p0) <- diag(p0)+del
  p0 <- rbind(b,p0)
  y <- rep(0,np+1)
  for (i in 1:(np+1)) y[i] <- fn(p0[i,],y0,delta0,x0,g0)
  ##
  psum <- apply(p0,2,sum)
  while (iter <= itmax) {
    o <- order(y) ## don't need a full sort, only smallest and two largest
    p0 <- p0[o,,drop=F]    ## so could be done more efficiently
    y <- y[o]
    ilo <- 1
    ihi <- np+1
    inhi <- np
    rtol <- 2*abs(y[ihi]-y[ilo])/(abs(y[ihi])+abs(y[ilo]))
    if (rtol < ftol) return(list(b=p0[ilo,],value=y[ilo],comp=c(iter=iter,error=0)))
    if (iter >= itmax) return(list(b=p0[ilo,],value=y[ilo],comp=c(iter=iter,error=1)))
    iter <- iter+2
    ## new point chosen by reflecting the worst current through the plane
    ## of the others
    z <- smptry(p0,y,psum,ihi,-1,y0,delta0,x0,g0)
    if (z[[1]] <= y[ilo]) { # new point is best--try going further
      z <- smptry(z[[4]],z[[2]],z[[3]],ihi,2,y0,delta0,x0,g0)
      y <- z[[2]]; psum <- z[[3]]; p0 <- z[[4]]
    } else if (z[[1]] >= y[inhi]) {
      ysave <- z[[2]][ihi] #new point is still worst, try smaller step
      z <- smptry(z[[4]],z[[2]],z[[3]],ihi,0.5,y0,delta0,x0,g0)
      y <- z[[2]]; psum <- z[[3]]; p0 <- z[[4]]
      if (z[[1]] >= ysave) { # still bad, shrink simplex
    for (i in (1:(np+1))[-ilo]) {
      psum <- (p0[i,]+p0[ilo,])/2
      p0[i,] <- psum
      y[i] <- fn(psum,y0,delta0,x0,g0)
    }
    iter <- iter+np
    psum <- apply(p0,2,sum)
      }
    } else {
      y <- z[[2]]; psum <- z[[3]]; p0 <- z[[4]]
      iter <- iter-1
    }
  }
}

### Splus Function : smptry : a subfuntion for simplex ###
smptry <- function(p, y, psum, fun, ihi, fac, ...)
        {
          ndim <- ncol(p)
          fac1 <- (1 - fac)/ndim
          fac2 <- fac1 - fac
          ptry <- psum * fac1 - p[ihi,  ] * fac2
          ytry <- fun(ptry, ...)
          if(ytry < y[ihi]) {
            y[ihi] <- ytry
            psum <- psum - p[ihi,  ] + ptry
            p[ihi,  ] <- ptry
          }
          list(ytry, y, psum, p)
        }

### Splus Function : simplex : minimiazation ###

simplex <- function(b, fun, del = 1, ftol = 1e-07,
                    itmax = 1000, ...)
  {
    ## minimization using the Nelder-Mead simplex method, based on
    ## num-rec amoeba output gives the minimizing parameter value (b),
    ## the minimum value of fun, and comp giving the number of function 
    ## evaluations and an error code (1=did not converge)
    ## b=initial values, fun=function to be minimized, called as fn(b,...)
    ## del=the step used from b in each coordinate direction to set up the
    ## initial simplex, ftol=convergence criterion--max relative difference
    ## between highest and lowest point in the simplex setup

    iter <- 0
    np <- length(b)
    p <- matrix(rep(b,np), nrow=np, ncol=np)
    dim(p) <- c(np, np)
    p <- t(p)
    diag(p) <- diag(p) + del
    p <- rbind(b, p)
    y <- rep(0, np + 1)
    for(i in 1:(np + 1))
      y[i] <- fun(p[i,  ], ...) ##
    psum <- apply(p, 2, sum)
    while(iter <= itmax) {
      o <- order(y) 
      ## don't need a full sort, only smallest and two largest
      p <- p[o, ,drop=F ]   ## so could be done more efficiently
      y <- y[o]
      ilo <- 1
      ihi <- np + 1
      inhi <- np
      rtol <- (2 * abs(y[ihi] - y[ilo]))/(abs(y[ihi]) + abs(y[ilo]))
      if(rtol < ftol)
        return(list(b = p[ilo,  ,drop=F  ], value = y[ilo], comp = c(
                                                     iter = iter, error = 0)))
      if(iter >= itmax)
        return(list(b = p[ilo,  ,drop=F  ], value = y[ilo], comp = c(
                                                     iter = iter, error = 1)))
      iter <- iter + 2  
      ## new point chosen by reflecting the worst current through the plane
      ## of the others
      z <- smptry(p, y, psum, fun, ihi, -1, ...)
      if(z[[1]] <= y[ilo]) {
        ## new point is best--try going further
        z <- smptry(z[[4]], z[[2]], z[[3]], fun, ihi, 2, ...)
        y <- z[[2]]
        psum <- z[[3]]
        p <- z[[4]]
      }
      else if(z[[1]] >= y[inhi]) {
        ysave <- z[[2]][ihi]    
    ## new point is still worst, try smaller step
        z <- smptry(z[[4]], z[[2]], z[[3]], fun, ihi, 0.5, ...)
        y <- z[[2]]
        psum <- z[[3]]
        p <- z[[4]]
        if(z[[1]] >= ysave) {
          ## still bad, shrink simplex
          for(i in (1:(np + 1))[ - ilo]) {
            psum <- (p[i,  ] + p[ilo,  ])/2
            p[i,  ] <- psum
            y[i] <- fun(psum, ...)
          }
          iter <- iter + np
          psum <- apply(p, 2, sum)
        }
      }
      else {
        y <- z[[2]]
        psum <- z[[3]]
        p <- z[[4]]
        iter <- iter - 1
      }
    }
  }


    
neq <- function(b,score.jcob, gtol=1e-8,iter=50,stepf=.5,steptol=1e-8,...)
  {
    ## Modified Newton's method for solving nonlinear equations
    ## b=initial parameter values
    ## gn=function to calculate vector of nonlinear equations, called as gn(b,...)
    ## jn=function to calc Jacobian of the nonlinear equations,
    ##   called as jn(b,...)  should return a matrix
    ## gtol solution identified if max(abs(gn(b,...)))<gtol
    ## iter=max # iterations (input), 
    ## stepf=the fraction by which the step size is reduced in each step of
    ##   the backtracking algorithm
    ## steptol--if step size is smaller than this then algorithm has stalled
    ## ... additional arguments to gn 
    ## returns a list with components b=approx solution, f=||gn(b)||^2/2,
    ##   and comp=a vector with components named
    ##   iter, giving the number of iterations used, an error 
    ##   code (error=0 no errors, =1 if error in directional search, =2 if 
    ##   max iterations exceeded, error=3 if iteration has stalled at a 
    ##   point that is not a solution), and steph giving the number of 
    ##   times the step length was reduced
    n <- length(b);  error <- 0; steph <- 0
    g0 <- score.jcob(b,...); j <- g0[[2]]; g0 <- g0[[1]]; f0 <- sqrt(sum(g0^2))/2
    for (ll in 1:iter)
      {
        sc <- -c(solve(j,g0,tol=1e-50)); bn <- b+sc;     
        g1 <- score.jcob(bn,...);   j <- g1[[2]]; g1 <- g1[[1]]; f1 <- sqrt(sum(g1^2))/2   
        i <- 0; lam <- -2*f0
        while (is.na(f1) || f1>f0+(1e-4)*lam) {
          i <- i+1; steph <- steph+1; sc <- sc*stepf; lam <- lam*stepf;
          bn <- b+sc; g1 <- score.jcob(bn,...)[[1]]; f1 <- sqrt(sum(g1^2))/2
          if (i>20) return(list(b=b,f=f0,comp=c(iter=ll,error=1,steph=steph)))}
        if (max(abs(g1))<gtol) # if true, iteration converged
          return(list(b=bn,f=f1,comp=c(iter=ll,error=0,steph=steph)))
        if (max(abs(b-bn)/pmax(abs(b),1))<steptol)
          return(list(b=bn,f=f1,comp=c(iter=ll,error=3,steph=steph)))
        b <- bn; g0 <- g1; f0 <- f1
      }
    ## max number of iterations exceeded
    list(b=bn,f=f1,comp=c(iter=ll,error=2,steph=steph))
  }


    
neq0 <- function(b,score,jcob, gtol=1e-6,iter=50,stepf=.5,steptol=1e-8,...)
  {
    ## Modified Newton's method for solving nonlinear equations
    ## b=initial parameter values
    ## gn=function to calculate vector of nonlinear equations, called as gn(b,...)
    ## jn=function to calc Jacobian of the nonlinear equations,
    ##   called as jn(b,...)  should return a matrix
    ## gtol solution identified if max(abs(gn(b,...)))<gtol
    ## iter=max # iterations (input), 
    ## stepf=the fraction by which the step size is reduced in each step of
    ##   the backtracking algorithm
    ## steptol--if step size is smaller than this then algorithm has stalled
    ## ... additional arguments to gn 
    ## returns a list with components b=approx solution, f=||gn(b)||^2/2,
    ##   and comp=a vector with components named
    ##   iter, giving the number of iterations used, an error 
    ##   code (error=0 no errors, =1 if error in directional search, =2 if 
    ##   max iterations exceeded, error=3 if iteration has stalled at a 
    ##   point that is not a solution), and steph giving the number of 
    ##   times the step length was reduced
    n <- length(b);  error <- 0; steph <- 0
    g0 <- score(b,...); f0 <- sum(g0*g0)/2
    for (ll in 1:iter)
      {
        j <- jcob(b,...); sc <- -c(solve(j,g0,tol=1e-50))
        bn <- b+sc;     g1 <- score(bn,...); f1 <- sum(g1^2)/2   
        i <- 0; lam <- -2*f0
        while (is.na(f1) || f1>f0+(1e-4)*lam) {
          i <- i+1; steph <- steph+1; sc <- sc*stepf; lam <- lam*stepf;
          bn <- b+sc; g1 <- score(bn,...); f1 <- sum(g1^2)/2
          if (i>20) return(list(b=b,f=f0,comp=c(iter=ll,error=1,steph=steph)))}
        if (max(abs(g1))<gtol) # if true, iteration converged
          return(list(b=bn,f=f1,comp=c(iter=ll,error=0,steph=steph)))
        if (max(abs(b-bn)/pmax(abs(b),1))<steptol)
          return(list(b=bn,f=f1,comp=c(iter=ll,error=3,steph=steph)))
        b <- bn; g0 <- g1; f0 <- f1
      }
    ## max number of iterations exceeded
    list(b=bn,f=f1,comp=c(iter=ll,error=2,steph=steph))
  }
