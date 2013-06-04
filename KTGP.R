
kernel.SE.ARGPTS <- function(x,y,sigma,l.sq,noise) {
  nx <- dim(x)[1]
  tx <- c(0:(nx-1))/l.sq[1]
  ny <- dim(y)[1]
  ty <- c(0:(ny-1))/l.sq[1]
  X <- outer(tx,ty,FUN="-")^2
  for (i in 1:ncol(X)) {
    x.lag <- x[,i]/l.sq[i+1]
    y.lag <- y[,i]/l.sq[i+1]
    X <- X + outer(x.lag,y.lag,"-")^2
  }
  K <- sigma^2*exp(-0.5*X)
  if (nrow(K) <= ncol(K))
    K <- K + cbind(diag(noise^2,nrow(K)),matrix(0,nrow=nrow(K),ncol=ncol(K)-nrow(K)))
  if (nrow(K) > ncol(K))
    K <- K + rbind(diag(noise^2,ncol(K)),matrix(0,nrow=nrow(K)-ncol(K),ncol=ncol(K)))
  return(K)
}

kernel.SE.ARGPTS.pred <- function(x,y,sigma,l.sq) {
  # assumes y is newdata, ie its time index begins and index(x)+1
  nx <- dim(x)[1]
  tx <- c(0:(nx-1))/l.sq[1]
  ny <- dim(y)[1]
  ty <- c(nx:(nx+ny-1))/l.sq[1]
  X <- outer(tx,ty,FUN="-")^2
  for (i in 1:ncol(x)) {
    x.lag <- x[,i]/l.sq[i+1]
    y.lag <- y[,i]/l.sq[i+1]
    X <- X + outer(x.lag,y.lag,"-")^2
  }
  K <- sigma^2*exp(-0.5*X)
  return(K)
}


NLL <- function(K,y) {
  stopifnot(is.matrix(K))
  n <- dim(K)[1]
  alpha <- solve(K,y)
  detK <- det(K)
  logdetK <- log(detK)
  if (!is.finite(logdetK)) logdetK <- -1.e10 
  logP <- -0.5 * t(y) %*% alpha - 0.5 * logdetK - n/2*log(2*pi)
  #    cat(-1*logP,"|",length(idx),"\n")
  return(-1*logP)
}

fn.SE.ARGPTS <- function(pars,dat,osu) {
  n <- ncol(dat)
  x <- dat[,-n]
  y <- as.matrix(dat[,n],ncol=1)
  sigma = pars[1]
  l.sq = pars[2:(length(pars)+1)]
  noise = pars[length(pars)]
  K <- kernel.SE.ARGPTS(x,x,sigma,l.sq,noise)
  if (osu) { nll <- NLL.OSU(K,y) }
  else {nll <- NLL(K,y)}
#   cat(nll,pars,"\n")
  return(nll)
}


fit.SE.ARGPTS <- function(x,y,par0=c(1,0.5,3,3,3,.1),lb=0.1*par0,ub=5*par0,scale=TRUE,OSU=FALSE) {
  if (is.zoo(y)) y <- coredata(y)
  if (is.zoo(x)) x <- coredata(x)
  if (scale) {
    x <- scale(x)
    y <- scale(y)
  }
  dat <- as.matrix(cbind(x,y))
  pars = par0
  ctrl <- list(trace=1)
  min.nll <- nlminb(pars,fn.SE.ARGPTS,gr=NULL,hessian=NULL,dat,OSU,control=ctrl,lower=lb,upper=ub)  
  return(min.nll)
}



gausspr.ARGPTS <- function(x,y,kpar,scale=TRUE) {
  # solve gpts on data y = f(x)
  stopifnot(nrow(x)==nrow(y))
  if (scale) {
    y <- scale(y)
    x <- scale(x)
  }
  n <- length(y)
  K <- kernel.SE.ARGPTS(x,x,sigma=kpar[1],l.sq=kpar[2:(length(kpar)-1)],noise=kpar[length(kpar)])
  K.inv <- solve(K)
  beta <- K.inv %*% y
  return(list("K.inv"=K.inv,"beta"=beta,"kpar"=kpar,"y"=y,"scale"=scale,"x"=x))
}


predict.ARGPTS <- function(gfit,newdata,variance=FALSE) {
  # re-use kernel parameters 
  sigma = gfit$kpar[1]
  l.sq = gfit$kpar[2:(length(gfit$kpar)-1)]
  noise = gfit$kpar[length(gfit$kpar)]
  n = nrow(gfit$x)
  sigma.star <- NULL
  # if no newdata, return evaluation on fitted data
  if (missing(newdata)) {
    K <- kernel.SE.ARGPTS(gfit$x,gfit$x,sigma,l.sq,noise=0)
    mu.star <- K %*% gfit$beta  
    #     if (variance) sigma.star <- rep(sigma^2+noise^2,n)
    if (variance) sigma.star <- rep(sigma^2,n)
  } else {
    if(gfit$scale) newdata <- (newdata - attr(gfit$y,"scaled:center"))*attr(gfit$y,"scaled:scale")
    K.star <- kernel.SE.ARGPTS.pred(gfit$x,y=newdata,sigma,l.sq)
    mu.star <- t(K.star) %*% gfit$beta
    if (variance) {
      K.star.star <- kernel.SE.ARGPTS.pred(x=newdata,y=newdata,sigma,l.sq)
      sigma.star <- diag(K.star.star - t(K.star) %*% gfit$K.inv %*% K.star)
    }
  }
  # rescale the predictions to match original y data
  if (gfit$scale) {
    mu.star <- mu.star*attr(gfit$y,"scaled:scale") + attr(gfit$y,"scaled:center")
    if (variance) sigma.star <- sigma.star * attr(gfit$y,"scaled:scale")^2
  }
  return(list("mean"=mu.star,"var"=sigma.star))
}
