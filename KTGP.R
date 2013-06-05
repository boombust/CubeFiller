
# -------------------------------------------------------------------
#  my TS kernels
# -------------------------------------------------------------------

kernel.SE <- function(x,y,sigma,l.sq,noise) {
  X <- outer(x,y,FUN="-")^2
  K <- sigma^2*exp(-0.5*X/l.sq^2)
  if (nrow(K) <= ncol(K))
    K <- K + cbind(diag(noise^2,nrow(K)),matrix(0,nrow=nrow(K),ncol=ncol(K)-nrow(K)))
  if (nrow(K) > ncol(K))
    K <- K + rbind(diag(noise^2,ncol(K)),matrix(0,nrow=nrow(K)-ncol(K),ncol=ncol(K)))
  return(K)
}

kernel.Periodic.Approx <- function(x,y,sigma,l.sq,omega,theta,noise) {
  X <- outer(x,y,FUN="-")^2
  K <- sigma^2*exp(-0.5*X/l.sq^2 - 2*sin(pi*omega*X)^2/theta^2)
#   K <- sigma^2*exp(-0.5*X/l.sq^2 - 2*sin(pi*X)^2/omega^2)
  if (nrow(K) <= ncol(K))
    K <- K + cbind(diag(noise^2,nrow(K)),matrix(0,nrow=nrow(K),ncol=ncol(K)-nrow(K)))
  if (nrow(K) > ncol(K))
    K <- K + rbind(diag(noise^2,ncol(K)),matrix(0,nrow=nrow(K)-ncol(K),ncol=ncol(K)))
  return(K)
}

kernel.Periodic <- function(x,y,sigma,omega,theta,noise) {
  X <- outer(x,y,FUN="-")^2
  K <- sigma^2*exp(-2*sin(pi*omega*X)^2/theta^2)
  #   K <- sigma^2*exp(-0.5*X/l.sq^2 - 2*sin(pi*X)^2/omega^2)
  if (nrow(K) <= ncol(K))
    K <- K + cbind(diag(noise^2,nrow(K)),matrix(0,nrow=nrow(K),ncol=ncol(K)-nrow(K)))
  if (nrow(K) > ncol(K))
    K <- K + rbind(diag(noise^2,ncol(K)),matrix(0,nrow=nrow(K)-ncol(K),ncol=ncol(K)))
  return(K)
}


kernel.RQ <- function(x,y,sigma,theta,l.sq,noise) {
  X <- outer(x,y,FUN="-")^2
#   K <- sigma^2*(1 + X/(2*theta*l.sq^2))^(-theta)
  logK <- 2*log(sigma) - theta*log(1 + X/(2*theta*l.sq^2))   # more stable??
  K <- exp(logK)
  if (nrow(K) <= ncol(K))
    K <- K + cbind(diag(noise^2,nrow(K)),matrix(0,nrow=nrow(K),ncol=ncol(K)-nrow(K)))
  if (nrow(K) > ncol(K))
    K <- K + rbind(diag(noise^2,ncol(K)),matrix(0,nrow=nrow(K)-ncol(K),ncol=ncol(K)))
  return(K)
}

kernel.complex.TS <- function(x,y,pars) {
  sigma <- pars[1]
  l.sq <- pars[2]
  K <- kernel.SE(x,y,sigma,l.sq,noise=0)
  sigma <- pars[3]
  omega <- pars[4]
  theta <- pars[5]
  K1 <- kernel.Periodic(x,y,sigma,omega,theta=1,noise=0)
  sigma <- pars[6]
  omega <- pars[7]
  theta <- pars[8]
  K2 <- kernel.Periodic(x,y,sigma,omega,theta=1,noise=0)
  K <- K + K1*K2
  noise = pars[length(pars)]
  if (nrow(K) <= ncol(K))
    K <- K + cbind(diag(noise^2,nrow(K)),matrix(0,nrow=nrow(K),ncol=ncol(K)-nrow(K)))
  if (nrow(K) > ncol(K))
    K <- K + rbind(diag(noise^2,ncol(K)),matrix(0,nrow=nrow(K)-ncol(K),ncol=ncol(K)))
  return(K)
}

kernel.complex.TS.2 <- function(x,y,pars) {
  sigma <- pars[1]
  l.sq <- pars[2]
  K <- kernel.SE(x,y,sigma,l.sq,noise=0)
  sigma <- pars[3]
  l.sq <- pars[4]
  omega <- pars[5]
  theta <- pars[6]
  K <- K + kernel.Periodic.Approx(x,y,sigma,l.sq,omega,theta=1,noise=0)
  noise = pars[length(pars)]
  if (nrow(K) <= ncol(K))
    K <- K + cbind(diag(noise^2,nrow(K)),matrix(0,nrow=nrow(K),ncol=ncol(K)-nrow(K)))
  if (nrow(K) > ncol(K))
    K <- K + rbind(diag(noise^2,ncol(K)),matrix(0,nrow=nrow(K)-ncol(K),ncol=ncol(K)))
  return(K)
}


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

# -------------------------------------------------------------------
#  log likelihood fitting of parameters
# -------------------------------------------------------------------

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

# -------------------------------------------------------------------
#   fitting the hyper-parameters
# -------------------------------------------------------------------
fn.SE.TS <- function(pars,yy,osu=FALSE) {
  sigma = pars[1]
  l.sq = pars[2]
  noise = pars[3]
  tt = c(0:(length(yy)-1))
  K <- kernel.SE(x=tt,y=tt,sigma,l.sq,noise)
  if (osu) { nll <- NLL.OSU(K,yy) }
  else {nll <- NLL(K,yy)}
  #   cat(nll,pars,"\n")
  return(nll)
}


fit.SE.TS <- function(y,par0=c(1,5,1),lb=c(0.1,1,0.1),ub=c(3,20,1),scale=TRUE,OSU=FALSE) {
  if (is.xts(y)) y <- coredata(y)
  if (scale) y <- scale(y)
  pars = par0
  ctrl <- list(trace=1)
  min.nll <- nlminb(pars,fn.SE.TS,gr=NULL,hessian=NULL,y,OSU,control=ctrl,lower=lb,upper=ub)  
  return(min.nll)
}

fn.complex.TS <- function(pars,yy,osu=FALSE) {
  tt <- c(0:(length(yy)-1))
  K <- kernel.complex.TS(x=tt,y=tt,pars)
  if (osu) { nll <- NLL.OSU(K,yy) }
  else {nll <- NLL(K,yy)}
  return(nll)
}

fn.complex.TS.2 <- function(pars,yy,osu=FALSE) {
  tt <- c(0:(length(yy)-1))
  K <- kernel.complex.TS.2(x=tt,y=tt,pars)
  if (osu) { nll <- NLL.OSU(K,yy) }
  else {nll <- NLL(K,yy)}
  return(nll)
}


fit.gpts <- function(y,fn=fn.complex.TS,par0=rep(1,length(y)),lb=0.1*par0,ub=5*par0,scale=TRUE,OSU=FALSE) {
  if (is.xts(y)) y <- coredata(y)
  if (scale) y <- scale(y)
  pars = par0
  ctrl <- list(trace=1)
  min.nll <- nlminb(pars,fn,gr=NULL,hessian=NULL,y,OSU,control=ctrl,lower=lb,upper=ub)  
  return(min.nll)
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

# -------------------------------------------------------------------
#  gaussian process fitting and prediction
# -------------------------------------------------------------------


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
    if (variance) sigma.star <- sigma^2+noise^2
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

# ----

gausspr.SE <- function(x,y,kpar,scale=TRUE) {
  # solve gpts on data y = f(x)
  stopifnot(length(x)==length(y))
  if (scale) y <- scale(y)
  n <- length(y)
  K <- kernel.SE(x,x,sigma=kpar[1],l.sq=kpar[2],noise=kpar[3])
  K.inv <- solve(K)
  beta <- K.inv %*% y
  return(list("K.inv"=K.inv,"beta"=beta,"kpar"=kpar,"y"=y,"scale"=scale))
}

gausspr.SE.TS <- function(y,kpar,scale=TRUE) {
  # solve gpts on data y = f(t)
  x.t <- c(0:(length(y)-1))
  gausspr.SE(x.t,y,kpar,scale)
}


predict.SE.TS <- function(gfit,newdata,variance=FALSE) {
  # re-use kernel parameters 
  sigma = gfit$kpar[1]
  l.sq = gfit$kpar[2]
  noise = gfit$kpar[3]
  n = length(gfit$y)
  sigma.star <- NULL
  # if no newdata, return evaluation on fitted data
  tt <- c(0:(n-1))
  if (missing(newdata)) {
    K <- kernel.SE(x=tt,y=tt,sigma=sigma,l.sq=l.sq,noise=0)
    mu.star <- K %*% gfit$beta  
    if (variance) sigma.star <- rep(sigma^2+noise^2,n)
  } else {
    K.star <- kernel.SE(x=tt,y=newdata,sigma=sigma,l.sq=l.sq,noise=0)
    mu.star <- t(K.star) %*% gfit$beta
    if (variance) {
      K.star.star <- kernel.SE(x=newdata,y=newdata,sigma=sigma,l.sq=l.sq,noise=0)
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


# ----

gausspr.complex <- function(kernel=kernel.complex.TS,x,y,kpar,scale=TRUE) {
  # solve gpts on data y = f(x)
  stopifnot(length(x)==length(y))
  if (scale) y <- scale(y)
  n <- length(y)
  K <- kernel(x,x,kpar)
  K.inv <- solve(K)
  beta <- K.inv %*% y
  return(list("K.inv"=K.inv,"beta"=beta,"kpar"=kpar,"y"=y,"x"=x,"scale"=scale,"kernel"=kernel))
}

gausspr.complex.TS <- function(kernel=kernel.complex.TS,y,kpar,scale=TRUE) {
  # solve gpts on data y = f(t)
  x.t <- c(0:(length(y)-1))
  gausspr.complex(kernel,x.t,y,kpar,scale)
}


predict.complex.TS <- function(gfit,newdata,variance=FALSE) {
  # re-use kernel parameters 
  pars <- gfit$kpar
  n <- length(gfit$x)
  tt <- c(0:(n-1))
  sigma.star <- NULL
  if (missing(newdata)) {
    # if no newdata, return evaluation on fitted data
    K <- gfit$kernel(x=tt,y=tt,pars)
    mu.star <- K %*% gfit$beta  
    if (variance) sigma.star <- rep(sum(pars[c(1,3,6,9)]),n)
  } else {
    # don't need to rescale newdata for a timeseries only model
    K.star <- gfit$kernel(x=tt,y=newdata,pars)
    mu.star <- t(K.star) %*% gfit$beta
    if (variance) {
      K.star.star <- gfit$kernel(x=newdata,y=newdata,pars)
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
