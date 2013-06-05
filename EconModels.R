source("KTGP.R")
# ------------------------------------------------------------------------------
#  Load Data
# ------------------------------------------------------------------------------

GDP.csv <- read.csv("../Data/5206007_Income_From_GDP.csv")
GDP.raw <- GDP.csv[-c(1:11),c(1,79)]
colnames(GDP.raw) <- c("Date","GPD")
GDP <- as.xts(zoo(as.numeric(GDP.raw[,2]),as.yearmon(GDP.raw[,1],format="%b-%Y")))["1992::"]
plot(GDP)

CPI.csv <- read.csv("../Data/RBA_CPI.csv",stringsAsFactors=FALSE) 
CPI.raw <- CPI.csv[-c(1:14),1:2]
colnames(CPI.raw) <- c("Date","CPI")
idx <- which(CPI.raw$Date=="")
CPI <- as.xts(zoo(as.numeric(CPI.raw[-idx,2]),as.yearmon(CPI.raw[-idx,1],format="%b-%Y")))

# ------------------------------------------------------------------------------
#   1D GPTS models
# ------------------------------------------------------------------------------
# 1. GDP (% change)
# 1a.  SE model
X.train <- GDP["1992::2006"]
X.test <- GDP["2007::"]
n.ahead <- 20
a0 <- fit.SE.TS(coredata(X.train),par0=c(1,10,0.1),ub=c(10,50,1),scale=TRUE,OSU=FALSE)
GDP.gpts <- gausspr.SE.TS(y=coredata(X.train),kpar=a0$par,scale=TRUE)
GDP.fit <- predict.SE.TS(GDP.gpts,variance=TRUE)
GDP.fit.xts <- xts(GDP.fit$mean,index(X.train))
GDP.ahead <- predict.SE.TS(GDP.gpts,newdata=length(X.train):(length(X.train)+n.ahead),variance=TRUE)
time.index.ahead <- yearmon(last(index(X.train)) + (index(X.train[2])-index(X.train[1]))*(1:(n.ahead+1)))
GDP.ahead.xts <- xts(GDP.ahead$mean,time.index.ahead)
plot(rbind(X.train,X.test),main="GDP")
lines(GDP.fit.xts,col=4)
lines(GDP.ahead.xts,col=6)

# 1b.  SE + Periodic model
X.train <- GDP["1995::2009"]
X.test <- GDP["2010::"]
n.ahead <- 20
a0 <- fit.gpts(coredata(X.train),fn=fn.complex.TS.2,par0=c(1.8,45,0.5,45,0.5,4,0.2),lb=c(0.1,2,0.1,2,.2,.1,0.1),ub=c(5,45,5,45,1,10,1),scale=TRUE,OSU=FALSE)
GDP.gpts <- gausspr.complex.TS(kernel.complex.TS.2,coredata(X.train),a0$par,scale=TRUE)
GDP.fit <- predict.complex.TS(GDP.gpts,variance=TRUE)
GDP.fit.xts <- xts(GDP.fit$mean,index(X.train))
GDP.ahead <- predict.complex.TS(GDP.gpts,newdata=length(X.train):(length(X.train)+n.ahead),variance=TRUE)
time.index.ahead <- yearmon(last(index(X.train)) + (index(X.train[2])-index(X.train[1]))*(1:(n.ahead+1)))
GDP.ahead.xts <- xts(GDP.ahead$mean,time.index.ahead)
plot(rbind(X.train,X.test),main="GDP")
lines(GDP.fit.xts,col=4)
lines(GDP.ahead.xts,col=6)

plot(ROC(rbind(X.train,X.test),1)["2000::"],main="GDP")
lines(ROC(GDP.fit.xts,1),col=4)
GDP.ahead.roc <- ROC(rbind(last(GDP.fit.xts),GDP.ahead.xts),1)
lines(GDP.ahead.roc,col=6)

# 2. Inflation
# 2a.  SE model
X.train <- CPI["1970::2006"]
X.test <- CPI["2007::"]
n.ahead <- 20
a0 <- fit.SE.TS(coredata(X.train),par0=c(1,10,0.1),ub=c(10,50,1),scale=TRUE,OSU=FALSE)
CPI.gpts <- gausspr.SE.TS(y=coredata(X.train),kpar=a0$par,scale=TRUE)
CPI.fit <- predict.SE.TS(CPI.gpts,variance=TRUE)
CPI.fit.xts <- xts(CPI.fit$mean,index(X.train))
CPI.ahead <- predict.SE.TS(CPI.gpts,newdata=length(X.train):(length(X.train)+n.ahead),variance=TRUE)
time.index.ahead <- yearmon(last(index(X.train)) + (index(X.train[2])-index(X.train[1]))*(1:(n.ahead+1)))
CPI.ahead.xts <- xts(CPI.ahead$mean,time.index.ahead)
plot(rbind(X.train,X.test),main="Headline Inflation")
lines(CPI.fit.xts,col=4)
lines(CPI.ahead.xts,col=6)

# 2b.  SE + Periodic model
X.train <- CPI["1970::2006"]
X.test <- CPI["2007::"]
n.ahead <- 20
a0 <- fit.gpts(coredata(X.train),fn=fn.complex.TS.2,par0=c(1.8,15,0.5,5,.4,4,0.2),lb=c(0.1,2,0.1,2,.2,.1,0.1),ub=c(5,45,5,45,10,10,1),scale=TRUE,OSU=FALSE)
CPI.gpts <- gausspr.complex.TS(kernel.complex.TS.2,coredata(X.train),a0$par,scale=TRUE)
CPI.fit <- predict.complex.TS(CPI.gpts,variance=TRUE)
CPI.fit.xts <- xts(CPI.fit$mean,index(X.train))
CPI.ahead <- predict.complex.TS(CPI.gpts,newdata=length(X.train):(length(X.train)+n.ahead),variance=TRUE)
time.index.ahead <- yearmon(last(index(X.train)) + (index(X.train[2])-index(X.train[1]))*(1:(n.ahead+1)))
CPI.ahead.xts <- xts(CPI.ahead$mean,time.index.ahead)
plot(rbind(X.train,X.test),main="CPI")
lines(CPI.fit.xts,col=4)
lines(CPI.ahead.xts,col=6)

plot(ROC(rbind(X.train,X.test),1)["2000::"],main="CPI")
lines(ROC(CPI.fit.xts,1),col=4)
CPI.ahead.roc <- ROC(rbind(last(CPI.fit.xts),CPI.ahead.xts),1)
lines(CPI.ahead.roc,col=6)


# ------------------------------------------------------------------------------
#  AutoRegressive Gaussian Process Time Series Models
# ------------------------------------------------------------------------------

