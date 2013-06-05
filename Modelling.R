library(quantmod)
# library(softImpute)
source("./KTGP.R")

# ------------------------------------------------------------------------------
#  load RData
# ------------------------------------------------------------------------------

load("GeographicData.RData")
load("TS_National.RData")
load("TS_WA.RData")

# ------------------------------------------------------------------------------
#   Make BBI
# ------------------------------------------------------------------------------

make.xts <- function(x) {
yrs <- as.yearmon(gsub("X","",colnames(x)),format="%Y")
xts(t(x),yrs)
}

Unemployment.State.xts <- make.xts(Unemployment.State)
TaxableIncome.State.xts <- make.xts(TaxableIncome.State)
MotorVehicles.State.xts <- make.xts(MotorVehicles.State)

BBI.State.xts <- scale(100-Unemployment.State.xts,center=TRUE)*scale(MotorVehicles.State.xts,center=FALSE)
plot.zoo(BBI.State.xts,main="BoomBustIndex by State")

# make xts of location data
MotorVehicles.Deltas.xts <- make.xts(MotorVehicles.Deltas[,-c(1:3)])
colnames(MotorVehicles.Deltas.xts) <- MotorVehicles.Deltas[,3]
Unemployment.Deltas.xts <- make.xts(Unemployment.csv[,-c(1:3)])
colnames(Unemployment.Deltas.xts) <- Unemployment.csv[,3]

BBI.Annual <- scale(100-Unemployment.Deltas.xts)*scale(MotorVehicles.Deltas.xts,center=FALSE)
write.csv(t(BBI.Annual),file="BBI_annual.csv")


# ------------------------------------------------------------------------------
#   Build State reference model
# ------------------------------------------------------------------------------

BBI.WA <- BBI.State.xts[,8]
X <- TS.WA[index(BBI.WA),]
XX <- TS.WA["2010-06-02::"]

a0 <- fit.SE.ARGPTS(coredata(X),coredata(BBI.WA),par0=c(1,rep(1,ncol(X)+1),1),lb=c(0.5,rep(1,ncol(X)+1),0),ub=c(5,rep(10,ncol(X)+1),2),scale=TRUE,OSU=FALSE)

gp1 <- gausspr.ARGPTS(x=coredata(X),y=coredata(BBI.WA),kpar=a0$par,scale=TRUE)
gp1.insample <- predict.ARGPTS(gp1,variance=TRUE)
gp1.insample.xts <- xts(gp1.insample$mean,index(BBI.WA))
gp1.insample.var <- cbind(gp1.insample.xts + sqrt(gp1.insample$var),gp1.insample.xts - sqrt(gp1.insample$var))
plot(BBI.WA,ylim=range(gp1.insample.var))
lines(gp1.insample.xts,col=4,type='b')
lines(gp1.insample.var[,1],col=3)
lines(gp1.insample.var[,2],col=3)
gp1.ahead <- predict.ARGPTS(gp1,newdata=coredata(XX),variance=TRUE)
gp1.ahead.xts <- xts(gp1.ahead$mean,index(XX))
plot(gp1.ahead.xts)

# ------------------------------------------------------------------------------
#  try softImpute to fill in the deltas
# ------------------------------------------------------------------------------

MotorVehicles.WA.xts <- make.xts(MotorVehicles.WA[,-c(1:3)])
colnames(MotorVehicles.WA.xts) <- MotorVehicles.WA[,3]

Motor.WA <- merge(TS.WA["2004::",8],MotorVehicles.WA.xts,fill=NA)

Motor.WA.imputed <- softImpute(coredata(Motor.WA),rank.max=100,lambda=30)

xna <- as(coredata(Motor.WA),"Incomplete")
fit2 <- softImpute(xna,rank=100,lambda=30)
Motor.WA.imputed <- fit2$u %*% diag(fit2$d) %*% t(fit2$v)
head(Motor.WA.imputed[,1:10],20)

# ------------------------------------------------------------------------------
#  model national quantities
# ------------------------------------------------------------------------------

# New Vehicle Sales
X <- TS.National["2003::2010",1]
X.oos <- TS.National["2011::",1]

par.names <- c("sigma","l.sq","sigma","omega","theta","sigma","omega","theta","noise")
par1 <- matrix(c(1,20,1,1,1,1,12,1,0.5),nrow=1)
colnames(par1) <- par.names
par.min <- matrix(c(0.1,1,0.5,1,0.5,0.5,12,.5,.2),nrow=1)
colnames(par.min) <- par.names
par.max <- matrix(c(10,100,10,1,100,10,12,100,5),nrow=1)
colnames(par.max) <- par.names
rbind(par1,par.min,par.max)

a0 <- fit.complex.TS(coredata(X),par0=par1,lb=par.min,ub=par.max,scale=TRUE,OSU=FALSE)
par.fit <- matrix(a0$par,nrow=1)
colnames(par.fit) <- par.names

gp.vehicles <- gausspr.complex.TS(y=coredata(X),kpar=par.fit,scale=TRUE)
vehic.fit <- predict.complex.TS(gp.vehicles,variance=TRUE)
vehic.xts <- xts(vehic.fit$mean,index(X))
vehic.var <- cbind(vehic.xts + sqrt(vehic.fit$var),vehic.xts - sqrt(vehic.fit$var))
plot(X)
lines(vehic.xts,col=4)
lines(vehic.var[,1],col=6)
lines(vehic.var[,2],col=6)

vehic.ahead <- predict.complex.TS(gp.vehicles,newdata=c(nrow(X):(nrow(X)+nrow(X.oos)-1)),variance=TRUE)
vehic.pred.xts <- xts(vehic.ahead$mean,index(X.oos))
vehic.pred.var <- cbind(vehic.pred.xts + 1.6*sqrt(vehic.ahead$var),vehic.pred.xts - 1.6*sqrt(vehic.ahead$var))
plot.xts(rbind(X,X.oos)["2006::"],main="New Vehicle Sales - National")
lines(vehic.pred.xts,col=4)
lines(vehic.pred.var[,1],col=3)
lines(vehic.pred.var[,2],col=3)


# Capex
X <- na.omit(TS.National["2003::2012",5])
X.oos <- na.omit(TS.National["2012-07::",5])

par.names <- c("sigma","l.sq","noise")
par1 <- matrix(c(1,20,0.5),nrow=1)
colnames(par1) <- par.names
par.min <- matrix(c(0.1,1,0.1),nrow=1)
colnames(par.min) <- par.names
par.max <- matrix(c(10,100,10),nrow=1)
colnames(par.max) <- par.names
rbind(par1,par.min,par.max)

a0 <- fit.SE.TS(coredata(X),par0=par1,lb=par.min,ub=par.max,scale=TRUE,OSU=FALSE)
par.fit <- matrix(a0$par,nrow=1)
colnames(par.fit) <- par.names

gp.capex <- gausspr.SE.TS(y=coredata(X),kpar=par.fit,scale=TRUE)
capex.fit <- predict.SE.TS(gp.capex,variance=TRUE)
capex.xts <- xts(capex.fit$mean,index(X))
capex.var <- cbind(capex.xts + sqrt(capex.fit$var),capex.xts - sqrt(capex.fit$var))
plot(X)
lines(capex.xts,col=4)
lines(capex.var[,1],col=6)
lines(capex.var[,2],col=6)

capex.ahead <- predict.SE.TS(gp.capex,newdata=c(nrow(X):(nrow(X)+max(nrow(X.oos),24)-1)),variance=TRUE)
time.index.ahead <- yearmon(last(index(X)) + (index(X[2])-index(X[1]))*(1:(max(nrow(X.oos),24))))
capex.pred.xts <- xts(capex.ahead$mean,time.index.ahead)
capex.pred.var <- cbind(capex.pred.xts + 1.64*sqrt(capex.ahead$var),capex.pred.xts - 1.64*sqrt(capex.ahead$var))
plot.xts(merge(rbind(X,X.oos)["2007::"],capex.pred.xts,all=TRUE),main="Total Capex - National")
lines(capex.pred.xts,col=4)
lines(capex.pred.var[,1],col=3)
lines(capex.pred.var[,2],col=3)


