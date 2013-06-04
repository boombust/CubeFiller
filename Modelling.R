library(quantmod)
library(softImpute)

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