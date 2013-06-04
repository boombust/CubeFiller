
# ------------------------------------------------------------------------------
#    load detailed geographic data
# ------------------------------------------------------------------------------

# PersonalIncome.csv <- read.csv("../ABS.Stat data/Estimates of Personal Income (year ended 30 June)- Average Total income (excl. Govt. pensions) ($).csv",header=TRUE)

Unemployment.csv <- read.csv("../ABS.Stat data/Estimates of Unemployment (June qtr)- Unemployment rate (Percent).csv")

TaxableIncome.csv <- read.csv("../ABS.Stat data/Taxation Statistics (year ended 30 June)- Taxable and non-taxable individuals- Avg taxable income ($).csv")

MotorVehicles.csv <- read.csv("../ABS.Stat data/Registered Motor Vehicles (at 31 March)- Total registered motor vehicles (Number).csv")

# ------------------------------------------------------------------------------
#  Make state totals and produce 
# ------------------------------------------------------------------------------

MotorVehicles.State <- apply(MotorVehicles.csv[,-c(1:3)],2,function(x) tapply(as.numeric(x),MotorVehicles.csv$State,FUN=sum))

MotorVehicles.Deltas <- MotorVehicles.csv
for (i in 1:nrow(MotorVehicles.csv)) {
  this.state <- MotorVehicles.csv[i,1]
  idx <- which(row.names(MotorVehicles.State)==this.state)
  state.total <- MotorVehicles.State[idx,-c(1:2)]
  MotorVehicles.Deltas[i,c(4:10)] <- MotorVehicles.csv[i,c(4:10)]/state.total    
}
write.csv(MotorVehicles.Deltas,file="MotorVehiclesDeltas.csv")

MotorVehicles.WA <- MotorVehicles.csv[MotorVehicles.csv$State=="WA",]
row.names(MotorVehicles.WA) <- NULL
# ---
# unemployment - already in percents so state "totals" are averages

Unemployment.State <- apply(Unemployment.csv[,-c(1:3)],2,function(x) tapply(as.numeric(x),Unemployment.csv$State,FUN=function(y) mean(na.omit(y))))

write.csv(Unemployment.csv,"UnemploymentDeltas.csv")

# ---
# Taxable income

TaxableIncome.State <- apply(TaxableIncome.csv[,-c(1:3)],2,function(x) tapply(as.numeric(x),TaxableIncome.csv$State,FUN=function(y) mean(na.omit(y))))


TaxableIncome.Deltas <- TaxableIncome.csv
for (i in 1:nrow(TaxableIncome.csv)) {
  this.state <- TaxableIncome.csv[i,1]
  idx <- which(row.names(TaxableIncome.State)==this.state)
  state.total <- TaxableIncome.State[idx,-c(1:2)]
  TaxableIncome.Deltas[i,c(4:7)] <- TaxableIncome.csv[i,c(4:7)]/state.total    
}

# write.csv(TaxableIncome.Deltas,"TaxableIncomeDeltas.csv")

# save(TaxableIncome.State,TaxableIncome.Deltas,Unemployment.State,Unemployment.csv,MotorVehicles.Deltas,MotorVehicles.State,file="GeographicData.RData")

# ------------------------------------------------------------------------------
#   load time series data
# ------------------------------------------------------------------------------
library(quantmod)

NewVehicleSales.csv <- read.csv("../Data/New Vehicle Sales - State.csv")
CPI.csv <- read.csv("../Data/CapitalCity CPI.csv")
colnames(CPI.csv)[10] <- "Australia"
WagePriceIndex.csv <- read.csv("../Data/WagePriceIndex - State.csv")
HousePrices.csv <- read.csv("../Data/HousePriceIndex.csv")
colnames(HousePrices.csv)[9] <- "Australia"

# Load Capex Data
Capex_ACT.csv <- read.csv("../Data/Capex_ACT.csv")
Capex_NT.csv <- read.csv("../Data/Capex_NT.csv")
Capex_Tas.csv <- read.csv("../Data/Capex_Tas.csv")
Capex_WA.csv <- read.csv("../Data/Capex_WA.csv")
Capex_SA.csv <- read.csv("../Data/Capex_SA.csv")
Capex_Qld.csv <- read.csv("../Data/Capex_Qld.csv")
Capex_Vic.csv <- read.csv("../Data/Capex_Vic.csv")
Capex_NSW.csv <- read.csv("../Data/Capex_NSW.csv")
Capex_National.csv <- read.csv("../Data/Capex_National.csv")

# turn into xts
NewVehicles.xts <- xts(NewVehicleSales.csv[,2:10],as.yearmon(NewVehicleSales.csv[,1],format="%b-%y"))
colnames(NewVehicles.xts) <- paste("NewVehicles",colnames(NewVehicles.xts),sep=".")
CPI.xts <- xts(CPI.csv[,2:10],as.yearmon(CPI.csv[,1],format="%b-%y"))
colnames(CPI.xts) <- paste("CPI",colnames(CPI.xts),sep=".")
WagePriceIndex.xts <- xts(WagePriceIndex.csv[,2:10],as.yearmon(WagePriceIndex.csv[,1],format="%b-%y"))
colnames(WagePriceIndex.xts) <- paste("WagePriceIndex",colnames(WagePriceIndex.xts),sep=".")
HousePrices.xts <- xts(HousePrices.csv[,2:9],as.yearmon(HousePrices.csv[,1],format="%b-%y"))
colnames(HousePrices.xts) <- paste("HousePrices",colnames(HousePrices.xts),sep=".")

Capex_ACT.xts <- xts(Capex_ACT.csv[,2],as.yearmon(Capex_ACT.csv[,1],format="%b-%y"))
colnames(Capex_ACT.xts) <- "Capex.ACT"
Capex_NT.xts <- xts(Capex_NT.csv[,2],as.yearmon(Capex_NT.csv[,1],format="%b-%y"))
colnames(Capex_NT.xts) <- "Capex.NT"
Capex_Tas.xts <- xts(Capex_Tas.csv[,2:5],as.yearmon(Capex_Tas.csv[,1],format="%b-%y"))
colnames(Capex_Tas.xts) <- paste("Capex","Tas",colnames(Capex_Tas.xts),sep=".")
Capex_WA.xts <- xts(Capex_WA.csv[,2:5],as.yearmon(Capex_WA.csv[,1],format="%b-%y"))
colnames(Capex_WA.xts) <- paste("Capex","WA",colnames(Capex_WA.xts),sep=".")
Capex_SA.xts <- xts(Capex_SA.csv[,2:5],as.yearmon(Capex_SA.csv[,1],format="%b-%y"))
colnames(Capex_SA.xts) <- paste("Capex","SA",colnames(Capex_SA.xts),sep=".")
Capex_Qld.xts <- xts(Capex_Qld.csv[,2:5],as.yearmon(Capex_Qld.csv[,1],format="%b-%y"))
colnames(Capex_Qld.xts) <- paste("Capex","Qld",colnames(Capex_Qld.xts),sep=".")
Capex_Vic.xts <- xts(Capex_Vic.csv[,2:5],as.yearmon(Capex_Vic.csv[,1],format="%b-%y"))
colnames(Capex_Vic.xts) <- paste("Capex","Vic",colnames(Capex_Vic.xts),sep=".")
Capex_NSW.xts <- xts(Capex_NSW.csv[,2:5],as.yearmon(Capex_NSW.csv[,1],format="%b-%y"))
colnames(Capex_NSW.xts) <- paste("Capex","NSW",colnames(Capex_NSW.xts),sep=".")
Capex_National.xts <- xts(Capex_National.csv[,2:5],as.yearmon(Capex_National.csv[,1],format="%b-%y"))
colnames(Capex_National.xts) <- paste("Capex","National",colnames(Capex_National.xts),sep=".")

# ---
# merge xts files and save as .RData
Capex.xts <- merge(Capex_ACT.xts,Capex_National.xts,Capex_NSW.xts,Capex_NT.xts,Capex_Qld.xts,Capex_SA.xts,Capex_Tas.xts,Capex_Vic.xts,Capex_WA.xts)
idx <- grep("Total",colnames(Capex.xts),value=FALSE)

TS.dat.quarterly <- merge(Capex.xts[,-idx],CPI.xts,WagePriceIndex.xts,HousePrices.xts)

TS.data <- merge(NewVehicles.xts,TS.dat.quarterly["1995::"])

TS.National <- na.approx(merge(NewVehicles.xts[,1],Capex_National.xts["1995::"],HousePrices.xts["1995::",8],WagePriceIndex.xts["1995::",9],CPI.xts["1995::",9],fill=NA))
save(TS.National,file="TS_National.RData")


TS.ACT <- na.approx(merge(NewVehicles.xts[,9],Capex_ACT.xts["1995::"],HousePrices.xts["1995::",7],WagePriceIndex.xts["1995::",8],CPI.xts["1995::",8]))
save(TS.ACT,file="TS_ACT.RData")

TS.WA <- na.approx(merge(NewVehicles.xts[,6],Capex_WA.xts["1995::"],HousePrices.xts[,4],WagePriceIndex.xts["1995::",5],CPI.xts["1995::",5]))
save(TS.WA,file="TS_WA.RData")

plot.zoo(TS.WA)