toInstall <- c("fGarch","zoo")
lapply(toInstall, require, character.only=T)

############################# Load data ############################
#y <- read.table("h15Hist.tsv", header=TRUE, sep="\t")
y <- read.table("frwTSY.txt", header=TRUE, sep="\t")
rdate=as.Date(as.character(y$date),"%Y%m%d")
# for Time-Series Plotting
yts  <- zoo(y[,2:11],rdate,frequency=365)
# extract yields of each maturity
ym <- y[,2:11]

###################### Run AR(2)+GARCH(1,1) ########################
nforecast <- 5
# each column is forecast of each maturity bond yield
mean <- matrix(, nrow=nforecast, ncol=10)
vol <- matrix(, nrow=nforecast, ncol=10)
#for (i in 1:10){
for (i in 9){
  series <- ym[,i]
  series <- series[!is.na(series)] # remove NA values
  #dates <- rdate[!is.na(series)]
  series <- series[series!=0] # remove 0
  #dates <- rdate[series!=0]
  slog <- ts(log(series/100), start=1)
  sdiff <- diff(slog)
  # AR(2)+GARCH(1,1)
  re <- garchFit(formula=~arma(2,0)+garch(1,1),
                 data=sdiff, cond.dist=c("norm"), trace=FALSE)
  #pred <- predict(re, n.ahead = nforecast, plot=TRUE, nx=50, mse=c("cond"))
  pred <- predict(re, n.ahead = nforecast, plot=F, nx=50, mse=c("cond"))
  mean[,i] <- pred[,1]
  vol[,i] <- pred[,2]
}
#################### Converting Back to Yield ###################
yld <- matrix(, nrow=nforecast, ncol=10)
sum <- rep(0,10)
for (i in 1:nrow(mean)){
  sum <- sum + mean[i,]
  print(exp(sum + log(tail(ym, n=1))))
}
