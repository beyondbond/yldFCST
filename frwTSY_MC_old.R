toInstall <- c("fGarch","zoo")
lapply(toInstall, require, character.only=T)
#################################### Utilities #################################### 
rateMC <- function(j, z, Nsim, Nforecast){
  data <- frwdiff[,j] # Use jth rate as benchmark
  # fit the log diff data to AR(2)/GARCH(1,1)
  fit <- garchFit(formula=~arma(2,0)+garch(1,1),
                  data=data, cond.dist=c("norm"), trace=FALSE)
  coef <- coef(fit) # coefficients of AR(2)/GARCH(1,1)
  n <- length(data)
  std <- tail(fit@sigma.t, n=1)
  unstd <- tail(fit@residuals, n=1) # error term of AR(2)
  # forecast matrix, ncol is the length of forecast, nrow is the number of simulations
  Fmatrix <- matrix(rep(0, Nsim*Nforecast), nrow=Nsim, ncol=Nforecast)
  # if i==1, AR(2) uses the last 2 observations
  # if i==2, AR(2) uses the last observation and the 1st forecast
  for (i in 1:Nforecast){
    zi <- z[,i]
    if (i==1){
      garch.std <- sqrt(coef[4]+coef[5]*unstd^2+coef[6]*std^2)
      garch.unstd <- garch.std*zi
      Fmatrix[,i] <- coef[1]+coef[2]*data[n]+coef[3]*data[n-1]+garch.unstd
    } else if (i==2){
      garch.std <- sqrt(coef[4]+coef[5]*(garch.unstd)^2+coef[6]*garch.std^2)
      garch.unstd <- garch.std*zi
      Fmatrix[,i] <- coef[1]+coef[2]*Fmatrix[,i-1]+coef[3]*data[n]+garch.unstd
    } else {
      garch.std <- sqrt(coef[4] + coef[5]*(garch.unstd)^2 + coef[6]*garch.std^2)
      garch.unstd <- garch.std*zi
      Fmatrix[,i] <- coef[1] + coef[2]*Fmatrix[,i-1] +coef[3]*Fmatrix[,i-2]+ garch.unstd
    }
  }
  return(Fmatrix)
}
# Monte Carlo simulation on rate with correlation coefficient with the benchmark rate
corMC <- function(j){
  data <- frwdiff[,j]
  fit <- garchFit(formula=~arma(2,0)+garch(1,1),
                  data=data, cond.dist=c("norm"), trace=FALSE)
  coef <- coef(fit)
  n <- length(data)
  std <- tail(fit@sigma.t, n=1)
  unstd <- tail(fit@residuals, n=1)
  rho <- cor(frwm[,c(j,9)])[1,2] #correlation
  Fmatrix <- matrix(rep(0, Nsim*Nforecast), nrow=Nsim, ncol=Nforecast)
  w <- matrix(rnorm(Nsim*Nforecast), nrow=Nsim, ncol=Nforecast)
  for (i in 1:Nforecast){
    zi <- z[,i]
    wi <- w[,i]
    if (i==1){
      garch.std <- sqrt(coef[4]+coef[5]*unstd^2+coef[6]*std^2)
      garch.unstd <- garch.std*rho*zi+garch.std*sqrt(1-rho^2)*wi
      Fmatrix[,i] <- coef[1]+coef[2]*data[n]+coef[3]*data[n-1]+garch.unstd
    } else if (i==2){
      garch.std <- sqrt(coef[4]+coef[5]*(garch.unstd)^2+coef[6]*garch.std^2)
      garch.unstd <- garch.std*rho*zi+garch.std*sqrt(1-rho^2)*wi
      Fmatrix[,i] <- coef[1]+coef[2]*Fmatrix[,i-1]+coef[3]*data[n]+garch.unstd
    } else {
      garch.std <- sqrt(coef[4] + coef[5]*(garch.unstd)^2 + coef[6]*garch.std^2)
      garch.unstd <- garch.std*rho*zi+garch.std*sqrt(1-rho^2)*wi
      Fmatrix[,i] <- coef[1] + coef[2]*Fmatrix[,i-1] +coef[3]*Fmatrix[,i-2]+ garch.unstd
    }
  }
  return(Fmatrix)
}
# convert log diff back to normal rates
diff2frw <- function(j, Fmatrix){
  n <- nrow(Fmatrix)
  m <- ncol(Fmatrix)
  frwFRCT <- matrix(rep(0, n*m), nrow=n, ncol=m)
  rowSum = rep(0,n)
  for (i in 1:m){
    rowSum <- rowSum + Fmatrix[,i]
    frwFRCT[,i] <- exp(rowSum + rep(log(tail(frwm[,j], n=1)), n))
  }
  
  return(frwFRCT)
}
##################################### Load data ####################################
path <- "Users/Jie/Google Drive/Work/"
filename <- "frwTSY.txt"
#frw <- read.table(paste(path,filename,sep=""), header=TRUE, sep="\t")
frw <- read.table(filename, header=TRUE, sep="\t")
frw[frw<0.1] <- 0.1 # replace rates < 0.1 with 0.1
frwm <- frw[,2:11]
frwdiff <- diff(ts(log(frwm/100), start=1))

############################## Monte Carlo simulation ##############################
################################# Benchmark rates ##################################
Nsim <- 1000 # number of simulations
Nforecast <- 2500 # how many days ahead to forecast
p <- 1 # use pth rate as benchmark
# generate Gaussian random numbers
z <- matrix(rnorm(Nsim*Nforecast), nrow=Nsim, ncol=Nforecast)
# forecast matrix of log diff
Fmatrix <- rateMC(p, z, Nsim, Nforecast)
# convert log diff to normal rates
frwFRCT <- diff2frw(p, Fmatrix)

########################## Other rates' forecast log diff ###########################
r <- c(1:10)
c <- r[!r==p] # rates excluding the benchmark rate
colnames <- colnames(frwm)
#for (j in c){
#  Fmatrix <- corMC(j)
#  filepath <- "/Users/Jie/Google Drive/Work/Output/"
#  write.table(Fmatrix, paste(filepath,"Fmatrix",colnames[j],".txt",sep=""),sep="\t")
#}

################################## Save quantiles ##################################
# Benchmark rate quantiles
mean <- apply(frwFRCT, 2, mean)
pct95 <- apply(frwFRCT,2,quantile,0.95)
pct05 <- apply(frwFRCT,2,quantile,0.05)
median <- apply(frwFRCT,2,quantile,0.5)
for (j in c){
  #Fmatrix <- read.table(paste(filepath,"Fmatrix",colnames[j],".txt",sep=""),sep="\t")
  Fmatrix <- corMC(j)
  frwFRCTj <- diff2frw(j, Fmatrix)
  mean <- cbind(mean, apply(frwFRCTj,2,mean))
  pct95 <- cbind(pct95, apply(frwFRCTj,2,quantile,0.95))
  pct05 <- cbind(pct05, apply(frwFRCTj,2,quantile,0.05))
  median <- cbind(median, apply(frwFRCTj,2,quantile,0.5))
}
colnames(mean) <- colnames(frwm)[c(9,1:8,10)]
colnames(pct95) <- colnames(frwm)[c(9,1:8,10)]
colnames(pct05) <- colnames(frwm)[c(9,1:8,10)]
colnames(median) <- colnames(frwm)[c(9,1:8,10)]

write.table(mean, "fwrTSY.FSCT.mean.txt", sep="\t")
write.table(pct95,"fwrTSY.FSCT.95pct.txt", sep="\t")
write.table(pct05,"fwrTSY.FSCT.5pct.txt" , sep="\t")
write.table(median,"fwrTSY.FSCT.median.txt" , sep="\t")
