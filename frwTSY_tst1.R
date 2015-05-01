toInstall <- c("fGarch","zoo")
lapply(toInstall, require, character.only=T)
#################################### Utilities #################################### 
rateMC <- function(p, z, Nsim, Nforecast){
  data <- frwdiff[,p] # Use jth rate as benchmark
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
  volbound=mean(fit@sigma.t)*4
  rho <- cor(frwm[,c(j,9)])[1,2] #correlation
  rho =0
  Fmatrix <- matrix(rep(0, Nsim*Nforecast), nrow=Nsim, ncol=Nforecast)
  Smatrix <- matrix(rep(0, Nsim*Nforecast), nrow=Nsim, ncol=Nforecast)
  #w <- matrix(rnorm(Nsim*Nforecast), nrow=Nsim, ncol=Nforecast)
  w1 <- matrix(rnorm(Nsim/2*Nforecast), nrow=Nsim/2, ncol=Nforecast)
  w2 <- matrix(rnorm(Nsim/2*Nforecast), nrow=Nsim/2, ncol=Nforecast)
  w = rbind(w1,w2)
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
      garch.std = pmin(garch.std,volbound)
      garch.unstd <- garch.std*rho*zi+garch.std*sqrt(1-rho^2)*wi
      garch.unstd = pmin(garch.unstd,volbound)
      Fmatrix[,i] <- coef[1] + coef[2]*Fmatrix[,i-1] +coef[3]*Fmatrix[,i-2]+ garch.unstd

#      x1=pmin(Fmatrix[,i],Fmatrix[,i-1]+volbound/2)
#      Fmatrix[,i] = pmax(Fmatrix[,i-1]-volbound/2,x1)
    }
    # Smatrix[,i] <- garch.std
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

# calc boundry
rightShift <- function(vx,boundry) {
   nx=length(vx)-1
   for(j in nx:1) {
    if(j>3 && vx[j]>1.0*boundry){
	vx[j]=vx[j+1] 
    }else if(j<=3 && vx[j]>boundry){
	vx[j]=vx[j+1] 
    }

#    if(vx[j]>vx[j+1]*2){
#	vx[j]=vx[j+1] 
#    }
   }
   return(vx)
 }

fqStats <- function(fp) {
  path <- "/apps/fafa/rx/tst/"
  #frw <- read.table(paste(path,filename,sep=""), header=TRUE, sep="\t")
  frw <- read.table(fp, header=TRUE, sep="\t")
  frw[frw<0.1] <- 0.1 # replace rates < 0.1 with 0.1
  frwm <- frw[,2:11]
  frwdiff <- diff(ts(log(frwm/100), start=1))
  summary(frwdiff)
  b <- sqrt(diag(cov(frwdiff)))
  return b
}

files <- c("frwSWP.Y.txt", "frwSWP.M.txt", "frwSWP.txt")
lapply(files, fqStats)
q()

##################################### Load data ####################################
# Load forecast date vector
fdate=as.Date(as.character(read.table("fcstDateTSY.dat",header=TRUE,sep="\t")$date),"%Y%m%d")
#fdate=as.Date(as.character(read.table("fcstDateSWP.dat",header=TRUE,sep="\t")$date),"%Y%m%d")

#path <- "Users/Jie/Google Drive/Work/"
path <- "/apps/fafa/rx/tst/"
filename <- "frwTSY.txt"
#filename <- "frwSWP.txt"
#frw <- read.table(paste(path,filename,sep=""), header=TRUE, sep="\t")
frw <- read.table(filename, header=TRUE, sep="\t")
frw[frw<0.1] <- 0.1 # replace rates < 0.1 with 0.1
frwm <- frw[,2:11]
frwdiff <- diff(ts(log(frwm/100), start=1))

############################## Monte Carlo simulation ##############################
################################# Benchmark rates ##################################
Nsim <- 1000 # number of simulations
Nforecast <- 2500 # how many days ahead to forecast
p <- 9 # use pth rate as benchmark
# generate Gaussian random numbers
z <- matrix(rnorm(Nsim*Nforecast), nrow=Nsim, ncol=Nforecast)
# forecast matrix of log diff
#Fmatrix <- rateMC(p, z, Nsim, Nforecast)
# convert log diff to normal rates
#frwFRCT <- diff2frw(p, Fmatrix)

########################## Other rates' forecast log diff ###########################
#r <- c(1:10)
#c <- r[!r==p] # rates excluding the benchmark rate
#for (j in c){
#  Fmatrix <- corMC(j)
#  filepath <- "/Users/Jie/Google Drive/Work/Output/"
#  write.table(Fmatrix, paste(filepath,"Fmatrix",colnames[j],".txt",sep=""),sep="\t")
#}

################################## Save quantiles ##################################
mean <- c(1:Nforecast)
pct95 <- c(1:Nforecast)
pct5 <- c(1:Nforecast)
median <- c(1:Nforecast)
# Benchmark rate quantiles
#mean <- apply(frwFRCT, 2, mean)
#pct95 <- apply(frwFRCT,2,quantile,0.95)
#pct5 <- apply(frwFRCT,2,quantile,0.05)
#median <- apply(frwFRCT,2,quantile,0.5)
for (j in 1:10){
  #Fmatrix <- read.table(paste(filepath,"Fmatrix",colnames[j],".txt",sep=""),sep="\t")
  if (j==p){
    Fmatrix <- rateMC(p, z, Nsim, Nforecast)
  } else {
    Fmatrix <- corMC(j)
  }
  frwFRCTj <- diff2frw(j, Fmatrix)
  mean <- cbind(mean, apply(frwFRCTj,2,mean))
  pct95 <- cbind(pct95, apply(frwFRCTj,2,quantile,0.95))
  pct5 <- cbind(pct5, apply(frwFRCTj,2,quantile,0.05))
  median <- cbind(median, apply(frwFRCTj,2,quantile,0.5))
}

for (i in 1:Nforecast){
  mean[i,2:11] = rightShift(mean[i,2:11],11) 
  pct95[i,2:11] = rightShift(pct95[i,2:11],11) 
  pct5[i,2:11] = rightShift(pct5[i,2:11],11) 
  median[i,2:11] = rightShift(median[i,2:11],11) 
}

colnames <- append(colnames(frwm),"Nsim",0)
colnames(mean) <- colnames
colnames(pct95) <- colnames
colnames(pct5) <- colnames
colnames(median) <- colnames

#f <- substring(filename,1,6)
#write.table(mean, paste(f,".FCST.mean",sep=""), sep="\t", row.names = F)
#write.table(pct95,paste(f,".FCST.95pct",sep=""), sep="\t", row.names = F)
#write.table(pct5,paste(f,".FCST.5pct",sep=""), sep="\t", row.names = F)
#write.table(median,paste(f,".FCST.median",sep=""), sep="\t", row.names = F)


#for (i in 2:11){
#  plot(mean[,i],type="l",main=paste(colnames[i],"mean"))
#  plot(pct5[,i],type="l",main=paste(colnames[i],"pct5"))
#  plot(median[,i],type="l",main=paste(colnames[i],"median"))
#}

Ndays <- Nforecast # how many days for plotting
#for (i in 2:11){
for (i in 10:10){
	lapply(as.data.frame(pct95[1:Ndays,i]),function(x) plot(fdate[1:Ndays],x,type='l',xlab='Date',ylab='Value',main=paste(colnames[i],"pct95")))
	lapply(as.data.frame(pct5[1:Ndays,i]),function(x) plot(fdate[1:Ndays],x,type='l',xlab='Date',ylab='Value',main=paste(colnames[i],"pct5")))
	lapply(as.data.frame(mean[1:Ndays,i]),function(x) plot(fdate[1:Ndays],x,type='l',xlab='Date',ylab='Value',main=paste(colnames[i],"mean")))
}

#write.table(mean, paste(filepath,"fwrSWP.FCST.mean.txt",sep="") , sep="\t", row.names = F)
#write.table(pct95, paste(filepath,"fwrSWP.FCST.95pct.txt",sep="") , sep="\t", row.names = F)
#write.table(pct5, paste(filepath,"fwrSWP.FCST.5pct.txt",sep="") , sep="\t", row.names = F)
#write.table(median, paste(filepath,"fwrSWP.FCST.median.txt",sep="") , sep="\t", row.names = F)
