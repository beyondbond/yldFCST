toInstall <- c("fGarch","zoo")
lapply(toInstall, require, character.only=T)

############################# Load data ############################
y <- read.table("h15Hist.tsv", header=TRUE, sep="\t")
b1m <- y[,2]
summary(b1m)
q()
rdate=as.Date(as.character(y$date),"%Y%m%d")
y$date=NULL
ym <- data.matrix(y)

#- for Time-Series Plotting
yt  <- zoo(ym,rdate,frequency=365)

ylog <- log(ym/100)
yDiff <- diff(ym)
yLDiff <- diff(ylog)
n <- nrow(yLDiff)
cor(yLDiff)
summary(yLDiff)
cor(yDiff)
summary(yDiff)
q()

###################### Run AR(2)+GARCH(1,1) ########################
re1 <- garchFit(formula=~arma(2,0)+garch(1,1),
                   data=yLDiff[,1], cond.dist=c("norm"))
coef1 <- coef(re1)
mean1 <- coef1[1]+coef1[2]*yLDiff[n,1]+coef1[3]*yLDiff[n-1,1] #AR(2)
std1 <- tail(re1@sigma.t, n=1)
unstd1 <- tail(re1@residuals, n=1)
garch.std1 <- coef1[4]+coef1[5]*unstd1+coef1[6]*std1

re2 <- garchFit(formula=~arma(2,0)+garch(1,1),
                data=yLDiff[,2], cond.dist=c("norm"))
coef2 <- coef(re2)
mean2 <- coef2[1]+coef2[2]*yLDiff[n,2]+coef2[3]*yLDiff[n-1,2] #AR(2)
std2 <- tail(re2@sigma.t, n=1)
unstd2 <- tail(re2@residuals, n=1)
garch.std2 <- coef2[4]+coef2[5]*unstd2+coef2[6]*std2

rho <- cor(yLDiff[,1:2])[1,2] # correlation

##################### Monte Carlo simulation #######################
# Generate Gaussian random number
Nsim <- 1000 # number of random numbers
Nforecast <- 1 # how many days ahead to forecast
z <- matrix(rnorm(Nsim*Nforecast), nrow=Nsim, ncol=Nforecast)
w <- matrix(rnorm(Nsim*Nforecast), nrow=Nsim, ncol=Nforecast)
# Forecast matrix
fB3M <- matrix(rep(mean1, Nsim*(Nforecast+1)), nrow=Nsim, ncol=(Nforecast+1))
fB6M <- matrix(rep(mean2, Nsim*(Nforecast+1)), nrow=Nsim, ncol=(Nforecast+1))
for (i in 2:(Nforecast+1)){
  fB3M[,i] <- fB3M[,i-1] + z*garch.std1
  fB6M[,i] <- fB6M[,i-1] + rho*z*garch.std1 + sqrt(1-rho^2)*w*garch.std2
}
