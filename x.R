toInstall <- c("tsfa","GPArotation", "setRNG", "tframe", "dse","EvalEst","CDNmoney")
#install.packages(toInstall, repos = "http://cran.r-project.org",dependencies=T)
lapply(toInstall, require, character.only=T)

z <- ts(rnorm(100), start=c(1982,1), frequency=12)
Dz <- tframed(diff(z), tfTruncate(tframe(z), start=2))
tframe(Dz)
IDz <- tframed(cumsum(c(0, Dz)), tfExpand(tframe(Dz), add.start=1))
tframe(IDz)
tframe(tfTruncate(z, start=5))
data("CanadianMoneyData.asof.28Jan2005", package="CDNmoney")
data("CanadianCreditData.asof.28Jan2005", package="CDNmoney")
cpi <- 100 * M1total / M1real
seriesNames(cpi) <- "CPI"
popm <- M1total / M1PerCapita
seriesNames(popm) <- "Population of Canada"
z <- tframed(tbind(
 MB2001,
 MB486 + MB452 + MB453 ,
 NonbankCheq,
 MB472 + MB473 + MB487p,
 MB475,
 NonbankNonCheq + MB454 + NonbankTerm + MB2046 + MB2047 + MB2048 +
 MB2057 + MB2058 + MB482),
names=c("currency", "personal cheq.", "NonbankCheq",
 "N-P demand & notice", "N-P term", "Investment" )
 )
TotalMoney <- tframed(rowSums(z), tframe(z))
z <- tbind (z, ConsumerCredit, ResidentialMortgage,
 ShortTermBusinessCredit, OtherBusinessCredit)
z <-tfwindow(z, start=c(1981,11), end=c(2004,11))
 scale <- tfwindow(1e8 /(popm * cpi), tf=tframe(z))
 MBandCredit <- sweep(z, 1, scale, "*")
tfplot(MBandCredit, graphs.per.page=3)
