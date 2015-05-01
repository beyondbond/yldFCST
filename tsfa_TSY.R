#-------------------------------------
toInstall <- c("tsfa","GPArotation", "setRNG", "tframe", "dse","EvalEst","CDNmoney","zoo")
#install.packages(toInstall, repos = "http://cran.r-project.org",dependencies=T)
lapply(toInstall, require, character.only=T)
#-------------------------------------

y <- read.table("frwTSY.txt",header=TRUE,sep="\t")
y$date=as.Date(as.character(y$date),"%Y%m%d")
ym <- data.matrix(y)
ym <- ym[,-1]
P  <- zoo(ym,y$date,frequency=365)
P_ln <- log(P/100)
Y_ln <- log(ym/100)
MBandCredit <- tfwindow(Y_ln, start=1,end=Tobs(Y_ln))
#MBandCredit =data.matrix(Y_ln, rownames.force = NA)
tfplot(MBandCredit, graphs.per.page=3)
tfplot(diff(MBandCredit), graphs.per.page=3)
#lapply(as.data.frame(P_ln),function(x) plot(y$date,x,type='l',xlab='Date',ylab='Value'))
start(MBandCredit)
end(MBandCredit)
Tobs(MBandCredit)
DX <- diff(MBandCredit, lag=1)
Tobs(MBandCredit)
nseries(MBandCredit)
colMeans(DX)
sqrt(diag(cov(DX)))
zz <- eigen(cor(diff(MBandCredit, lag=1)), symmetric=TRUE)[["values"]]
print(zz)
par(omi=c(0.1,0.1,0.1,0.1),mar=c(4.1,4.1,0.6,0.1))
plot(zz, ylab="Value", xlab="Eigenvalue Number", pch=20:20,cex=1,type="o")
z <- FAfitStats(MBandCredit)
print(z, digits=3)
c2withML <- estTSF.ML(MBandCredit, 2)
z <- matrix(0,10,3)
z[matrix(c( 1,6,2,1:3),3,2)] <- c(10, 56, 41)
c3withML <- estTSF.ML(MBandCredit, 3, BpermuteTarget=z)
z <- matrix(0,10,4)
z[matrix(c( 1,6,2,7,1:4),4,2)] <- c(13, 54, 37, 24)
c4withML <- estTSF.ML(MBandCredit, 4, BpermuteTarget=z)
z <- matrix(0,10,5)
z[matrix(c( 1,6,2,7,5,1:5),5,2)] <- c(13, 67, 34, 30, 72)
c5withML <- estTSF.ML(MBandCredit, 5, BpermuteTarget=z)
print(DstandardizedLoadings(c4withML) )
print(c4withML$Phi, digits=3)
print(1 - c4withML$stats$uniquenesses)
print(1 - c2withML$stats$uniquenesses)
print(1 - c3withML$stats$uniquenesses)
print(1 - c5withML$stats$uniquenesses)
print(loadings(c4withML) )
tfplot(ytoypc(factors(c4withML)),
	Title= "Factors from 4 factor model (year-to-year growth rate)",
	lty=c("solid"),
	col=c("black"),
	xlab=c(""),ylab=c("factor 1","factor 2","factor 3","factor 4"),
	par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
	reset.screen=TRUE)
tfplot(factors(c4withML),
	Title="Factors from 4 factor model",
	lty=c("solid"),
	col=c("black"),
	xlab=c(""),ylab=c("factor 1","factor 2","factor 3","factor 4"),
	par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
	reset.screen=TRUE)
z <- explained(c4withML)
tfplot(ytoypc(MBandCredit), ytoypc(z), series=1:5, graphs.per.page=5,
	lty=c("solid", "dashed"),
	col=c("black", "red"),
	ylab=c("currency", "personal cheq.", "NonbankCheq",
	"N-P demand & notice", "N-P term"),
	ylim=list(NULL,NULL,c(-100,100),NULL,NULL),
	Title=
	"Explained money indicator 1-5 (year-to-year growth rate) using 4 factors",
	par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
	reset.screen=TRUE)
tfplot(ytoypc(MBandCredit), ytoypc(explained(c4withML)), series=6:10,
	graphs.per.page=5,
	lty=c("solid", "dashed"),
	col=c("black", "red"),
	ylab=c("","","","","",
	"Investment","Consumer Credit", "Residential Mortgage",
	"Short Term Business Credit", "Other Business Credit"),
	Title=
	"Explained money indicator 6-10 (year-to-year growth rate)using 4 factors",
	par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
	reset.screen=TRUE)
tfplot( MBandCredit, explained(c4withML), series=1:5, graphs.per.page=5,
	lty=c("solid", "dashed"),
	col=c("black", "red"),
	Title= "Explained money indicators 1-5 using 4 factors",
	par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
	reset.screen=TRUE)
tfplot( MBandCredit, explained(c4withML), series=6:10, graphs.per.page=5,
	lty=c("solid", "dashed"),
	col=c("black", "red"),
	Title= "Explained money indicator 6-10 using 4 factors",
	par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
	reset.screen=TRUE)
tfplot( diff(MBandCredit), diff(explained(c4withML)), graphs.per.page=3)
summary(MBandCredit)
DstandardizedLoadings(c2withML)
print(c2withML$Phi, digits=3)
DstandardizedLoadings(c3withML)
print(c3withML$Phi, digits=3)
print(DstandardizedLoadings(c5withML), digits=3)
print(c5withML$Phi, digits=3)
DstandardizedLoadings(c4withML)
print(c2withML$Phi, digits=3)
tfplot(ytoypc(factors(c4withML)), ytoypc(factors(c2withML)),
	ytoypc(factors(c3withML)),
	ytoypc(factors(c5withML)), series=1:2,
	xlab=c(""),ylab=c("factor 1","factor 2"),
	lty=c("solid", "dotdash", "dashed", "dotted"),
	col=c("black","green","red","blue"),
	Title= paste("Factors transaction and long ",
	"(year-to-year growth rate) using 2, 3, 4 and 5 factor models", sep=""),
	par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
	reset.screen=TRUE)
tfplot(ytoypc(factors(c4withML)),
	ytoypc(factors(c3withML)),
	ytoypc(factors(c5withML)), series=3,
	lty=c("solid", "dashed", "dotted"),
	xlab=c(""),ylab=c("","","factor 3"),
	col=c("black","red","blue"),
	Title= paste("Factor near ",
	"(year-to-year growth rate) using 3, 4 and 5 factor models", sep=""),
	par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
	reset.screen=TRUE)
c4withMLg0.5 <- estTSF.ML(MBandCredit, 4, BpermuteTarget=loadings(c4withML),
	rotation="oblimin", rotationArgs=list(gam=0.5))
loadings(c4withMLg0.5)
DstandardizedLoadings(c4withMLg0.5)
DstandardizedLoadings(c4withMLg0.5) - DstandardizedLoadings(c4withML)
summary(c4withMLg0.5)
c4withMLgneg0.5 <- estTSF.ML(MBandCredit, 4, BpermuteTarget=loadings(c4withML),
	rotation="oblimin", rotationArgs=list(gam=-0.5))
loadings(c4withMLgneg0.5)
DstandardizedLoadings(c4withMLgneg0.5)
DstandardizedLoadings(c4withMLgneg0.5) - DstandardizedLoadings(c4withML)
summary(c4withMLgneg0.5)
c4withMLgneg1.0 <- estTSF.ML(MBandCredit, 4, BpermuteTarget=loadings(c4withML),
	rotation="oblimin", rotationArgs=list(gam=-1.0))
loadings(c4withMLgneg1.0)
DstandardizedLoadings(c4withMLgneg1.0)
DstandardizedLoadings(c4withMLgneg1.0) - DstandardizedLoadings(c4withML)
summary(c4withMLgneg1.0)
c4withMLbQ <- estTSF.ML(MBandCredit, 4, rotation="bentlerQ",
	BpermuteTarget=loadings(c4withML))
loadings(c4withMLbQ)
DstandardizedLoadings(c4withMLbQ)
DstandardizedLoadings(c4withMLbQ) - DstandardizedLoadings(c4withML)
summary(c4withMLbQ)
tfplot(ytoypc(factors(c4withML)), ytoypc(factors(c4withMLg0.5)),
	ytoypc(factors(c4withMLgneg0.5)), ytoypc(factors(c4withMLgneg1.0)),
	ytoypc(factors(c4withMLbQ)),
	xlab=c(""),ylab=c("factor 1","factor 2","factor 3","factor 4"),
	lty=c("solid", "dashed", "dotted", "dotdash", "longdash"),
	col=c("black","red","blue","green","pink"),
	Title= paste(
	"Factors from various 4 factor models (year-to-year growth rate)",
	"\n and oblimin with gam=0 (solid)"),
	par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
	reset.screen=TRUE)
c4withMLgm <- estTSF.ML(MBandCredit, 4, rotation="geominQ",
	BpermuteTarget=loadings(c4withML))
loadings(c4withMLgm)
DstandardizedLoadings(c4withMLgm)
DstandardizedLoadings(c4withMLgm) - DstandardizedLoadings(c4withML)
	summary(c4withMLgm)
	tfplot(ytoypc(factors(c4withML)), ytoypc(factors(c4withMLgm)),
	xlab=c(""),ylab=c("factor 1","factor 2","factor 3","factor 4"),
	lty=c("solid", "dashed"),
	col=c("black","red"),
	Title= paste(
	"Factors from geomin (dashed) 4 factor model (year-to-year growth rate)",
	"\n and oblimin with gam=0 (solid)"),
	par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
	reset.screen=TRUE)
c4withMLnotNorm <- estTSF.ML(MBandCredit, 4, normalize=FALSE,
	BpermuteTarget=loadings(c4withML))
DstandardizedLoadings(c4withMLnotNorm)
DstandardizedLoadings(c4withML) - DstandardizedLoadings(c4withMLnotNorm)
z <- matrix(0,10,4)
z[matrix(c( 1,6,2,7,1:4),4,2)] <- c(11, 104, 20, 13)
c4withMLbefore90 <- estTSF.ML(tfwindow(MBandCredit, end=c(1989,12)), 4,
	BpermuteTarget=z)
c4withMLafter95 <- estTSF.ML(tfwindow(MBandCredit, start=c(1995,1)), 4,
	BpermuteTarget=loadings(c4withML))
z <- matrix(0,10,4)
z[matrix(c( 1,6,2,7,1:4),4,2)] <- c(11, 104, 20, 13)
c4withMLbefore95 <- estTSF.ML(tfwindow(MBandCredit, end=1000), 4,
	BpermuteTarget=z)
c4withMLafter00 <- estTSF.ML(tfwindow(MBandCredit, start=1000), 4,
	BpermuteTarget=loadings(c4withML))
c4withML90to00 <- estTSF.ML(tfwindow(MBandCredit, start=2000), 4,
	BpermuteTarget=loadings(c4withML))
tfplot(ytoypc(factors(c4withML)), ytoypc(factors(c4withMLbefore90)),
ytoypc(factors(c4withMLbefore95)), ytoypc(factors(c4withMLafter95)),
	ytoypc(factors(c4withMLafter00)), ytoypc(factors(c4withML90to00)),
	xlab=c(""),ylab=c("factor 1","factor 2","factor 3","factor 4"),
	ylim=list(NULL,c(-20,20),c(-25,40),NULL),
	graphs.per.page=4,
	lty=c("dashed", "dotted", "dotdash", "longdash", "dotted",
	"twodash"),
	col=c("red","blue","green","pink","violet","brown"),
	Title= paste(
	"Factors (year to year growth) using full sample and sub-samples\n",
	"ML estimation with quartimin rotation objective", sep=""),
	par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
	reset.screen=TRUE)
tfplot(ytoypc(MBandCredit), ytoypc(explained(c4withML)),
	ytoypc(explained(c4withMLbefore90)), ytoypc(explained(c4withMLbefore95)),
	ytoypc(explained(c4withMLafter95)), ytoypc(explained(c4withMLafter00)),
	ytoypc(explained(c4withML90to00)), series=1:5, graphs.per.page=5,
	ylab=c("currency", "personal cheq.", "NonbankCheq",
	"N-P demand & notice", "N-P term"),
	ylim=list(NULL,NULL,c(-70,70),NULL,c(-70,70)),
	lty=c("solid", "dashed", "dotted", "dotdash", "longdash", "dotted",
	"twodash"),
	col=c("black", "red","blue","green","pink","violet","brown"),
	Title= paste("Explained money indicators 1-5 (year to year growth)\n",
	"using 4 factors, full sample and sub-samples", sep=""),
	par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
	reset.screen=TRUE)
tfplot(ytoypc(MBandCredit), ytoypc(explained(c4withML)),
	ytoypc(explained(c4withMLbefore90)), ytoypc(explained(c4withMLbefore95)),
	ytoypc(explained(c4withMLafter95)), ytoypc(explained(c4withMLafter00)),
	ytoypc(explained(c4withML90to00)), series=6:10, graphs.per.page=5,
	ylab=c("","","","","",
	"Investment","Consumer Credit", "Residential Mortgage",
	"Short Term Business Credit", "Other Business Credit"),
	lty=c("solid", "dashed", "dotted", "dotdash", "longdash", "dotted",
	"twodash"),
	col=c("black", "red","blue","green","pink","violet","brown"),
	Title= paste("Explained money indicators 6-10 (year to year growth)\n",
	"using 4 factors, full sample and sub-samples", sep=""),
	par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
	reset.screen=TRUE)
