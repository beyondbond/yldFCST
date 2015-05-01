doInstall <- FALSE  # Change to FALSE if you don't want packages installed.
toInstall <- c("tsfa","GPArotation", "setRNG", "tframe", "dse","EvalEst","CDNmoney")
if(doInstall){install.packages(toInstall, repos = "http://cran.r-project.org")}
lapply(toInstall, library, character.only = TRUE)


#sv1 <- read.table(pipe("cut -f2- frwTSY.txt"),header=TRUE,sep="\t")
sv1 <- read.table("frwTSY.txt",header=TRUE,sep="\t")
rdate=as.Date(as.character(sv1$date),"%Y%m%d")
sv1$date=rdate
sv1$date=NULL
y <- as.matrix(sv1)
y <- log(y/100)
Dv1 <- diff(y,lag=1)
summary(Dv1)
cov(Dv1)
q()
cor(Dv1)
eigen(cor(Dv1))

Dv1.pca <- princomp(Dv1)
summary(Dv1.pca )
plot(Dv1.pca)

Dv1.fa1 <- factanal(Dv1,factors=4,rotation="none",scores="Bartlett",trace=T)
Dv1.fa1
attributes(Dv1.fa1)

Dv1.fa2 <- factanal(Dv1,factors=4,rotation="varimax",scores="Bartlett",trace=T)
coef <- solve(Dv1.fa2$correlation) %*% Dv1.fa2$loadings
means=cbind(lapply(as.data.frame(Dv1),mean))
sds=cbind(diag(cov(Dv1))^.5)
#scale(Dv1, means, sds) %*% coef
predicts=scale(Dv1, center=T,scale=T) %*% coef
head(predicts)


#Dv1.fa2
#attributes(Dv1.fa2)

#print(Dv1.fa2$loadings, cutoff = .30, sort = TRUE);
