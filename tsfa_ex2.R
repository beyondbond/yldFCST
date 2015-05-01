#-------------------------------------
toInstall <- c("tsfa","GPArotation", "setRNG", "tframe", "dse","EvalEst","CDNmoney")
#install.packages(toInstall, repos = "http://cran.r-project.org",dependencies=T)
lapply(toInstall, require, character.only=T)
#-------------------------------------

data("WansbeekMeijer", package="GPArotation")
fa.unrotated <- estFAmodel(NetherlandsTV, 2, n.obs=2150, rotation="none" )
fa.varimax <- estFAmodel(NetherlandsTV, 2, n.obs=2150, rotation="Varimax" )
fa.eiv <- estFAmodel(NetherlandsTV, 2, n.obs=2150, rotation="eiv" )
fa.oblimin <- estFAmodel(NetherlandsTV, 2, n.obs=2150, rotation="oblimin" )
cbind(loadings(fa.unrotated), loadings(fa.varimax), loadings(fa.oblimin), loadings(fa.eiv))

