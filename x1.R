toInstall <- c("MARSS","nlme","mvtnorm","KFAS","stats","utils","graphics","Hmisc","maps","xtable","stringr")
#install.packages(toInstall, repos = "http://cran.r-project.org",dependencies=T)
lapply(toInstall, require, character.only=T)
