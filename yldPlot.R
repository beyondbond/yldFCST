#UNFINISHED
toInstall <- c("fGarch","zoo")
lapply(toInstall, require, character.only=T)
##################################### Load data ####################################
zxy frw2spot.egx 'sBdate=20150319;outType="yld";inFile="frwTSY.FCST.5pct";  outFile="yldTSY.FCST.5pct";'
zxy frw2spot.egx 'sBdate=20150319;outType="yld";inFile="frwTSY.FCST.95pct"; outFile="yldTSY.FCST.95pct";'
zxy frw2spot.egx 'sBdate=20150319;outType="yld";inFile="frwTSY.FCST.mean";  outFile="yldTSY.FCST.mean";'
zxy frw2spot.egx 'sBdate=20150319;outType="yld";inFile="frwTSY.FCST.median";outFile="yldTSY.FCST.median";'

zxy frw2spot.egx 'sBdate=20150403;outType="yld";inFile="frwSWP.FCST.5pct";  outFile="yldSWP.FCST.5pct";'
zxy frw2spot.egx 'sBdate=20150403;outType="yld";inFile="frwSWP.FCST.95pct"; outFile="yldSWP.FCST.95pct";'
zxy frw2spot.egx 'sBdate=20150403;outType="yld";inFile="frwSWP.FCST.mean";  outFile="yldSWP.FCST.mean";'
zxy frw2spot.egx 'sBdate=20150403;outType="yld";inFile="frwSWP.FCST.median";outFile="yldSWP.FCST.median";'

# Load forecast date vector
fdate=as.Date(as.character(read.table("fcstDateTSY.dat",header=TRUE,sep="\t")$date),"%Y%m%d")
#fdate=as.Date(as.character(read.table("fcstDateSWP.dat",header=TRUE,sep="\t")$date),"%Y%m%d")

#path <- "Users/Jie/Google Drive/Work/"
path <- "/apps/fafa/rx/tst/"
filename <- "yldTSY.FCST.5pct"
#filename <- "frwSWP.txt"
#frw <- read.table(paste(path,filename,sep=""), header=TRUE, sep="\t")
frw <- read.table(filename, header=TRUE, sep="\t")
frw[frw<0.1] <- 0.1 # replace rates < 0.1 with 0.1
frwm <- frw[,2:11]
frwdiff <- diff(ts(log(frwm/100), start=1))

#for (i in 2:11){
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
