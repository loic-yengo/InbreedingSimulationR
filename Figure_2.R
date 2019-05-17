load("Figure_2.RData")
mininmumROHlengthMB <- 1.5
png("Figure_2.png",width=1000,height=900,res=200)
ax <- hist(segLengthMB-mininmumROHlengthMB,nclass=30,probability = TRUE,border="grey",col="lightblue",
           xlab=paste0("ROH length (in Mb) - ",mininmumROHlengthMB," Mb"),
           ylab="Density",main="",cex.lab=0.8)
x  <- sort(segLengthMB-mininmumROHlengthMB)
## md <- getParMixExp(x) has been pre-calculated

y1 <- dexp(x,rate=1/md["MeanMinor.rate2"])
y2 <- dexp(x,rate=1/md["MeanMajor.rate1"])
y  <- md["ProbMinor.prob2"]*y1+md["ProbMajor.prob1"]*y2

lines(x,y,col="coral1",lwd=2,lty=3)
legend(10,0.08,title="Fitted Mixture of Exponential Distributions",
       legend=expression(paste(p(x)," = 0.84 [0.06 exp(-0.06x)] + 0.16 [1.4 exp(-1.4x)]")), ## parameters from Table 1
       lty=3,col="coral1",lwd=2,box.lty=0,cex=0.7)
dev.off()
