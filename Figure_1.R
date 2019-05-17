load("Figure_1.RData")
FROHthresh1 <- max(xsens1[which(ysens1 >= yspec1)])
ibest1      <- which.min(abs(xsens1-FROHthresh1))

FROHthresh2 <- max(xsens2[which(ysens2 >= yspec2)])
ibest2      <- which.min(abs(xsens2-FROHthresh2))

png("Figure_1.png",width=2000,height=1800,res=300)
Cols <- c("lightblue","lightgreen")
op <- par(mfrow=c(2,2))
par(mar=c(5,5,2,2))
plot(c(0,0),c(1,1),xlim=c(0,1),ylim=c(0,1),
     xlab="Specificity",ylab="Sensitivity",axes=FALSE,type="n",
     main="MT1 vs. MT2")
axis(1,at=seq(0,1,by=.2),labels=c("1.0","0.8","0.6","0.4","0.2","0.0"))
axis(2)
lines(1-yspec1,ysens1,col="coral1",lwd=2)
segments(x0=0,y0=0,x1=1,y1=1,col="goldenrod",lty=2)
arrows(x0=0.3,y0=0.9,x1=0.1,y1=0.9,col="blue",length = 0.05,lty=1,lwd=1)
text(0.51,0.9,expression(paste(F[ROH]," = 0.174")),col="blue",cex=1.0)
segments(x0=1-yspec1[ibest1],y0=0,x1=1-yspec1[ibest1],y1=ysens1[ibest1],col="grey",lty=3)
segments(x0=0,y0=ysens1[ibest1],x1=1-yspec1[ibest1],y1=ysens1[ibest1],col="grey",lty=3)
text(0.32,0.6,paste0("AUC = ",round(auc1,3)),font=2)
mtext("a",at=0,side=3,line=0,font=2,cex=1.2)

par(mar=c(5,4,2,2))
plot(c(0,0),c(1,1),xlim=c(0,1),ylim=c(0,1),xlab="Specificity",
     ylab="Sensitivity",axes=FALSE,type="n",main="MT1/MT2 vs. MT3")
axis(1,at=seq(0,1,by=.2),labels=c("1.0","0.8","0.6","0.4","0.2","0.0"))
axis(2)
lines(1-yspec2,ysens2,col="coral1",lwd=2)
segments(x0=0,y0=0,x1=1,y1=1,col="goldenrod",lty=2)
arrows(x0=0.3,y0=0.86,x1=0.13,y1=0.86,col="blue",length = 0.05,lty=1,lwd=1)
text(0.5,0.86,expression(paste(F[ROH]," = 0.087")),col="blue",cex=1.0) ## That value is obtained as FROHthresh2

segments(x0=1-yspec2[ibest2],y0=0,x1=1-yspec2[ibest2],y1=ysens2[ibest2],col="grey",lty=3)
segments(x0=0,y0=ysens2[ibest2],x1=1-yspec2[ibest2],y1=ysens2[ibest2],col="grey",lty=3)
text(0.34,0.6,paste0("AUC = ",round(auc2,3)),font=2)
mtext("b",at=0,side=3,line=0, font=2,cex=1.2)

par(mar=c(5,5,0,2))
matplot(cbind(xsens1,xspec1),cbind(ysens1,yspec1),type="l",col=Cols,ylim=c(0,1),lwd=2,
        xlab=expression(paste(F[ROH]," threshold")),ylab="Sensitivity / Specificity",axes=FALSE)
axis(1,at=c(0.1,0.2,0.3,0.4));axis(2,at=c(0,.2,0.4,0.6,0.8,1.0))
#legend(0.25,0.35,legend=c("Sensitivity","Specificity"),box.lty=0,cex=0.9,lty=1,col=Cols,lwd=3)
abline(v=FROHthresh1,col="red",lwd=2,lty=1)
arrows(xsens1[ibest1]+0.1,ysens1[ibest1],xsens1[ibest1]+0.01,ysens1[ibest1],col="blue",length = 0.1,lty=1,lwd=1.2)
text(xsens1[ibest1]+0.175,ysens1[ibest1],expression(paste(F[ROH]," = 0.174")),col="blue")
mtext("c",at=0.1,side=3,line=0, font=2,cex=1.2)

par(mar=c(5,4,0,2))
matplot(cbind(xsens2,xspec2),cbind(ysens2,yspec2),type="l",col=Cols,ylim=c(0,1),lwd=2,
        xlab=expression(paste(F[ROH]," threshold")),ylab="Sensitivity / Specificity",axes=FALSE)
axis(1,at=c(0,0.1,0.2,0.3,0.4));axis(2,at=c(0,.2,0.4,0.6,0.8,1.0))
arrows(xsens2[ibest2]+0.1,ysens2[ibest2],xsens2[ibest2]+0.01,ysens2[ibest2],col="blue",length = 0.1,lty=1,lwd=1.2)
text(xsens2[ibest2]+0.19,ysens2[ibest2],expression(paste(F[ROH]," = 0.087")),col="blue") ## That value is obtained as FROHthresh2

abline(v=c(FROHthresh2),col=c("red"),lwd=2,lty=c(1))
legend(0.25,0.35,legend=c("Specificity","Sensitivity"),box.lty=0,cex=0.9,lty=1,col=Cols,lwd=3)
mtext("d",at=0,side=3,line=0, font=2,cex=1.2)
par(op)
dev.off()
