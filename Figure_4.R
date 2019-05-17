load("Figure_4.RData")
Cols <- c("gray45","dodgerblue","coral1","black")
Pchs <- c(19,1,5,2)

png("Figure_4.png",width=2100,height=1200,res=300)
par(mar=c(5,6,3,2))
plot(R[,"mobs"],pch=Pchs[1],ylim=c(-1.53,1.53),axes=FALSE,col=Cols[1],
     ylab="Phenotypic reduction (unit: trait SD)\nin EI cases vs. EI controls",
     xlab="Inbreeding associated traits",cex.lab=0.9)
dx <- 0.1
for(i in 1:10){
  points(i-dx,R[i,"mpred.fit"],col=Cols[2],pch=Pchs[2])
  segments(i-dx,R[i,"pred_mse_low.fit"],i-dx,R[i,"pred_mse_up.fit"],col=Cols[2],lty=2)
  points(i+dx,Rfroh[i,"mpred.fit"],col=Cols[3],pch=Pchs[3])
  segments(i+dx,Rfroh[i,"pred_mse_low.fit"],i+dx,Rfroh[i,"pred_mse_up.fit"],col=Cols[3],lty=3)
  
  segments(i,R[i,"mobs_low"],i,R[i,"mobs_up"],col=Cols[1],lty=1)
}
abline(h=0,col="grey",lty=4)
axis(1,at=1:10,cex.axis=0.75,tick=FALSE,
     labels=Code[rownames(R)],
     col=Cols[1])
axis(2,cex.axis=0.9,at=seq(-1.5,1.5,by=0.5),las=2)
legend(1,1.6,legend=c("Observed in EI cases",
                      expression(paste("Predicted (",F[UNI],")  [99.5% CI]")),
                      expression(paste("Predicted (",F[ROH],") [99.5% CI]"))),
       lty=c(1,2,3),
       col=Cols,box.lty=0,pch=Pchs,cex=0.8)
dev.off()

