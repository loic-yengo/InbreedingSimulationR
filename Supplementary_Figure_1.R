load("Supplementary_Figure_1.RData")
png("Supplementary_Figure_1.png",width=1700,height=1500,res=300)
dxp <- 0.1
par(mar=c(4,5,3,2))
ax <- boxplot(LRT[,"Froh_1.5_XC"]~LRT[,"Z"],axes=FALSE,border=0,ylim=c(0,1),
              ylab=expression(paste(F[ROH]," from X-chromosome")))
points(LRT[,"Z"],LRT[,"Froh_1.5_XC"],pch=19,col="grey",cex=0.8)

for(k in 1:2){
  segments(k,m[k,"Froh_1.5_XC"]-1.96*s[k,"Froh_1.5_XC"],k,m[k,"Froh_1.5_XC"]+1.96*s[k,"Froh_1.5_XC"],col=2)
  segments(k-dxp,m[k,"Froh_1.5_XC"],k+dxp,m[k,"Froh_1.5_XC"],col=2,lwd=5)
}
axis(2)
axis(1,at=1:2,labels=c(paste0("Group 1\nMore likely offspring\nof brother-sister mating\n(N = ",sum(LRT$Z==1),")"),
                       paste0("Group 2\nMore likely offspring\nof parent-offspring mating\n(N = ",sum(LRT$Z==2),")")),tick = FALSE,line=+1.5)
abline(h=c(0.25,0.5),lty=2:3,col="dodgerblue")

legend(0.5,1.0,title="Theoretical expectations",title.col = 4,legend=c("Under brother-sister mating","Under parent-offspring mating"),box.lty=0,col="dodgerblue",lty=2:3,cex=0.8)
legend(1.9,1.0,legend=expression(paste("Mean" %+-%" 95% CI")),box.lty=0,col="red",lty=1,lwd=5,cex=0.7)
dev.off()
