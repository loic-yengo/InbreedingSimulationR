load("Supplementary_Figure_2.RData")
cm1  <- 58.483247
cm2  <- 61.994151
Cols <- c("dodgerblue","coral1","goldenrod","khaki3","lightblue")
png("Supplementary_Figure_2.png",width=1800,height=1100,res=300)
par(mar=c(5,2,2,0))
plot(c(0,200),c(0,1),ylim=c(0,1.05),type="n",axes=FALSE,xlab="Genomic position (Mb) on the X-chromosome",ylab="",
     main="X-chromosome ROH",cex.main=1)
for(i in 1:nrow(ROHx)){
  segments(ROHx[i,"POS1"],ROHx[i,"Height"],ROHx[i,"POS2"],ROHx[i,"Height"],col=Cols[ROHx[i,"Grp"]],lwd=3)
}
rect(xleft=cm1,ybottom=0,xright=cm2,ytop=1,col="grey",border=0)
mains <- c("[1.5 - 5)",
           "[5 - 10)",
           "[10 - 20)",
           "[20 - 50)",
           "[50 - 100)")
text(x=0.5*(cm1+cm2),1.06,"Centromere",cex = 0.5)
arrows(x0=0.5*(cm1+cm2),x1 = 0.5*(cm1+cm2),y0=1.03,y1=1,lwd=1,col=2,len=0.025)
legend(160,1,legend=mains,horiz==TRUE,box.lty=0,border=0,fill=Cols,cex=0.7,title="ROH Length (Mb)",title.col = 2)
axis(1,at=seq(0,160,by=40))
mtext(text=paste0("Individuals (N=",length(unique(ROHx$IID)),")"),line = 0,side = 2,cex=1.0)
dev.off()
