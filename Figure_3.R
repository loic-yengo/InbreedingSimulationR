load("Figure_3.RData")
myGray      <- "gray50"
Cols        <- c("dodgerblue","coral1","goldenrod","khaki3","lightblue","red")

mains <- c("ROH<5Mb\n1258 segments (32.6%)",
           "5Mb<ROH<10Mb\n761 segments (19.7%)",
           "10Mb<ROH<20Mb\n846 segments (21.9%)",
           "20Mb<ROH<50Mb\n847 segments (21.9%)",
           "50Mb<ROH<100Mb\n145 segments (3.8%)",
           "ROH>100Mb\nsegments (0.2%)")
png("Figure_3.png",width=1800,height=1600,res=300)
op <- par(mfrow=c(2,3))
for(k in 1:6){
  dt    <- segs[which(segs[,"Grp"]==k),]
  dt$id <- paste0(letters[dt$CHR],"-",dt$IID)
  dt[,"Height"] <- as.numeric(as.factor(dt$id))
  dt[,"Height"] <- dt[,"Height"] / max(dt[,"Height"])
  dt            <- dt[order(as.numeric(dt$CHR)),]
  dy            <- 0.5/nrow(dt)
  par(mar=c(5,2,2,2))
  plot(c(0,250),c(0,1),ylim=c(0,1.02),type="n",axes=FALSE,xlab="Genomic position (Mb)",ylab="",main=mains[k],cex.main=1)
  for(i in 1:nrow(dt)){
    if(k<6){
      segments(dt[i,"POS1"],dt[i,"Height"],dt[i,"POS2"],dt[i,"Height"],col=Cols[dt[i,"Grp"]],lwd=2)
    }else{
      segments(dt[i,"POS1"],dt[i,"Height"]-0.5*dy,dt[i,"POS2"],dt[i,"Height"]-0.5*dy,col=Cols[dt[i,"Grp"]],lwd=2)
    }
  }
  axis(1)
  mtext(text=paste0("Individuals (N=",length(unique(dt$IID)),")"),line = 0,side = 2,cex=0.8)
  p <- par('usr')
  text(290, mean(p[3:4]), labels = 'Chromosomes', xpd = NA, srt = -90,col=myGray)
  mtext(letters[k],at=-15,side=3,line=-0.5, font=2,cex=1)
  
  ## Add chromosome information
  chroms <- unique(dt$CHR)
  for(chrom in chroms){
    tmp  <- dt[which(dt$CHR==chrom),]
    ymin <- min(tmp[,"Height"])-dy
    ymax <- max(tmp[,"Height"])+dy
    mdy  <- 0.5*(ymin+ymax)
    segments(0,ymin,250,ymin,col=myGray,lty=1,lwd=0.3)
    segments(0,ymax,250,ymax,col=myGray,lty=1,lwd=0.3)
    mtext(text=chrom,at=mdy,line = 0,side = 4,cex=0.5,las=2,col=myGray)
  }
  if(k==6){
    dt <- dt[order(dt[,"Height"],decreasing = TRUE),]
    text(0.5*(dt[1,"POS1"]+dt[1,"POS2"]),dt[1,"Height"],expression(paste(F[ROH]," = 0.14")),cex=1)
    text(0.5*(dt[2,"POS1"]+dt[2,"POS2"]),dt[2,"Height"],expression(paste(F[ROH]," = 0.20")),cex=1)
    text(0.5*(dt[3,"POS1"]+dt[3,"POS2"]),dt[3,"Height"],expression(paste(F[ROH]," = 0.25")),cex=1)
    text(0.5*(dt[4,"POS1"]+dt[4,"POS2"]),dt[4,"Height"],expression(paste(F[ROH]," = 0.24")),cex=1)
    text(0.5*(dt[5,"POS1"]+dt[5,"POS2"]),dt[5,"Height"],expression(paste(F[ROH]," = 0.28")),cex=1)
    text(0.5*(dt[6,"POS1"]+dt[6,"POS2"]),dt[6,"Height"],expression(paste(F[ROH]," = 0.13")),cex=1)
  }
}
par(op)
dev.off()
