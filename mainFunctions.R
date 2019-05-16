setwd("/home/l.yengo/loic/ID/extreme/Revision_NC/simulations")
dirMap  <- "/home/l.yengo/loic/ID/extreme/bcftools/maps" # can be downloaded here: ftp://ngs.sanger.ac.uk/production/samtools/genetic-map.tgz
dirHaps <- "/home/l.yengo/loic/ID/extreme/Revision_NC/simulations/base_population/phasing/pool-haps"

## Hap must be formatted in R as a list of 3 elements: 
## "SNPID" (which must match the mapfiles in plink format), "H1" and "H2", which matrices with M SNPs (as in mapfiles) and N columns of 0/1
## N is the number of haplotypes (SHAPEIT2 format) = number of phased indviduals.

## Load Error parameters
load("errors-parameters.RData") ## Provided on github

## Sampling recombination breakpoints
breakPoints <- function(pms,pme,rpbMap){
  rce  <- rbinom(length(rpbMap),size=1,prob=rpbMap)
  nbrk <- sum(rce)
  if(nbrk>0){
    l   <- which(rce==1)
    brk <- sort( runif(n=nbrk,min=pms[l],max=pme[l]) )
  }else{
    brk <- NULL
  }
  return(brk)
}

## Induce errors
induceError <- function(x){
  if(x==0){
    y <- ifelse(rbinom(1,size=1,prob=probErrorsGenotypes[1,1])==1,1,2)
  }
  if(x==1){
    y <- ifelse(rbinom(1,size=1,prob=probErrorsGenotypes[2,1])==1,2,0)
  }
  if(x==2){
    y <- ifelse(rbinom(1,size=1,prob=probErrorsGenotypes[3,1])==1,1,0)
  }
  return(y)
}

## Recombine
recomb <- function(x,y,posChip,pms,pme,rpbMap){
  brks <- breakPoints(pms,pme,rpbMap)
  nbks <- length(brks)
  if(nbks>0){
    boundaries <- unique( c(0,brks,max(c(brks,posChip))) )
  }else{
    boundaries <- c(0,max(posChip))
  }
  nsegs  <- length(boundaries)-1
  offset <- rbinom(n=1,size=1,prob=0.5)
  ((1:nsegs) + offset)%%2
  
  z <- rep(NA,length(x))
  for(j in 1:nsegs){
    segment_j    <- which(posChip>=boundaries[j] & posChip<=boundaries[j+1])
    w            <- (j+offset)%%2
    if(w==1){
      z[segment_j] <- x[segment_j]
    }else{
      z[segment_j] <- y[segment_j]
    }
  }
  return(z)
}

NonInbred  <- function(MatingType,posChip,pms,pme,rpbMap,p1h,p2h,m1h,m2h,chr){
  p <- length(posChip)
  ## P and M
  p1   <- rep("p1",p)
  p2   <- rep("p2",p)
  m1   <- rep("m1",p)
  m2   <- rep("m2",p)
  o12a <- recomb(p1,p2,posChip,pms,pme,rpbMap)
  o12b <- recomb(m1,m2,posChip,pms,pme,rpbMap)

  ## Check IBD
  comp  <- o12a==o12b
  nsg   <- length(table(diff(which(comp))))
  
  if(nsg>0){
    ecomp <- c(0,comp,0)
    epos  <- c(posChip[1],posChip,posChip[p])/1e6
    strs  <- epos[which(diff(ecomp)==+1)]
    ends  <- epos[which(diff(ecomp)==-1)]
    segs  <- cbind(CHR=chr,START=strs,END=ends,LENGTHinMB=ends-strs)
    SIBD  <- sum(ends-strs)  
  }else{
    segs  <- NULL
    SIBD  <- 0     
  }
  
  ## Simulate genotypes
  oa  <- rep(NA,p)
  ob  <- rep(NA,p)
  for(i in c("m1","m2","p1","p2")){
    li <- which(o12a==i)
    if(length(li)>0){
      oa[li] <- get(paste0(i,"h"))[li]
    }
    li <- which(o12b==i)
    if(length(li)>0){
      ob[li] <- get(paste0(i,"h"))[li]
    }
  }
  go <- oa + ob
  
  ## Add errors here!
  lambdaPoisson <- errorsChroms[chr,"Mean"]
  numErrorChrom <- rpois(1,lambda = lambdaPoisson)
  who           <- sample(1:p,numErrorChrom)
  goNoised      <- go
  goNoised[who] <- sapply(go[who],induceError)
  
  ## Inbreeding measures
  CL   <- (posChip[p]/1e6)
  FIBD <- SIBD / CL
  
  ## Return
  return( list(genotype=go,
               genotype_noised=goNoised,
               segments=segs,
               numErrorChrom=numErrorChrom,
               fstats=c(nSEG=nsg,SIBDinMB=SIBD,Fibd=FIBD,CL=CL)
  )
  )
}



FirstDegree  <- function(MatingType,posChip,pms,pme,rpbMap,p1h,p2h,m1h,m2h,chr){
  p <- length(posChip)
  ## P and M
  p1 <- rep("p1",p)
  p2 <- rep("p2",p)
  m1 <- rep("m1",p)
  m2 <- rep("m2",p)
  
  if(MatingType=="FS"){
    ## Sib A
    p12a <- recomb(p1,p2,posChip,pms,pme,rpbMap)
    m12a <- recomb(m1,m2,posChip,pms,pme,rpbMap)
    
    ## Sib B
    p12b <- recomb(p1,p2,posChip,pms,pme,rpbMap)
    m12b <- recomb(m1,m2,posChip,pms,pme,rpbMap)
    
    ## A + B Offspring 
    o12a <- recomb(p12a,m12a,posChip,pms,pme,rpbMap)
    o12b <- recomb(p12b,m12b,posChip,pms,pme,rpbMap)
  }
  
  if(MatingType=="PO"){
    ## O = Offspring(P,M)
    p12a <- recomb(p1,p2,posChip,pms,pme,rpbMap)
    m12a <- recomb(m1,m2,posChip,pms,pme,rpbMap)
    
    ## Offspring of O = Offspring(O,P)
    o12a <- recomb(p1,p2,posChip,pms,pme,rpbMap)
    o12b <- recomb(p12a,m12a,posChip,pms,pme,rpbMap)
  }
  
  ## Check IBD
  comp  <- o12a==o12b
  nsg   <- length(table(diff(which(comp))))
 
  if(nsg>0){
   ecomp <- c(0,comp,0)
   epos  <- c(posChip[1],posChip,posChip[p])/1e6
   strs  <- epos[which(diff(ecomp)==+1)]
   ends  <- epos[which(diff(ecomp)==-1)]
   segs  <- cbind(CHR=chr,START=strs,END=ends,LENGTHinMB=ends-strs)
   SIBD  <- sum(ends-strs)  
  }else{
   segs  <- NULL
   SIBD  <- 0     
  }
  
  ## Simulate genotypes
  oa  <- rep(NA,p)
  ob  <- rep(NA,p)
  for(i in c("m1","m2","p1","p2")){
    li <- which(o12a==i)
    if(length(li)>0){
      oa[li] <- get(paste0(i,"h"))[li]
    }
    li <- which(o12b==i)
    if(length(li)>0){
      ob[li] <- get(paste0(i,"h"))[li]
    }
  }
  go <- oa + ob
  
  ## Add errors here!
  lambdaPoisson <- errorsChroms[chr,"Mean"]
  numErrorChrom <- rpois(1,lambda = lambdaPoisson)
  who           <- sample(1:p,numErrorChrom)
  goNoised      <- go
  goNoised[who] <- sapply(go[who],induceError)
  
  ## Inbreeding measures
  CL   <- (posChip[p]/1e6)
  FIBD <- SIBD / CL
  
  ## Return
  return( list(genotype=go,
               genotype_noised=goNoised,
               segments=segs,
               numErrorChrom=numErrorChrom,
               fstats=c(nSEG=nsg,SIBDinMB=SIBD,Fibd=FIBD,CL=CL)
  )
  )
}

SecondDegree <- function(MatingType,posChip,pms,pme,rpbMap,p1h,p2h,m1h,m2h,x1h,x2h,chr){
  p <- length(posChip)
  ## P and M
  p1 <- rep("p1",p); p2 <- rep("p2",p)
  m1 <- rep("m1",p); m2 <- rep("m2",p)
  x1 <- rep("x1",p); x2 <- rep("x2",p)
  
  if(MatingType=="HS"){
    ## H-Sib A: o(p,m)
    p12a <- recomb(p1,p2,posChip,pms,pme,rpbMap)
    m12a <- recomb(m1,m2,posChip,pms,pme,rpbMap)
    
    ## H-Sib B: o(p,x)
    p12b <- recomb(p1,p2,posChip,pms,pme,rpbMap)
    x12b <- recomb(x1,x2,posChip,pms,pme,rpbMap)
    
    ## A + B Offspring: o(A,B) 
    o12a <- recomb(p12a,m12a,posChip,pms,pme,rpbMap)
    o12b <- recomb(p12b,x12b,posChip,pms,pme,rpbMap)
  }
  
  if(MatingType=="AV"){
    ## Sib A
    p12a <- recomb(p1,p2,posChip,pms,pme,rpbMap)
    m12a <- recomb(m1,m2,posChip,pms,pme,rpbMap)
    
    ## Sib B
    p12b <- recomb(p1,p2,posChip,pms,pme,rpbMap)
    m12b <- recomb(m1,m2,posChip,pms,pme,rpbMap)
    
    ## Sib B mate with X to produce Y
    u12x <- recomb(p12b,m12b,posChip,pms,pme,rpbMap)
    v12x <- recomb(x1,x2,posChip,pms,pme,rpbMap)
    
    ## Sib A mate with Y (A is Y uncle/aunt)
    o12a <- recomb(p12a,m12a,posChip,pms,pme,rpbMap)
    o12b <- recomb(u12x,v12x,posChip,pms,pme,rpbMap)
  }
  
  if(MatingType=="GP"){
    ## Sib A
    p12a <- recomb(p1,p2,posChip,pms,pme,rpbMap)
    m12a <- recomb(m1,m2,posChip,pms,pme,rpbMap)
    
    ## Offspring of A and X is Y
    p12b <- recomb(p12a,m12a,posChip,pms,pme,rpbMap)
    x12b <- recomb(x1,x2,posChip,pms,pme,rpbMap)
    
    ## Sib B mate with X to produce Y
    o12a <- recomb(p12b,x12b,posChip,pms,pme,rpbMap)
    o12b <- recomb(p1,p2,posChip,pms,pme,rpbMap)
  }
  
  ## Check IBD
  comp  <- o12a==o12b
  nsg   <- length(table(diff(which(comp))))
  
  if(nsg>0){
   ecomp <- c(0,comp,0)
   epos  <- c(posChip[1],posChip,posChip[p])/1e6
   strs  <- epos[which(diff(ecomp)==+1)]
   ends  <- epos[which(diff(ecomp)==-1)]
   segs  <- cbind(CHR=chr,START=strs,END=ends,LENGTHinMB=ends-strs)
   SIBD  <- sum(ends-strs)
  }else{
   segs <- NULL
   SIBD <-   0
  }
  
  ## Simulate genotypes
  oa  <- rep(NA,p)
  ob  <- rep(NA,p)
  for(i in c("m1","m2","p1","p2","x1","x2")){
    li <- which(o12a==i)
    if(length(li)>0){
      oa[li] <- get(paste0(i,"h"))[li]
    }
    li <- which(o12b==i)
    if(length(li)>0){
      ob[li] <- get(paste0(i,"h"))[li]
    }
  }
  go <- oa + ob
  
  ## Add errors here!
  lambdaPoisson <- errorsChroms[chr,"Mean"]
  numErrorChrom <- rpois(1,lambda = lambdaPoisson)
  who           <- sample(1:p,numErrorChrom)
  goNoised      <- go
  goNoised[who] <- sapply(go[who],induceError)
  
  ## Inbreeding measures
  CL   <- (posChip[p]/1e6)
  FIBD <- SIBD / CL
  
  ## Return
  return( list(genotype=go,
               genotype_noised=goNoised,
               segments=segs,
               numErrorChrom=numErrorChrom,
               fstats=c(nSEG=nsg,SIBDinMB=SIBD,Fibd=FIBD,CL=CL)
  )
  )
}

DFC  <- function(MatingType,posChip,pms,pme,rpbMap,p1h,p2h,m1h,m2h,x1h,x2h,y1h,y2h,chr){
  p <- length(posChip)
  ## P and M
  p1 <- rep("p1",p); p2 <- rep("p2",p)
  m1 <- rep("m1",p); m2 <- rep("m2",p)
  x1 <- rep("x1",p); x2 <- rep("x2",p)
  y1 <- rep("y1",p); y2 <- rep("y2",p)
  
  ## Offspring of p and m
  ## Sib Apm
  p12a <- recomb(p1,p2,posChip,pms,pme,rpbMap)
  m12a <- recomb(m1,m2,posChip,pms,pme,rpbMap)
  
  ## Sib Bpm
  p12b <- recomb(p1,p2,posChip,pms,pme,rpbMap)
  m12b <- recomb(m1,m2,posChip,pms,pme,rpbMap)
  
  ## Offspring of x and y
  ## Sib Axy
  x12a <- recomb(x1,x2,posChip,pms,pme,rpbMap)
  y12a <- recomb(y1,y2,posChip,pms,pme,rpbMap)
  
  ## Sib Bxy
  x12b <- recomb(x1,x2,posChip,pms,pme,rpbMap)
  y12b <- recomb(y1,y2,posChip,pms,pme,rpbMap) 
  
  ## A(p,m) mating with A(x,y)
  ua1  <- recomb(p12a,m12a,posChip,pms,pme,rpbMap)
  ua2  <- recomb(x12a,y12a,posChip,pms,pme,rpbMap)
  
  ## B(p,m) mating with B(x,y)
  ub1  <- recomb(p12b,m12b,posChip,pms,pme,rpbMap)
  ub2  <- recomb(x12b,y12b,posChip,pms,pme,rpbMap)
  
  ## Offspring of O = Offspring(uA,uB)
  o12a <- recomb(ua1,ua2,posChip,pms,pme,rpbMap)
  o12b <- recomb(ub1,ub2,posChip,pms,pme,rpbMap)

  ## Check IBD
  comp  <- o12a==o12b
  nsg   <- length(table(diff(which(comp))))
  
  if(nsg>0){
    ecomp <- c(0,comp,0)
    epos  <- c(posChip[1],posChip,posChip[p])/1e6
    strs  <- epos[which(diff(ecomp)==+1)]
    ends  <- epos[which(diff(ecomp)==-1)]
    segs  <- cbind(CHR=chr,START=strs,END=ends,LENGTHinMB=ends-strs)
    SIBD  <- sum(ends-strs)  
  }else{
    segs  <- NULL
    SIBD  <- 0     
  }
  
  ## Simulate genotypes
  oa  <- rep(NA,p)
  ob  <- rep(NA,p)
  for(i in c("m1","m2","p1","p2","x1","x2","y1","y2")){
    li <- which(o12a==i)
    if(length(li)>0){
      oa[li] <- get(paste0(i,"h"))[li]
    }
    li <- which(o12b==i)
    if(length(li)>0){
      ob[li] <- get(paste0(i,"h"))[li]
    }
  }
  go <- oa + ob
  
  ## Add errors here!
  lambdaPoisson <- errorsChroms[chr,"Mean"]
  numErrorChrom <- rpois(1,lambda = lambdaPoisson)
  who           <- sample(1:p,numErrorChrom)
  goNoised      <- go
  goNoised[who] <- sapply(go[who],induceError)
  
  ## Inbreeding measures
  CL   <- (posChip[p]/1e6)
  FIBD <- SIBD / CL
  
  ## Return
  return( list(genotype=go,
               genotype_noised=goNoised,
               segments=segs,
               numErrorChrom=numErrorChrom,
               fstats=c(nSEG=nsg,SIBDinMB=SIBD,Fibd=FIBD,CL=CL)
               )
          )
}

FirstCousin <- function(MatingType,posChip,pms,pme,rpbMap,p1h,p2h,m1h,m2h,x1h,x2h,y1h,y2h,chr){
  p <- length(posChip)
  ## P and M
  p1 <- rep("p1",p); p2 <- rep("p2",p)
  m1 <- rep("m1",p); m2 <- rep("m2",p)
  x1 <- rep("x1",p); x2 <- rep("x2",p)
  y1 <- rep("y1",p); y2 <- rep("y2",p)
  
  ## Offspring of p and m
  ## Sib Apm
  p12a <- recomb(p1,p2,posChip,pms,pme,rpbMap)
  m12a <- recomb(m1,m2,posChip,pms,pme,rpbMap)
  
  ## Sib Bpm
  p12b <- recomb(p1,p2,posChip,pms,pme,rpbMap)
  m12b <- recomb(m1,m2,posChip,pms,pme,rpbMap)
  
  ## Offspring of A12 and x
  x12a <- recomb(x1,x2,posChip,pms,pme,rpbMap)
  y12a <- recomb(p12a,m12a,posChip,pms,pme,rpbMap)
  
  ## Offspring of B12 and x
  x12b <- recomb(p12b,m12b,posChip,pms,pme,rpbMap)
  y12b <- recomb(y1,y2,posChip,pms,pme,rpbMap) 
  
  ## Offspring of O = Offspring(uA,uB)
  o12a <- recomb(x12a,y12a,posChip,pms,pme,rpbMap)
  o12b <- recomb(x12b,y12b,posChip,pms,pme,rpbMap)
  
  ## Check IBD
  comp  <- o12a==o12b
  nsg   <- length(table(diff(which(comp))))
  
  if(nsg>0){
    ecomp <- c(0,comp,0)
    epos  <- c(posChip[1],posChip,posChip[p])/1e6
    strs  <- epos[which(diff(ecomp)==+1)]
    ends  <- epos[which(diff(ecomp)==-1)]
    segs  <- cbind(CHR=chr,START=strs,END=ends,LENGTHinMB=ends-strs)
    SIBD  <- sum(ends-strs)  
  }else{
    segs  <- NULL
    SIBD  <- 0     
  }
  
  ## Simulate genotypes
  oa  <- rep(NA,p)
  ob  <- rep(NA,p)
  for(i in c("m1","m2","p1","p2","x1","x2","y1","y2")){
    li <- which(o12a==i)
    if(length(li)>0){
      oa[li] <- get(paste0(i,"h"))[li]
    }
    li <- which(o12b==i)
    if(length(li)>0){
      ob[li] <- get(paste0(i,"h"))[li]
    }
  }
  go <- oa + ob
  
  ## Add errors here!
  lambdaPoisson <- errorsChroms[chr,"Mean"]
  numErrorChrom <- rpois(1,lambda = lambdaPoisson)
  who           <- sample(1:p,numErrorChrom)
  goNoised      <- go
  goNoised[who] <- sapply(go[who],induceError)
  
  ## Inbreeding measures
  CL   <- (posChip[p]/1e6)
  FIBD <- SIBD / CL
  
  ## Return
  return( list(genotype=go,
               genotype_noised=goNoised,
               segments=segs,
               numErrorChrom=numErrorChrom,
               fstats=c(nSEG=nsg,SIBDinMB=SIBD,Fibd=FIBD,CL=CL)
  )
  )
} 

## --- ##

simulateNonInbred <- function(indexParent1,indexParent2,MatingType,Chromosome,verbose=FALSE){
  chr <- Chromosome
  ## Read genomic map
  if(verbose){
    cat(paste0("\tRead genomic map from chromosome ",ifelse(chr<10,paste0("0",chr),chr),".\n"))
  }
  
  mapfile <- paste0(dirMap,"/genetic_map_chr",chr,"_combined_b37.txt")
  map     <- read.table(mapfile,header=TRUE)
  rpbMap  <- diff(map[,3]) * 0.01
  #pms,pme  <- map[1:(nrow(map)-1),1] + 0.5*diff(map[,1])
  pms     <- map[1:(nrow(map)-1),1]
  pme     <- map[2:nrow(map),1]
  
  ## Loading haplotypes
  if(verbose){
    cat("\tLoading haplotypes.\n")
  }
  load(paste0(dirHaps,"/chrom",chr,".haps.RData"))
  posChip  <- as.numeric( gsub(paste0("chr",chr,":"),"",haps$ids,fixed=TRUE) )
  
  ## Get haplotypes
  p1h     <- haps$H1[,indexParent1]
  p2h     <- haps$H2[,indexParent1]
  m1h     <- haps$H1[,indexParent2]
  m2h     <- haps$H2[,indexParent2]
  
  results <- NonInbred(MatingType,posChip,pms,pme,rpbMap,p1h,p2h,m1h,m2h,chr)
  return(results)
}


simulateInbredGenotype1 <- function(indexParent1,indexParent2,MatingType,Chromosome,verbose=FALSE){
  chr <- Chromosome
  ## Read genomic map
  if(verbose){
    cat(paste0("\tRead genomic map from chromosome ",ifelse(chr<10,paste0("0",chr),chr),".\n"))
  }
  
  mapfile <- paste0(dirMap,"/genetic_map_chr",chr,"_combined_b37.txt")
  map     <- read.table(mapfile,header=TRUE)
  rpbMap  <- diff(map[,3]) * 0.01
  #posMap  <- map[1:(nrow(map)-1),1] + 0.5*diff(map[,1])
  pms     <- map[1:(nrow(map)-1),1]
  pme     <- map[2:nrow(map),1]
  
  ## Loading haplotypes
  if(verbose){
    cat("\tLoading haplotypes.\n")
  }
  load(paste0(dirHaps,"/chrom",chr,".haps.RData"))
  posChip  <- as.numeric( gsub(paste0("chr",chr,":"),"",haps$ids,fixed=TRUE) )
  
  ## Get haplotypes
  p1h     <- haps$H1[,indexParent1]
  p2h     <- haps$H2[,indexParent1]
  m1h     <- haps$H1[,indexParent2]
  m2h     <- haps$H2[,indexParent2]
  
  results <- FirstDegree(MatingType,posChip,pms,pme,rpbMap,p1h,p2h,m1h,m2h,chr)
  return(results)
}

simulateInbredGenotype2 <- function(indexParent1,indexParent2,indexParent3,MatingType,Chromosome,verbose=FALSE){
  chr <- Chromosome
  ## Read genomic map
  if(verbose){
    cat(paste0("\tRead genomic map from chromosome ",ifelse(chr<10,paste0("0",chr),chr),".\n"))
  }
  mapfile <- paste0(dirMap,"/genetic_map_chr",chr,"_combined_b37.txt")
  map     <- read.table(mapfile,header=TRUE)
  rpbMap  <- diff(map[,3]) * 0.01
  #posMap  <- map[1:(nrow(map)-1),1] + 0.5*diff(map[,1])
  pms     <- map[1:(nrow(map)-1),1]
  pme     <- map[2:nrow(map),1]
  
  ## Loading haplotypes
  if(verbose){
    cat("\tLoading haplotypes.\n")
  }
  load(paste0(dirHaps,"/chrom",chr,".haps.RData"))
  posChip  <- as.numeric( gsub(paste0("chr",chr,":"),"",haps$ids,fixed=TRUE) )
  
  ## Get haplotypes
  p1h     <- haps$H1[,indexParent1]
  p2h     <- haps$H2[,indexParent1]
  m1h     <- haps$H1[,indexParent2]
  m2h     <- haps$H2[,indexParent2]
  x1h     <- haps$H1[,indexParent3]
  x2h     <- haps$H2[,indexParent3]
  
  results <- SecondDegree(MatingType,posChip,pms,pme,rpbMap,p1h,p2h,m1h,m2h,x1h,x2h,chr)
  return(results)
}

simulateInbredGenotype3 <- function(indexParent1,indexParent2,indexParent3,indexParent4,MatingType,Chromosome,verbose=FALSE){
  chr <- Chromosome
  ## Read genomic map
  if(verbose){
    cat(paste0("\tRead genomic map from chromosome ",ifelse(chr<10,paste0("0",chr),chr),".\n"))
  }
  mapfile <- paste0(dirMap,"/genetic_map_chr",chr,"_combined_b37.txt")
  map     <- read.table(mapfile,header=TRUE)
  rpbMap  <- diff(map[,3]) * 0.01
  #posMap  <- map[1:(nrow(map)-1),1] + 0.5*diff(map[,1])
  pms     <- map[1:(nrow(map)-1),1]
  pme     <- map[2:nrow(map),1]
  
  ## Loading haplotypes
  if(verbose){
    cat("\tLoading haplotypes.\n")
  }
  load(paste0(dirHaps,"/chrom",chr,".haps.RData"))
  posChip  <- as.numeric( gsub(paste0("chr",chr,":"),"",haps$ids,fixed=TRUE) )
  
  ## Get haplotypes
  p1h     <- haps$H1[,indexParent1]
  p2h     <- haps$H2[,indexParent1]
  m1h     <- haps$H1[,indexParent2]
  m2h     <- haps$H2[,indexParent2]
  x1h     <- haps$H1[,indexParent3]
  x2h     <- haps$H2[,indexParent3]
  y1h     <- haps$H1[,indexParent4]
  y2h     <- haps$H2[,indexParent4]
  
  results <- DFC(MatingType,posChip,pms,pme,rpbMap,p1h,p2h,m1h,m2h,x1h,x2h,y1h,y2h,chr)
  return(results)
}

simulateInbredGenotype4 <- function(indexParent1,indexParent2,indexParent3,indexParent4,MatingType,Chromosome,verbose=FALSE){
  chr <- Chromosome
  ## Read genomic map
  if(verbose){
    cat(paste0("\tRead genomic map from chromosome ",ifelse(chr<10,paste0("0",chr),chr),".\n"))
  }
  mapfile <- paste0(dirMap,"/genetic_map_chr",chr,"_combined_b37.txt")
  map     <- read.table(mapfile,header=TRUE)
  rpbMap  <- diff(map[,3]) * 0.01
  pms     <- map[1:(nrow(map)-1),1]
  pme     <- map[2:nrow(map),1]
  
  ## Loading haplotypes
  if(verbose){
    cat("\tLoading haplotypes.\n")
  }
  load(paste0(dirHaps,"/chrom",chr,".haps.RData"))
  posChip  <- as.numeric( gsub(paste0("chr",chr,":"),"",haps$ids,fixed=TRUE) )
  
  ## Get haplotypes
  p1h     <- haps$H1[,indexParent1]
  p2h     <- haps$H2[,indexParent1]
  m1h     <- haps$H1[,indexParent2]
  m2h     <- haps$H2[,indexParent2]
  x1h     <- haps$H1[,indexParent3]
  x2h     <- haps$H2[,indexParent3]
  y1h     <- haps$H1[,indexParent4]
  y2h     <- haps$H2[,indexParent4]
  
  results <- FirstCousin(MatingType,posChip,pms,pme,rpbMap,p1h,p2h,m1h,m2h,x1h,x2h,y1h,y2h,chr)
  return(results)
}


