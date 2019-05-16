## This script can be run in command line using Rscript as
## "Rscript simulation.R MatingType", where Mating type implemented are listed below
## PO: Parent-Offspring
## FS: Fullsibs
## HS: Halfsibs
## AV: Avuncular (i.e. uncle-niece or aunt-nephew)
## GP: Grandparent-grandchild
## DC: Double-first cousin
## FC: First-cousin
## UN: Unrelated parents
## Part of this script runs PLink (v1.9) so requires to specify a path


## The file arcihtecture of "pathTorunSimulation" must contain directories named mapfiles, rohs, fstats and segments
pathTorunSimulation <- ""
plinkPath           <- "/home/l.yengo/exe/plink"

source("mainFunctions.R")
options    <- commandArgs(trailingOnly = TRUE)
MatingType <- options[1]
SimID      <- as.numeric(options[2])

baseROHcmd <- "--homozyg --homozyg-density 50 --homozyg-gap 1000 --homozyg-kb 1500 --homozyg-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-snp 50 --out"

if(MatingType%in%c("FS","PO","HS","AV","GP","DC","UN","FC")){
  ## Size of the pool
  Ne <- 972
  if(MatingType=="UN"){
    ## First degree relatives mating
    indexParents <- sample(1:Ne ,2)
    indexParent1 <- indexParents[1]
    indexParent2 <- indexParents[2]
    simFun       <- function(chr){
      simulateNonInbred(indexParent1,indexParent2,MatingType,Chromosome=chr)
    }
    ResultsAllChroms <- lapply(1:22,simFun)
  }
  if(MatingType%in%c("FS","PO")){
    ## First degree relatives mating
    indexParents <- sample(1:Ne ,2)
    indexParent1 <- indexParents[1]
    indexParent2 <- indexParents[2]
    simFun       <- function(chr){
      simulateInbredGenotype1(indexParent1,indexParent2,MatingType,Chromosome=chr)
    }
    ResultsAllChroms <- lapply(1:22,simFun)
  }
  if(MatingType%in%c("HS","AV","GP")){
    ## Second degree relatives mating
    indexParents <- sample(1:Ne ,3)
    indexParent1 <- indexParents[1]
    indexParent2 <- indexParents[2]
    indexParent3 <- indexParents[3]
    simFun       <- function(chr){
      simulateInbredGenotype2(indexParent1,indexParent2,indexParent3,MatingType,Chromosome=chr)
    }
    ResultsAllChroms <- lapply(1:22,simFun)
  }
  if(MatingType=="DC"){
    ## Second degree relatives mating
    indexParents <- sample(1:Ne ,4)
    indexParent1 <- indexParents[1]
    indexParent2 <- indexParents[2]
    indexParent3 <- indexParents[3]
    indexParent4 <- indexParents[4]
    simFun       <- function(chr){
      simulateInbredGenotype3(indexParent1,indexParent2,indexParent3,indexParent4,MatingType,Chromosome=chr)
    }
    ResultsAllChroms <- lapply(1:22,simFun)
  }
  if(MatingType=="FC"){
    ## Second degree relatives mating
    indexParents <- sample(1:Ne ,4)
    indexParent1 <- indexParents[1]
    indexParent2 <- indexParents[2]
    indexParent3 <- indexParents[3]
    indexParent4 <- indexParents[4]
    simFun       <- function(chr){
      simulateInbredGenotype4(indexParent1,indexParent2,indexParent3,indexParent4,MatingType,Chromosome=chr)
    }
    ResultsAllChroms <- lapply(1:22,simFun)
  }
  
  ## Output results
  iid <- ifelse(SimID<10,paste0("00",SimID),ifelse(SimID<100,paste0("0",SimID),SimID))
  IID <- paste0(MatingType,iid)
  IDt <- paste0(MatingType,iid,"T")
  IDn <- paste0(MatingType,iid,"N")
  grp <- c("A\tA","A\tG","G\tG")
  
  segs   <- NULL
  fstats <- NULL
  genoRd <- NULL
  roh    <- NULL
  for(chr in 1:22){
    ## write plink file
    gen <- rbind(grp[1+ResultsAllChroms[[chr]]$genotype],grp[1+ResultsAllChroms[[chr]]$genotype_noised])
    ped <- cbind(FID=c(IDt,IDn),IID=c(IDt,IDn),PID=c(0,0),MID=c(0,0),SEX=c(1,1),PHENO=c(0,0),gen)
    ped <- c(paste(ped[1,],collapse = "\t"),paste(ped[2,],collapse = "\t"))
    pedfilename <- paste0("pedfiles/",IID,"_chrom",chr,".ped")
    write(ped,pedfilename)
    
    ## Run ROH
    cmd <- paste0(plinkPath," --silent --ped ",pedfilename," --map ",paste0("mapfiles/chrom",chr,".map")," ",baseROHcmd,paste0(" rohs/",IID,"_chrom",chr))
    system(cmd)
    roh     <- rbind(roh,read.table(paste0("rohs/",IID,"_chrom",chr,".hom"),header=TRUE,stringsAsFactors=FALSE))
     
    ## compile segments
    segs    <- rbind(segs,ResultsAllChroms[[chr]]$segments)
    genoRd  <- c(genoRd,list(cbind(True=ResultsAllChroms[[chr]]$genotype,Noised=ResultsAllChroms[[chr]]$genotype_noised)))
    fstats  <- rbind(fstats,c(ResultsAllChroms[[chr]]$fstats,numError=ResultsAllChroms[[chr]]$numErrorChrom))
    
    ## clean files
    system(paste0("rm ",pedfilename))
    system(paste0("rm rohs/",IID,"_chrom",chr,".hom"))
    system(paste0("rm rohs/",IID,"_chrom",chr,".hom.indiv"))
    system(paste0("rm rohs/",IID,"_chrom",chr,".hom.summary"))
    system(paste0("rm rohs/",IID,"_chrom",chr,".log"))
  }
  save(list=c("segs","roh"),file=paste0("segments/",IID,".RData"))
  save(list="fstats",file=paste0("fstats/",IID,".RData"))
  ## save(list="genoRd",file=paste0("genoRData/",IID,".RData"))
}else{
  cat("\tIncorrect mating type! Must be UN, FS, PO, HS, AV, GP, FC or DC.\n")
  cat("\tPlease try again.\n\n")
}

