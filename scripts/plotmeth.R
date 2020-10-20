#!/usr/bin/env Rscript

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager',repos='http://cran.us.r-project.org')
if (!requireNamespace('optparse', quietly = TRUE))
  install.packages('optparse',repos='http://cran.us.r-project.org')
if (!requireNamespace('data.table', quietly = TRUE))
  install.packages('data.table',repos='http://cran.us.r-project.org')
if (!requireNamespace('karyoploteR', quietly = TRUE))
  BiocManager::install('karyoploteR')

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(karyoploteR))

#functions for estimation of HMM parameters

ParamEstSeq <- function(DataMatrix,omega) {

  T=ncol(DataMatrix)
  NExp<-nrow(DataMatrix)
  sigmax<-c()
  mi<-c()
  
  for (i in 1:NExp)
  {
    mi[i]<-0
    sigmax[i]<-sd(DataMatrix[i,])^2
  }
  
  smu<-sqrt(omega*sigmax)
  sepsilon<-sqrt((1-omega)*sigmax)
  Results<-list()
  Results$mi<-mi
  Results$smu<-smu
  Results$sepsilon<-sepsilon
  
  Results

}


MukEst <- function(DataMatrix,mw)
{
  NExp<-dim(DataMatrix)[1]
  if (NExp==1)
  {
    muk<-rbind(seq(-1,1,by=0.1))
  }
  if (NExp>1)
  {
    DataMatrix[which(DataMatrix>1)]<-1
    DataMatrix[which(DataMatrix< -1)]<- -1
    
    DataMatrixA<-c()
    for (i in 1:NExp)
    {
      DataMatrixA<-rbind(DataMatrixA,SMA(DataMatrix[i,], n=mw))
    }
    
    DataMatrixA<-DataMatrixA[,2:length(DataMatrixA[1,])]
    
    binsize=0.2
    binVec<-c(seq(-1,-0.2,by=binsize),0,seq(0.2,1,by=binsize))
    binmean<-c(seq(-0.9,-0.3,by=binsize),0,0,seq(0.3,0.9,by=binsize))
    
    DataQuant<-DataMatrixA
    
    for (i in 1:(length(binVec)-1))
    {
      DataQuant[which(DataMatrixA>binVec[i] & DataMatrixA<=binVec[i+1])]<-binmean[i]
    }
    
    muk<-unique(DataQuant,  MARGIN = 2)
    muk<-muk[,-1]
  }
  muk
}


JointSegIn <- function(DataMatrix,muk,mi,smu,sepsilon,Pos,omega,eta,stepeta)
{
  CovPos<-diff(Pos)
  CovPosNorm<-CovPos/stepeta
  etavec<-eta+((1-eta)*exp(log(eta)/CovPosNorm))
  
  NCov<-length(etavec)
  K0<-ncol(muk)
  etav<-log(rep(1,K0)*(1/K0))
  T=ncol(DataMatrix)
  NExp<-nrow(DataMatrix)
  
  TruncCoef<-matrix(0,nrow=NExp,ncol=K0)
  for (jj in 1:NExp)
  {
    for (kk in 1:length(muk))
    {
      TruncCoef[jj,kk]<-sepsilon[jj]*(pnorm(1, mean = muk[jj,kk], sd = sepsilon[jj])-pnorm(-1, mean = muk[jj,kk], sd = sepsilon[jj]))
    }
  }
  
  P<-matrix(data=0,nrow=K0,ncol=(K0*NCov))
  G<-matrix(data=0,nrow=K0,ncol=K0)
  emission<-matrix(data=0,nrow=K0,ncol=T)
  out<-.Fortran("transemisi",as.vector(muk),as.vector(mi),as.double(etavec),as.integer(NCov),as.matrix(DataMatrix),as.integer(K0),as.integer(NExp),as.vector(smu),as.vector(sepsilon),as.integer(T),as.matrix(G),as.matrix(P),as.matrix(emission),as.matrix(TruncCoef))
  
  P<-out[[12]]
  emission<-out[[13]]
  
  psi<-matrix(data=0,nrow=K0,ncol=T)
  path<-c(as.integer(rep(0,T)))
  out2<-.Fortran("bioviterbii",as.vector(etav),as.matrix(P),as.matrix(emission),as.integer(T),as.integer(K0),as.vector(path),as.matrix(psi))
  s<-out2[[6]]
   
  sortResult <- SortState(s)
  TotalPredBreak<-sortResult[[3]]
  TotalPredBreak
}


SortState <- function(s)
{
  l<-1
  seg<-c()
  brek<-c()
  t<-1
  for (k in 1:(length(s)-1))
  {
    if (s[k]!=s[k+1])
    {
      brek[t]<-k
      t<-t+1
      if (length(which(seg==s[k]))==0)
      {
        seg[l]<-s[k]
        l<-l+1
      }
    }
  }
  brek<-c(0,brek,length(s))
  if (length(which(seg==s[length(s)]))==0)
  {
    seg<-c(seg,s[length(s)])
  }
  
  s0<-c()
  for (k in 1:length(seg))
  {
    s0[which(s==seg[k])]<-k
  }
  
  SortResult<-list()
  SortResult[[1]]<-s0
  SortResult[[2]]<-seg
  SortResult[[3]]<-brek
  SortResult
}


SegResults <- function(DataSeq,TotalPredBreak)
{
  
  regioncall<-c()
  TotalPred<-c()
  NExp<-nrow(DataSeq)
  for (j in 1:NExp)
  {
    s<-rep(0,ncol(DataSeq))
    for (i in 1:(length(TotalPredBreak)-1))
    {
      s[(TotalPredBreak[i]+1):TotalPredBreak[i+1]]<-median(DataSeq[j,(TotalPredBreak[i]+1):TotalPredBreak[i+1]])
    }
    TotalPred<-rbind(TotalPred,s)
  }
  
  Result<-TotalPred
  Result
}


FilterSeg <- function(TotalPredBreak,FW)
{
  controllength<-diff(TotalPredBreak)
  
  indF<-which(controllength<=FW)
  if (length(indF)!=0)
  {
    if (indF[1]==1)
    {
      indF[1]<-2
      indF<-unique(indF)
      TotalPredBreak1 <- TotalPredBreak[-(indF)]
    }
    if (indF[1]!=1)
    {
      TotalPredBreak1 <- TotalPredBreak[-(indF)]
    }
  }
  if (length(indF)==0)
  {
    TotalPredBreak1<-TotalPredBreak
  }
  TotalPredBreak1
}


#initial parameters

omega <- 0.1
eta <- 1e-5
stepeta <- 10e5
FW <- 1
mw <- 1

defaultW <- getOption("warn")
options(warn = -1)

option_list = list(
  make_option(c('-p', '--haplotype1'), action='store', type='character', help='.tsv file containing methylation frequencies for haplotype 1 [required]'),
  make_option(c('-q', '--haplotype2'), action='store', type='character', help='.tsv file containing methylation frequencies for haplotype 2', default=NULL),
  make_option(c('-o', '--outputdir'), action='store', type='character', help='output directory [required]'),
  make_option(c('-b', '--bed'), action='store', type='character', help='.bed file with regions to restrict the analysis to [required]'),
  make_option(c('-r', '--release'), action='store', type='character', help='genome release [hg38]', default = 'hg38')
)

opt = parse_args(OptionParser(option_list=option_list))

#get .Rscript path

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- '--file='
script.name <- sub(file.arg.name, '', initial.options[grep(file.arg.name, initial.options)])
script.basename <- normalizePath(dirname(script.name))
FortranLibrary<-file.path(script.basename, 'JointSLMTruncl.so')

if (!= file.exists(fortranin)) {

  now<-Sys.time()
  stop('[',now,'][Error] Could not find ', FortranLibrary, ' file. Did Fortran compilation fail?')

} else {

  dyn.load(FortranLibrary)

}

now<-Sys.time()
message('[',now,'][Message] Reading input BED file')
BED<-fread(file.path(opt$bed), sep='\t', header=FALSE)

now<-Sys.time()
message('[',now,'][Message] Reading methylation frequency files')
M1<-fread(file.path(opt$haplotype1), sep='\t', header=TRUE)

if (! is.null(opt$haplotype2)) {

  M2<-fread(file.path(opt$haplotype2), sep='\t', header=TRUE)

} else {

  M2<-NULL

}

#create output folder if it does not exist

dir.create(file.path(opt$outputdir), showWarnings = FALSE)

#process the meth frequencies based on BED lines

for (row in 1:nrow(BED)) {

  chrom<-BED$V1[row]

  if (isFALSE(startsWith(chrom, 'chr'))) {

    subchrom<-paste0('chr', chrom)

  } else {

    subcrom<-chrom

  }

  start<-BED$V2[row]
  end<-BED$V3[row]
  region<-paste0(chrom,':',start,'-',end)

  now<-Sys.time()
  message('[',now,'][Message] Processing ', region)

  #grep region

  subM1<-data.table(M1[grep(region,M1$chromosome),])

  if (nrow(subM1) == 0) {

    now<-Sys.time()
    stop('[',now,'][Error] Missing region in haplotype 1 methylation frequency file')

  }

  if (! is.null( M2)) {

    subM2<-data.table(M2[grep(region,M2$chromosome),])

    if (nrow(subM2) == 0) {

      now<-Sys.time()
      stop('[',now,'][Error] Missing region in haplotype 2 methylation frequency file')

    }

  } else {

    subM2<-NULL

  }


  #haplotype1, merge if multiple chromosomes in the same matrix

  if (length(unique(subM1$chromosome)) != 1) {

    M1res <- subM1[,.(chromosome=region,start=unique(start), end=unique(start), methylated_frequency=weighted.mean(x=methylated_frequency,w=c(called_sites/num_cpgs_in_group)/sum(c(called_sites/num_cpgs_in_group)))),by=start]

  } else {

    M1res<-subM1[,.(chromosome=region, start=start, end=start, methylated_frequency=methylated_frequency), by=start]

  }

  #calculate segmentation

  now<-Sys.time()
  message('[',now,'][Message] Performing segmentation on haplotype 1 ', region)


  ParamList <- ParamEstSeq(rbind(as.numeric(M1res$methylated_frequency)),omega)
  mi <- ParamList$mi
  smu <- ParamList$smu
  sepsilon <- ParamList$sepsilon
  muk <- MukEst(rbind(as.numeric(M1res$methylated_frequency)),mw)

  calls <- JointSegIn(rbind(as.numeric(M1res$methylated_frequency)),muk,mi,smu,sepsilon,as.numeric(M1res$start),omega,eta,stepeta)
  calls <- FilterSeg(calls,FW)
  callsh1 <- SegResults(rbind(as.numeric(M1res$methylated_frequency)),calls)
  M1res$segmentation<- as.numeric(t(callsh1)[,1])

  M1pos<-rescale(M1res$start, c(start,end))
  M1methvals<-data.frame(chromosome=subchrom, start=M1pos, end=M1pos, y=M1res$methylated_frequency)
  M1methseg<-data.frame(chromosome=subchrom, start=M1pos, end=M1pos, y=M1res$segmentation)

  #convert to GRanges for visualization

  GRMV1<-makeGRangesFromDataFrame(M1methvals, keep.extra.columns = TRUE)
  GRMS1<-makeGRangesFromDataFrame(M1methseg, keep.extra.columns = TRUE)

  #haplotype2, merge if multiple chromosomes in the same matrix

  if (! is.null(subM2)) {

    if (length(unique(subM2$chromosome)) != 1) {

      M2res <- subM2[,.(chromosome=region, start=unique(start), end=unique(start) methylated_frequency=weighted.mean(x=methylated_frequency,w=c(called_sites/num_cpgs_in_group)/sum(c(called_sites/num_cpgs_in_group)))),by=start]

    } else {

      M2res <- subM2[,.(chromosome=region, start=start, end=start, methylated_frequency=methylated_frequency),by=start]

    }

    #calculate segmentation

    ParamList <- ParamEstSeq(rbind(as.numeric(M2res$methylated_frequency)),omega)
    mi <- ParamList$mi
    smu <- ParamList$smu
    sepsilon <- ParamList$sepsilon
    muk <- MukEst(rbind(as.numeric(M2res$methylated_frequency)),mw)

    calls <- JointSegIn(rbind(as.numeric(M2res$methylated_frequency)),muk,mi,smu,sepsilon,as.numeric(M2res$start),omega,eta,stepeta)
    calls <- FilterSeg(calls,FW)
    callsh2 <- SegResults(rbind(as.numeric(M2res$methylated_frequency)),calls)
    M2res$segmentation<- as.numeric(t(callsh2)[,1])

    M2pos<-rescale(M2res$start, c(start,end))
    M2methvals<-data.frame(chromosome=subchrom, start=M2pos, end=M2pos, y=M2res$methylated_frequency)
    M2methseg<-data.frame(chromosome=subchrom, start=M2pos, end=M2pos, y=M2res$segmentation)

    #convert to GRanges for visualization

    GRMV2<-makeGRangesFromDataFrame(M2methvals, keep.extra.columns = TRUE)
    GRMS2<-makeGRangesFromDataFrame(M2methseg, keep.extra.columns = TRUE)

  } else {

    PLOT2<-NULL

  }

  zoom.region <- toGRanges(data.frame(subchrom, start-10, end+10))

  #plot

  pdf(file.path(opt$outputdir, paste0(region, '.pdf')), height=10, width=15)

  if (!is.null(PLOT2)) {

    kp <- plotKaryotype(genome=opt$release, plot.type = 2, zoom=zoom.region)
    kpDataBackground(kp, data.panel=1, col='white', r0=0.1, r1=0.9)
    kpDataBackground(kp, data.panel = 2,col='white',r0=0.1, r1=0.9)
    kpAddBaseNumbers(kp, add.units=TRUE,minor.ticks=TRUE,minor.tick.dist=10,major.ticks=TRUE, tick.dist=100)
    kpAxis(kp, data.panel=1,r0=0.1, r1=0.9)
    kpAxis(kp, data.panel=2,r0=0.1, r1=0.9)
    kpPoints(kp, data.panel=1, data=GRMV1, cex=.2,r0=0.1, r1=0.9)
    kpPoints(kp, data.panel=2, data=GRMV2, cex=.2,r0=0.1, r1=0.9)
    kpLines(kp, data.panel=1, data=GRMS1, col='darkred', lwd = 1.5,r0=0.1, r1=0.9)
    kpLines(kp, data.panel=2, data=GRMS2,col='darkblue', lwd = 1.5,r0=0.1, r1=0.9)
    kpAddMainTitle(kp, paste0('Methylation profile of ', region))
    kpAddLabels(kp, data.panel=1,labels = expression('Methylation frequency '['(hap1)']), r0=0.85, r1=1.05,cex=.6,label.margin = .07, srt=90)
    kpAddLabels(kp, data.panel=2,labels = expression('Methylation frequency '['(hap2)']), r0=-0.05, r1=0.15,cex=.6,label.margin = .07, srt=90)
  
  } else {

    kp <- plotKaryotype(genome=opt$release, plot.type = 1, zoom=zoom.region)
    kpDataBackground(kp, data.panel=1, col='white', r0=0.1, r1=0.9)
    kpAddBaseNumbers(kp, add.units=TRUE,minor.ticks=TRUE,minor.tick.dist=10,major.ticks=TRUE, tick.dist=100)
    kpAxis(kp, data.panel=1,r0=0.1, r1=0.9)
    kpPoints(kp, data.panel=1, data=GRMV1, cex=.2,r0=0.1, r1=0.9)
    kpLines(kp, data.panel=1, data=GRMS1, col='darkred', lwd = 1.5,r0=0.1, r1=0.9)
    kpAddMainTitle(kp, paste0('Methylation profile of ', region))
    kpAddLabels(kp, data.panel=1,labels = expression('Methylation frequency '['(hap1)']), r0=0.75, r1=1.05,cex=.8,label.margin = .07, srt=90)
  
  }

  dev.off()

}

options(warn = defaultW)