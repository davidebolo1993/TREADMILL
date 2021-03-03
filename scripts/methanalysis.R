#!/usr/bin/env Rscript

if (!requireNamespace('optparse', quietly = TRUE))
  install.packages('optparse',repos='http://cran.us.r-project.org')
if (!requireNamespace('data.table', quietly = TRUE))
  install.packages('data.table',repos='http://cran.us.r-project.org')
if (!requireNamespace('changepoint', quietly = TRUE))
  install.packages('changepoint',repos='http://cran.us.r-project.org')
if (!requireNamespace('tseries', quietly = TRUE))
  install.packages('tseries',repos='http://cran.us.r-project.org')
if (!requireNamespace('gtools', quietly = TRUE))
  install.packages('gtools',repos='http://cran.us.r-project.org')
if (!requireNamespace('scales', quietly = TRUE))
  install.packages('scales',repos='http://cran.us.r-project.org')

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(changepoint))
suppressPackageStartupMessages(library(tseries))
suppressPackageStartupMessages(library(gtools))


#functions that will be used later

defaultW <- getOption("warn")
options(warn = -1)

option_list = list(
  make_option(c('-f', '--frequencies'), action='store', type='character', help='.tsv file containing methylation frequencies as from TREADMILL/scripts/callmethylation.sh [required]'),
  make_option(c('-o', '--outputdir'), action='store', type='character', help='output directory [required]'),
  make_option(c('-b', '--bed'), action='store', type='character', help='.bed file with one or more regions to restrict the analysis to [required]'),
  make_option(c('-x', '--overwrite'), action='store', type='character', help='.tsv file with different minseglen that users want to apply to segmentation', default=NULL)
  )

opt = parse_args(OptionParser(option_list=option_list))


now<-Sys.time()
message('[',now,'][Message] Reading input BED file')

#read BED. This msut be the same input given to TREADMILL RACE
BED<-fread(file.path(opt$bed), sep='\t', header=FALSE)

now<-Sys.time()
message('[',now,'][Message] Reading methylation frequency files')

#read frequency TSV. This must have the same format of the TSV from callmethylation.sh script
M<-fread(file.path(opt$frequencies), sep='\t', header=TRUE)

#create output folder if it does not exist
dir.create(file.path(opt$outputdir), showWarnings = FALSE)

if (! is.null(opt$overwrite)) {

  F<-fread(file.path(opt$overwrite), sep='\t', header=TRUE)

}

overviewfile<-list()

#process the meth frequencies based on BED lines
for (row in 1:nrow(BED)) {

  overviewregion<-list()

  chrom<-BED$V1[row]
  start<-BED$V2[row]
  end<-BED$V3[row]
  region<-paste0(chrom,':',start,'-',end)

  now<-Sys.time()
  message('[',now,'][Message] Processing ', region)

  #grep region

  subM<-data.table(M[grep(region,M$chrom),])
  flanking<-as.numeric(unlist(strsplit(subM$chrom[1],"_"))[7])

  #stop if region not in table

  if (nrow(subM) == 0) {

    now<-Sys.time()
    error('[',now,'][Error] Missing region in methylation frequency file')

  }

  #split by allele

  subs<-split(subM, by="allele_name")
  subs<-setNames(lapply(sort(names(subs)), FUN = function(n) subs[[n]]), mixedsort(names(subs))) #sort so that the first 2 are main genotype and others are subgroups

  #iterate over alleles
  counter<-1
  region_to_plot<-gsub(":", "_", region)
  pdf(file.path(opt$outputdir, paste0(region_to_plot, ".pdf")), height=15, width=20, onefile=FALSE)
  layout(matrix(1:(length(subs)*2), length(subs), 2, byrow = TRUE), widths=c(2,4))

  for (subM1 in subs) {

    allelename<-paste0("allele #",counter)
    mergedM1<-list()
    splittedM1<-split(subM1, by="chrom")
    last<-as.numeric(unlist(strsplit(unique(splittedM1[[length(splittedM1)]]$chrom),"_"))[9])

    #this allows to normalize coordinates among different chromosomes having different lengths

    for (tab in splittedM1) {

      chrom_long<-unique(tab$chrom)
      ie<-as.numeric(unlist(strsplit(chrom_long, "_"))[9])
      toadd<-last-ie
      index_to_mod<-which(tab$pos > ie)
      tab$pos[index_to_mod] <- tab$pos[index_to_mod] + toadd
      mergedM1[[chrom_long]]<-tab

    }

    subM1<-do.call(rbind,mergedM1)

    if (length(unique(subM1$chrom)) != 1) {

      M1res <- subM1[,.(chromosome=region,start=unique(pos), end=unique(pos), methylated_frequency=weighted.mean(x=modification_frequency,w=coverage)),by=pos]
      M1res<-M1res[complete.cases(M1res), ]
      setorderv(M1res, c("pos"), c(1))

    } else {

      M1res<-subM1[,.(chromosome=region, start=pos, end=pos, methylated_frequency=modification_frequency), by=pos]

    }

    #calculate p for trend Stationarity

    t1<-kpss.test(M1res$methylated_frequency, null="Trend")
    kpss1_2<-paste0("p = ", t1$p.value)
    kpss1_3<-paste0("KPSS trend = ", t1$statistic)
    now<-Sys.time()
    message('[',now,'][Message] KPSS Test for Trend Stationarity on ', allelename, ', region ', region, ': ', kpss1_2, '. ', kpss1_3)
    message('[',now,'][Message] Performing changepoint detection analysis on ', allelename, ', region ', region)

    if (! is.null(opt$overwrite)) {

      subF<-data.table(F[grep(region,F$region),])
      subFA<-data.table(F[grep(allelename,F$allele),])
      minseglen<-as.numeric(subFA$minseglen)

    } else {

      minseglen<-20

    }

    overviewregion[[counter]]<-cbind(region,allelename,minseglen)

    segmentate<-cpt.meanvar(M1res$methylated_frequency,method="PELT",penalty = "MBIC", minseglen = minseglen)

    indexes<-c(0, segmentate@cpts)
    segments<-rep(segmentate@param.est$mean,diff(indexes))

    groups<-rep(1,length(segments))
    num<-segments[1]
    thisgroup<-1

    for (j in 2:length(segments)) {

      if (segments[j] != num) {

        thisgroup<-thisgroup+1
        num<-segments[j]

      }

      groups[j]<-thisgroup

    }

    #plot methylation of each group

    plot(segmentate@data.set,log2(M1res$start),pch=20,col=scales::alpha(rainbow(length(unique(groups)))[as.numeric(groups)],0.5), xlab="",ylab="",xaxt = "n", ylim=c(round(min(log2(M1res$start)))-1, round(max(log2(M1res$start)))+1), xlim=c(0,1))
    text(c(0,0.2,0.4,0.6,0.8,1.0),round(min(log2(M1res$start)))-1,c(0,0.2,0.4,0.6,0.8,1.0))
    axis(1, at=c(0,0.2,0.4,0.6,0.8,1.0), NA, tck=.01)
    title(xlab="Metylation frequencies", line=0.3)
    title(main=paste0("Methylation groups in ", allelename), line=0.3)
    mtext("Index (decoy reference, log2)",side=4,line=0.3,cex=.7)

    #plot segmentation

    indices<-seq(1, length(M1res$start),length.out=5)
    numbers<-M1res$start[round(indices)]

    plot(segmentate, xaxt = "n", ylim=c(-0.1,1.1), xlab="",  ylab="", main="")
    text(indices, -0.1, numbers)
    axis(1, at=indices, NA, tck=.01)
    title(xlab="Index (decoy reference)", line=0.3)
    title(main=paste0("Methylation profile of ", allelename), line=0.3)
    mtext("Methylation frequencies",side=4,line=0.3,cex=.7)


    counter<-counter+1

  }

  dev.off()

  regtab<-do.call(rbind,overviewregion)
  overviewfile[[region]]<-regtab

}

alltab<-data.table(do.call(rbind,overviewfile))
colnames(alltab)<-c("region", "allele", "minseglen")

fwrite(alltab, file=file.path(opt$outputdir, "overwrite.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names = FALSE)

now<-Sys.time()
message('[',now,'][Message] Done')

options(warn = defaultW)