#!/usr/bin/env Rscript

if (!requireNamespace('optparse', quietly = TRUE))
  install.packages('optparse',repos='http://cran.us.r-project.org')
if (!requireNamespace('data.table', quietly = TRUE))
  install.packages('data.table',repos='http://cran.us.r-project.org')
if (!requireNamespace('scales', quietly = TRUE))
  install.packages('scales',repos='http://cran.us.r-project.org')
if (!requireNamespace('changepoint', quietly = TRUE))
  install.packages('changepoint',repos='http://cran.us.r-project.org')
if (!requireNamespace('EnvCpt', quietly = TRUE))
  install.packages('EnvCpt',repos='http://cran.us.r-project.org')
if (!requireNamespace('tseries', quietly = TRUE))
  install.packages('tseries',repos='http://cran.us.r-project.org')
if (!requireNamespace('gtools', quietly = TRUE))
  install.packages('gtools',repos='http://cran.us.r-project.org')

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(changepoint))
suppressPackageStartupMessages(library(EnvCpt))
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


modelc<-"meancpt"

now<-Sys.time()
message('[',now,'][Message] Reading input BED file')

#read BED. This should be the same given as input to TREADMILL RACE
BED<-fread(file.path(opt$bed), sep='\t', header=FALSE)

now<-Sys.time()
message('[',now,'][Message] Reading methylation frequency files')

#read frequency TSV. This should have the same format of the TSV from callmethylation.sh script
M<-fread(file.path(opt$frequencies), sep='\t', header=TRUE)

#create output folder if it does not exist
dir.create(file.path(opt$outputdir), showWarnings = FALSE)

if (! is.null(opt$overwrite)) {

  F<-fread(file.path(opt$overwrite), sep='\t', header=TRUE)

}

penaltyfile<-list()

#process the meth frequencies based on BED lines
for (row in 1:nrow(BED)) {

  penaltyregion<-list()

  chrom<-BED$V1[row]

  if (isFALSE(startsWith(chrom, 'chr'))) {

    subchrom<-paste0('chr', chrom)

  } else {

    subchrom<-chrom

  }

  start<-BED$V2[row]
  end<-BED$V3[row]
  region<-paste0(chrom,':',start,'-',end)

  now<-Sys.time()
  message('[',now,'][Message] Processing ', region)

  #grep region

  subM<-data.table(M[grep(region,M$chrom),])
  flanking<-as.numeric(unlist(strsplit(subM$chrom[1],"_"))[7])
  start<-start-flanking
  end<-end+flanking

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

      minseglen<-100

    }

    penaltyregion[[counter]]<-cbind(region,allelename,modelc,minseglen)

    cpts<-envcpt(M1res$methylated_frequency, minseglen=minseglen)
    model<-cpts[[modelc]]

    now<-Sys.time()
    message('[',now,'][Message] Plotting')

    indices<-seq(0, length(M1res$start),length(M1res$start)/(5-1))
    numbers<-seq(min(M1res$start),max(M1res$start),(max(M1res$start)-min(M1res$start))/(5-1))

    plot(cpts, xaxt = "n", main="Changepoint analysis (overview)")
    #axis(1,at=indices, labels=round(numbers))
    plot(model, ylab="Methylation frequency", xaxt = "n", main=paste0("Methylation profile of ", allelename), ylim=c(-0.1,1.1))
    #axis(1,at=indices, labels=round(numbers))
    text(indices, -0.1, round(numbers))
    axis(1, at=indices, NA, tck=.01)

    counter<-counter+1

  }

  dev.off()

  regtab<-do.call(rbind,penaltyregion)
  penaltyfile[[region]]<-regtab

}

alltab<-data.table(do.call(rbind,penaltyfile))
colnames(alltab)<-c("region", "allele", "applied_model", "minseglen")

fwrite(alltab, file=file.path(opt$outputdir, "changepoints.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names = FALSE)

now<-Sys.time()
message('[',now,'][Message] Done')

options(warn = defaultW)