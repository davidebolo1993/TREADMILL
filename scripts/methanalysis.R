#!/usr/bin/env Rscript

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager',repos='http://cran.us.r-project.org')
if (!requireNamespace('optparse', quietly = TRUE))
  install.packages('optparse',repos='http://cran.us.r-project.org')
if (!requireNamespace('data.table', quietly = TRUE))
  install.packages('data.table',repos='http://cran.us.r-project.org')
if (!requireNamespace('scales', quietly = TRUE))
  install.packages('scales',repos='http://cran.us.r-project.org')
if (!requireNamespace('changepoint', quietly = TRUE))
  install.packages('changepoint',repos='http://cran.us.r-project.org')
if (!requireNamespace('tseries', quietly = TRUE))
  install.packages('tseries',repos='http://cran.us.r-project.org')
if (!requireNamespace('gtools', quietly = TRUE))
  install.packages('gtools',repos='http://cran.us.r-project.org')
if (!requireNamespace('RColorBrewer', quietly = TRUE))
  install.packages('RColorBrewer',repos='http://cran.us.r-project.org')
if (!requireNamespace('grid', quietly = TRUE))
  install.packages('grid',repos='http://cran.us.r-project.org')
if (!requireNamespace('ggplot2', quietly = TRUE))
  install.packages('ggplot2',repos='http://cran.us.r-project.org')

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(karyoploteR))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(changepoint))
suppressPackageStartupMessages(library(tseries))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))

#functions that will be used later

cptfn <- function(data, pen) { ##https://rpubs.com/richkt/269908
  ans <- cpt.mean(data, test.stat="Normal", method = "PELT", penalty = "Manual", pen.value = pen) 
  length(cpts(ans)) +1
}

define_region <- function(row, col) {
  viewport(layout.pos.row = row, layout.pos.col = col)
}

defaultW <- getOption("warn")
options(warn = -1)

option_list = list(
  make_option(c('-f', '--frequencies'), action='store', type='character', help='.tsv file containing methylation frequencies as from TREADMILL/scripts/callmethylation.sh [required]'),
  make_option(c('-o', '--outputdir'), action='store', type='character', help='output directory [required]'),
  make_option(c('-b', '--bed'), action='store', type='character', help='.bed file with one or more regions to restrict the analysis to [required]'),
  make_option(c('-x', '--overwrite'), action='store', type='character', help='.tsv file with penalties that users want to apply to segmentation', default=NULL)
  )

opt = parse_args(OptionParser(option_list=option_list))

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
  col_vector<- brewer.pal(n = 8, name = "Dark2") #no more than 8 elements, but should be enough
  region_to_plot<-gsub(":", "_", region)
  pdf(file.path(opt$outputdir, paste0(region_to_plot, ".pdf")), height=10, width=20, onefile=FALSE)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(length(subs), 20)))


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

    pen.vals <- seq(0, 10,.1)
    elbowplotData <- unlist(lapply(pen.vals, function(p) cptfn(data = M1res$methylated_frequency, pen = p)))
    elbowframe<-data.frame(x=pen.vals,y=elbowplotData, stringsAsFactors = FALSE)

    allpenalties<-pen.vals[which(diff(elbowplotData) == -1)]

    #We use first non-large step we see in the elbowplot Data (-1)

    if (! is.null(opt$overwrite)) {

      subF<-data.table(F[grep(region,F$region),])
      subFA<-data.table(F[grep(allelename,F$allele),])
      penalty<-subFA$applied_penalty

    } else {

      penalty<-pen.vals[which(diff(elbowplotData) == -1)[1]]

    }

    penaltyregion[[counter]]<-cbind(region,allelename,penalty,paste(allpenalties,collapse=","))

    #Compute changepoint

    cptm_CP <- cpt.mean(M1res$methylated_frequency, penalty='Manual',pen.value=penalty,method='PELT')
    indexes<-c(0,cptm_CP@cpts)
    vals<-cptm_CP@param.est$mean
    svals<-rep(vals,diff(indexes))

    #tables for plotting

    M1pos<-rescale(M1res$start, c(start,end))
    M1methvals<-data.frame(chromosome=subchrom, start=M1pos, end=M1pos, y=M1res$methylated_frequency, z=svals)

    now<-Sys.time()
    message('[',now,'][Message] Plotting')

    p1<-ggplot(M1methvals) + geom_point(aes(x=start,y=y)) + geom_line(aes(x=start,y=z),col=col_vector[counter]) + theme_classic() + scale_y_continuous(limits=c(0,1)) + labs(x="Genomic coordinates (rescaled)", y= "Methylation frequency") + ggtitle(paste0("Methylation profile of ", allelename)) + theme(plot.title=element_text(hjust=0.5))
    p2<-ggplot(elbowframe, aes(x=x,y=y)) + geom_point() + geom_line() + theme_bw() + labs(x="PELT penalty parameter", y= "# cpts") + geom_vline(xintercept=penalty, col=col_vector[counter]) + geom_hline(yintercept=elbowplotData[which(pen.vals==penalty)], col=col_vector[counter]) + annotate(geom="text", x=9,y=max(elbowplotData)-10, label=paste0("# cpts = ", elbowplotData[which(pen.vals==penalty)])) +  annotate(geom="text", x=5,y=max(elbowplotData)-10, label=paste0("penalty = ", penalty)) + ggtitle(paste0("Elbow plot for ", allelename)) + theme(plot.title=element_text(hjust=0.5)) 

    print(p1, vp=define_region(counter, 1:15))
    print(p2, vp=define_region(counter, 16:20))

    counter<-counter+1

  }

  dev.off()

  regtab<-do.call(rbind,penaltyregion)
  penaltyfile[[region]]<-regtab

}

alltab<-data.table(do.call(rbind,penaltyfile))
colnames(alltab)<-c("region", "allele", "applied_penalty", "alternative_penalties")

fwrite(alltab, file=file.path(opt$outputdir, "penalties.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names = FALSE)

now<-Sys.time()
message('[',now,'][Message] Done')

options(warn = defaultW)