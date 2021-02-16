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
if (!requireNamespace('ggplotify', quietly = TRUE))
  install.packages('ggplotify',repos='http://cran.us.r-project.org')
if (!requireNamespace('ggpubr', quietly = TRUE))
  install.packages('ggpubr',repos='http://cran.us.r-project.org')
if (!requireNamespace('karyoploteR', quietly = TRUE))
  BiocManager::install('karyoploteR')

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(karyoploteR))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(changepoint))
suppressPackageStartupMessages(library(tseries))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplotify))
suppressPackageStartupMessages(library(ggpubr))


cptfn <- function(data, pen) { ##https://rpubs.com/richkt/269908
  ans <- cpt.mean(data, test.stat="Normal", method = "PELT", penalty = "Manual", pen.value = pen) 
  length(cpts(ans)) +1
}


defaultW <- getOption("warn")
options(warn = -1)

option_list = list(
  make_option(c('-f', '--frequencies'), action='store', type='character', help='.tsv file containing methylation frequencies as from TREADMILL/scripts/callmethylation.sh [required]'),
  make_option(c('-o', '--outputdir'), action='store', type='character', help='output directory [required]'),
  make_option(c('-b', '--bed'), action='store', type='character', help='.bed file with one or more regions to restrict the analysis to [required]'),
  make_option(c('-r', '--release'), action='store', type='character', help='genome release [hg19]', default = 'hg19')
)

opt = parse_args(OptionParser(option_list=option_list))

now<-Sys.time()
message('[',now,'][Message] Reading input BED file')
BED<-fread(file.path(opt$bed), sep='\t', header=FALSE)

now<-Sys.time()
message('[',now,'][Message] Reading methylation frequency files')
M<-fread(file.path(opt$frequencies), sep='\t', header=TRUE)

#create output folder if it does not exist

dir.create(file.path(opt$outputdir), showWarnings = FALSE)

#process the meth frequencies based on BED lines

for (row in 1:nrow(BED)) {

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

  if (nrow(subM) == 0) {

    now<-Sys.time()
    error('[',now,'][Error] Missing region in methylation frequency file')

  }

  #split by allele

  subs<-split(subM, by="allele_name")
  subs<-setNames(lapply(sort(names(subs)), FUN = function(n) subs[[n]]), mixedsort(names(subs))) #sort so that the first 2 are main genotype and others are mosaics, if any

  #iterate over alleles
  counter<-0
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  plotlist<-list()

  for (subM1 in subs) {

    counter<-counter+1
    allelename<-paste0("allele ",counter)
    mergedM1<-list()
    splittedM1<-split(subM1, by="chrom")
    last<-as.numeric(unlist(strsplit(unique(splittedM1[[length(splittedM1)]]$chrom),"_"))[9])

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

    now<-Sys.time()
    message('[',now,'][Message] Performing changepoint detection analysis on ', allelename, ', region ', region)

    pen.vals <- seq(0, 20,.1)
    elbowplotData <- unlist(lapply(pen.vals, function(p) cptfn(data = M1res$methylated_frequency, pen = p)))
    penalty<-pen.vals[which(diff(elbowplotData) == -1)[1]]
    cptm_CP <- cpt.mean(M1res$methylated_frequency, penalty='Manual',pen.value=penalty,method='PELT', class=TRUE) 
    indexes<-c(0,cptm_CP@cpts)
    vals<-cptm_CP@param.est$mean
    svals<-rep(vals,diff(indexes))
    t1<-kpss.test(M1res$methylated_frequency, null="Trend")
    kpss1_2<-paste0("p = ", t1$p.value)
    kpss1_3<-paste0("KPSS trend = ", t1$statistic)

    now<-Sys.time()
    message('[',now,'][Message] KPSS Test for Trend Stationarity on ', allelename, ': ', kpss1_2, '. ', kpss1_3)

    M1pos<-rescale(M1res$start, c(start,end))
    M1methvals<-data.frame(chromosome=subchrom, start=M1pos, end=M1pos, y=M1res$methylated_frequency)
    M1methseg<-data.frame(chromosome=subchrom, start=M1pos, end=M1pos, y=svals)

    #convert to GRanges for visualization

    GRMV1<-makeGRangesFromDataFrame(M1methvals, keep.extra.columns = TRUE)
    GRMS1<-makeGRangesFromDataFrame(M1methseg, keep.extra.columns = TRUE)

    zoom.region <- toGRanges(data.frame(subchrom, start-50, end+50))

    now<-Sys.time()
    message('[',now,'][Message] Plotting')

    #plot

    plotlist[[counter]]<-as.ggplot(expression(kp <- plotKaryotype(genome=opt$release, plot.type = 1, zoom=zoom.region),
    kpDataBackground(kp, data.panel=1, col='grey95', r0=0.1, r1=0.9),
    kpAddBaseNumbers(kp, add.units=TRUE,minor.ticks=TRUE,minor.tick.dist=50,major.ticks=TRUE, tick.dist=500),
    kpAxis(kp, data.panel=1,r0=0.1, r1=0.9),
    kpPoints(kp, data.panel=1, data=GRMV1, cex=.4,pch=19, r0=0.1, r1=0.9),
    kpLines(kp, data.panel=1, data=GRMS1, col=col_vector[counter], lwd = 1.5,r0=0.1, r1=0.9),
    kpAddMainTitle(kp, paste0('Methylation profile of ', allelename)),
    kpAddLabels(kp, data.panel=1,labels = 'Methylation frequency', r0=0.65, r1=0.95,cex=.8,label.margin = .07, srt=90)))

  }

  #combine multiple alleles plot for the same regions

  region_to_plot<-gsub(":", "_", region)
  ptot<-ggarrange(plotlist=plotlist, labels=LETTERS[1:counter], nrow=counter)
  ggsave(file.path(opt$outputdir, paste0(region_to_plot, ".pdf")), width=15, height=7*counter)

}

now<-Sys.time()
message('[',now,'][Message] Done')

if (file.exists("Rplots.pdf")) {

  invisible(file.remove("Rplots.pdf"))

}

options(warn = defaultW)