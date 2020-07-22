#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",repos='http://cran.us.r-project.org')
if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse",repos='http://cran.us.r-project.org')
if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table",repos='http://cran.us.r-project.org')
if (!requireNamespace("annotatr", quietly = TRUE))
  BiocManager::install("annotatr")
if (!requireNamespace("karyoploteR", quietly = TRUE))
  BiocManager::install("karyoploteR")


suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(annotatr))
suppressPackageStartupMessages(library(karyoploteR))

option_list = list(
  make_option(c("-p", "--haplotype1"), action="store", type='character', help=".tsv file containing methylation frequencies for haplotype 1 [required]"),
  make_option(c("-q", "--haplotype2"), action="store", type='character', help=".tsv file containing methylation frequencies for haplotype 2 [required]"),
  make_option(c("-o", "--outputdir"), action="store", type='character', help="output directory [required]"),
  make_option(c("-r", "--release"), action="store", type='character', help="genome release [hg38]", default = "hg38"),
  make_option(c("-w", "--window"), action="store", type='numeric', help="window size for binning on methylation frequencies [50000]", default = 50000),
  make_option(c("-b", "--bed"), action="store", type='character', help=".bed file with regions to restrict the analysis to [NULL]", default = NULL)
)

opt = parse_args(OptionParser(option_list=option_list))


#hg38 release

if (opt$release == "hg38") {
  
  chromosomes<-c(paste0("chr", c(1:22)), "chrX", "chrY")
  starts<-rep(1,24)
  ends<-c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
  ref<-as.list(setNames(c(1:24), chromosomes))
  dimrange<-makeGRangesFromDataFrame(data.frame(chromosome=chromosomes,start=starts,end=ends))
  annots <-'hg38_cpg_islands'
  cpgs <-build_annotations(genome = 'hg38', annotations = annots)
} else if (opt$release == "hg19") {
  
  chromosomes<-c(paste0("chr", c(1:22)), "chrX", "chrY")
  starts<-rep(1,24)
  ends<-c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
  ref<-as.list(setNames(c(1:24), chromosomes))
  dimrange<-makeGRangesFromDataFrame(data.frame(chromosome=chromosomes,start=starts,end=ends))
  annots <-'hg19_cpg_islands'
  cpgs <-build_annotations(genome = 'hg19', annotations = annots)
} else {
  
  stop("Cannot process ", opt$release, ". Supported releases are: hg38,hg19")
}

#read methylation frequencies

#haplotype 1
meth1<-fread(file.path(opt$haplotype1), sep="\t", header=TRUE)
if ( ! startsWith(meth1$chromosome[1], "chr")) {
  
  meth1$chromosome<-paste0("chr",meth1$chromosome)
}

m1rg<-makeGRangesFromDataFrame(meth1)

#haplotype 2
meth2<-fread(file.path(opt$haplotype2), sep="\t", header=TRUE)

if ( ! startsWith(meth2$chromosome[1], "chr")) {
  
  meth2$chromosome<-paste0("chr",meth2$chromosome)
}

m2rg<-makeGRangesFromDataFrame(meth2)

#get regions

if (is.null(opt$bed)) { #plot entire chromosomes
  
  chrs<-unique(c(meth1$chromosome, meth2$chromosome))
  
  for (chr in chrs) {
    
    intervals<-tile(dimrange[ref[[chr]],], width=opt$window)[[1]]
    
    dflist1<-list()
    dflist2<-list()
    
    for (i in 1:length(intervals)) {
      
      ov1<-data.frame(findOverlaps(intervals[i], m1rg))
      ov2<-data.frame(findOverlaps(intervals[i], m2rg))
      
      int1<-data.frame(intervals[i])
      int1mod<-data.frame(chromosome=int1$seqnames, start=int1$start, end=int1$end)
      int2mod<-int1mod
      
      #haplotype1
      
      if (length(ov1$queryHits) != 0) {
        
        meanmeth1<-mean(as.numeric(meth1[ov1$subjectHits,]$methylated_frequency))
        
      } else {
        
        meanmeth1<-NA
        
      }
      
      #haplotype2
      
      if (length(ov2$queryHits) != 0) {
        
        meanmeth2<-mean(as.numeric(meth2[ov2$subjectHits,]$methylated_frequency))
        
      } else {
        
        meanmeth2<-NA
        
      }
      
      #store in list
      
      int1mod$y<-meanmeth1
      dflist1[[i]]<-int1mod
      
      int2mod$y<-meanmeth2
      dflist2[[i]]<-int2mod
      
      
    }
    
    dfmeth1<-do.call(rbind,dflist1)
    dfmeth2<-do.call(rbind,dflist2)
    
    grmeth1<-makeGRangesFromDataFrame(dfmeth1, keep.extra.columns = TRUE)
    grmeth2<-makeGRangesFromDataFrame(dfmeth2, keep.extra.columns = TRUE)
    
    pp <- getDefaultPlotParams(plot.type=2)
    pp$data1height <- 500
    pp$data2height <- 60
    
    pdf(file.path(opt$output,paste0(chr,".pdf")), width=15, height=10)
    
    kp <- plotKaryotype(chromosomes=chr, plot.type=2, plot.params = pp,genome = opt$release, main="Haplotype-specific methylation frequencies")
    kpAddBaseNumbers(kp, tick.dist = 10000000, minor.tick.dist = 1000000, add.units = TRUE, cex=0.4, digits = 2)
     
    #lower panel
    
    kpPlotRegions(kp, data=cpgs, col="#AAAAAA", border="#AAAAAA", data.panel=2)
    kpPlotDensity(kp, data=cpgs, col="#AA88FF", window.size = opt$window, r0=0, r1=1,data.panel=2)

    #upper panels. h2 lower
    
    kpRect(kp, data=dimrange , y0=0, y1=1, col="#FFDDDD", border=NA, r0=0, r1=0.45)
    kpAxis(kp, ymin = 0, ymax = 1, r0=0, r1=0.45, numticks = 5, col="#666666", cex=0.5, data.panel = 1)
    kpPoints(kp, data=grmeth2, pch=16, cex=0.3, col="black", r0=0, r1=0.45, data.panel = 1)
    #add kpLines with segmentation
    kpText(kp, chr=data.frame(regionrg)$seqnames, labels="haplotype 2", y=1.05, x=(data.frame(dimrange[ref[[chr]],])$end-data.frame(dimrange[ref[[chr]],])$start)/2, r0=0, r1=0.45, data.panel = 1, cex=.6)
    #h1 upper
    
    kpRect(kp, data=dimrange , y0=0, y1=1, col="#ADD8E6", border=NA, r0=0.50, r1=0.95)
    kpAxis(kp, ymin = 0, ymax = 1, r0=0.50, r1=0.95, numticks = 5, col="#666666", cex=0.5, data.panel = 1)
    kpPoints(kp, data=grmeth2, pch=16, cex=0.3, col="black", r0=0.50, r1=0.95, data.panel = 1)
    #add kpLines with segmentation
    kpText(kp, chr=data.frame(regionrg)$seqnames, labels="haplotype 1", y=1.05, x=(data.frame(dimrange[ref[[chr]],])$end-data.frame(dimrange[ref[[chr]],])$start)/2, r0=0.50, r1=0.95, data.panel = 1, cex=.6)
    
    dev.off()
    
  }
} else {
  
  bedregion<-fread(file.path(opt$bed), sep="\t", header=FALSE)
  if ( ! startsWith(bedregion$V1[1], "chr")) {
    
    bedregion$V1<-paste0("chr",bedregion$V1)
  }
  
  for (i in 1:nrow(bedregion)) {
    
    region<-bedregion[i,]
    regionrg<-makeGRangesFromDataFrame(data.frame(chromosome=region$V1,start=region$V2,end=region$V3))
    
    intervals<-tile(regionrg, width=opt$window)[[1]]
    
    dflist1<-list()
    dflist2<-list()
    
    for (i in 1:length(intervals)) {
      
      ov1<-data.frame(findOverlaps(intervals[i], m1rg))
      ov2<-data.frame(findOverlaps(intervals[i], m2rg))
      
      int1<-data.frame(intervals[i])
      int1mod<-data.frame(chromosome=int1$seqnames, start=int1$start, end=int1$end)
      int2mod<-int1mod
      
      #haplotype1
      
      if (length(ov1$queryHits) != 0) {
        
        meanmeth1<-mean(as.numeric(meth1[ov1$subjectHits,]$methylated_frequency))
        
      } else {
        
        meanmeth1<-NA
        
      }
      
      #haplotype2
      
      if (length(ov2$queryHits) != 0) {
        
        meanmeth2<-mean(as.numeric(meth2[ov2$subjectHits,]$methylated_frequency))
        
      } else {
        
        meanmeth2<-NA
        
      }
      
      #store in list
      
      int1mod$y<-meanmeth1
      dflist1[[i]]<-int1mod
      
      int2mod$y<-meanmeth2
      dflist2[[i]]<-int2mod
      
      
    }
    
    dfmeth1<-do.call(rbind,dflist1)
    dfmeth2<-do.call(rbind,dflist2)
    
    grmeth1<-makeGRangesFromDataFrame(dfmeth1, keep.extra.columns = TRUE)
    grmeth2<-makeGRangesFromDataFrame(dfmeth2, keep.extra.columns = TRUE)
    
    pp <- getDefaultPlotParams(plot.type=2)
    pp$data1height <- 500
    pp$data2height <- 60
    
    pdf(file.path(opt$output,paste0(region$V1,":",region$V2,"-",region$V3,".pdf")), width=15, height=10)
    
    kp <- plotKaryotype(zoom=regionrg, plot.type=2, plot.params = pp,genome = opt$release, main="Haplotype-specific methylation frequencies")
    kpAddBaseNumbers(kp, tick.dist = 10000000, minor.tick.dist = 1000000, add.units = TRUE, cex=0.4, digits = 2)
    #lower panel
    
    kpPlotRegions(kp, data=cpgs, col="#AAAAAA", border="#AAAAAA", data.panel=2, r0=0.0, r1=1)
    kpPlotDensity(kp, data=cpgs, col="#AA88FF", window.size = opt$window, r0=0.0, r1=1,data.panel=2)
    #upper panels. h2 lower
    
    kpRect(kp, data=dimrange , y0=0, y1=1, col="#FFDDDD", border=NA, r0=0, r1=0.45)
    kpAxis(kp, ymin = 0, ymax = 1, r0=0, r1=0.45, numticks = 5, col="#666666", cex=0.5, data.panel = 1)
    kpPoints(kp, data=grmeth2, pch=16, cex=0.3, col="black", r0=0, r1=0.45, data.panel = 1)
    #add kpLines with segmentation
    kpText(kp, chr=data.frame(regionrg)$seqnames, labels="haplotype 2", y=1.05, x=(data.frame(regionrg)$start+data.frame(regionrg)$end)/2, r0=0, r1=0.45, data.panel = 1, cex=.6)
    #h1 upper
    
    kpRect(kp, data=dimrange , y0=0, y1=1, col="#ADD8E6", border=NA, r0=0.50, r1=0.95)
    kpAxis(kp, ymin = 0, ymax = 1, r0=0.50, r1=0.95, numticks = 5, col="#666666", cex=0.5, data.panel = 1)
    kpPoints(kp, data=grmeth2, pch=16, cex=0.3, col="black", r0=0.50, r1=0.95, data.panel = 1)
    #add kpLines with segmentation
    kpText(kp, chr=data.frame(regionrg)$seqnames, labels="haplotype 1", y=1.05, x=(data.frame(regionrg)$start+data.frame(regionrg)$end)/2, r0=0.50, r1=0.95, data.panel = 1, cex=.6)
    
    dev.off()    
  }
}
