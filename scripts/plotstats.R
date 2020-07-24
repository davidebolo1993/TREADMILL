#!/usr/bin/env Rscript

if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse",repos='http://cran.us.r-project.org')
if (!requireNamespace("rjson", quietly = TRUE))
  install.packages("rjson",repos='http://cran.us.r-project.org')
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2",repos='http://cran.us.r-project.org')
if (!requireNamespace("plyr", quietly = TRUE))
  install.packages("plyr",repos='http://cran.us.r-project.org')
if (!requireNamespace("ggrepel", quietly = TRUE))
  install.packages("ggrepel",repos='http://cran.us.r-project.org')
if (!requireNamespace("ggforce", quietly = TRUE))
  install.packages("ggforce",repos='http://cran.us.r-project.org')


suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(grid))


vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

option_list = list(
  make_option(c("-m", "--minimap2"), action="store", type='character', help=".json file with informations from a minimap2-generated alignment"),
  make_option(c("-n", "--ngmlr"), action="store", type='character', help=".json file with informations from a ngmlr-generated alignment"),
  make_option(c("-l", "--last"), action="store", type='character', help=".json file with informations from a last-generated alignment"),
  make_option(c("-o", "--output"), action="store", type='character', help="output .pdf file with calculated stats")
)

opt = parse_args(OptionParser(option_list=option_list))


listall<-list()
listerr<-list()
listcov<-list()

if (! is.null(opt$minimap2)) {
  
  if (file.exists(file.path(opt$minimap2))) {
    mjson<-fromJSON(file = file.path(opt$minimap2))
    mkeyall<-c("unmapped", "supplementary", "secondary", "primary", "on-target", "off-target")
    mvalueall<-c(mjson$BAM_UNMAP,mjson$BAM_SUPP,mjson$BAM_SEC,mjson$BAM_PRIM,mjson$BAM_ONTARGET,mjson$BAM_OFFTARGET)
    malignerall<-rep('minimap2',length(mkeyall))
    
    mall<-data.frame(x=mkeyall,y=mvalueall,z=malignerall,stringsAsFactors = FALSE)
    listall[['m']]<-mall
    
    
    mkeyerr<-c("match", "mismatch", "deletion", "insertion", "soft-clipped")
    mvalueerr<-c(mjson$BAM_CMATCH, mjson$BAM_CDIFF, mjson$BAM_CDEL, mjson$BAM_CINS, mjson$BAM_CSOFT_CLIP)
    malignererr<-rep('minimap2',length(mkeyerr))
    
    merror<-data.frame(x=malignererr,y=mvalueerr,z=mkeyerr,stringsAsFactors = FALSE)
    listerr[['m']]<-merror
    
    mcovkeys<-grep(":", names(mjson), value=TRUE)
    mtmplist<-list()
    
    for (mk in mcovkeys) {
      
      y<-as.numeric(as.character(unlist(mjson[mk])))
      x<-rep(mk, length(y))
      z<-rep('minimap2',length(y))       
      mtmplist[[mk]]<-data.frame(x=x,y=y,z=z,stringsAsFactors = FALSE)
    }
    
    listcov[['m']]<-do.call(rbind,mtmplist)
  }
}

if (! is.null(opt$ngmlr)) {
  
  
  if (file.exists(file.path(opt$ngmlr))) {
    njson<-fromJSON(file = file.path(opt$ngmlr))
    nkeyall<-c("unmapped", "supplementary", "secondary", "primary", "on-target", "off-target")
    nvalueall<-c(njson$BAM_UNMAP,njson$BAM_SUPP,njson$BAM_SEC,njson$BAM_PRIM,njson$BAM_ONTARGET,njson$BAM_OFFTARGET)
    nalignerall<-rep('ngmlr',length(nkeyall))
    
    nall<-data.frame(x=nkeyall,y=nvalueall,z=nalignerall,stringsAsFactors = FALSE)
    listall[['n']]<-nall
    
    nkeyerr<-c("match", "mismatch", "deletion", "insertion", "soft-clipped")
    nvalueerr<-c(njson$BAM_CMATCH, njson$BAM_CDIFF, njson$BAM_CDEL, njson$BAM_CINS, njson$BAM_CSOFT_CLIP)
    nalignererr<-rep('ngmlr',length(nkeyerr))
    
    nerror<-data.frame(x=nalignererr,y=nvalueerr,z=nkeyerr,stringsAsFactors = FALSE)
    listerr[['n']]<-nerror
    
    ncovkeys<-grep(":", names(njson), value=TRUE)
    ntmplist<-list()
    
    for (nk in ncovkeys) {
      
      y<-as.numeric(as.character(unlist(njson[nk])))
      x<-rep(nk, length(y))
      z<-rep('ngmlr',length(y))       
      ntmplist[[nk]]<-data.frame(x=x,y=y,z=z,stringsAsFactors = FALSE)
    }
    
    listcov[['n']]<-do.call(rbind,ntmplist)
    
  }
}

if (! is.null(opt$last)) {


  if (file.exists(file.path(opt$last))) {
    ljson<-fromJSON(file = file.path(opt$minimap2))
    lkeyall<-c("unmapped", "supplementary", "secondary", "primary", "on-target", "off-target")
    lvalueall<-c(ljson$BAM_UNMAP,ljson$BAM_SUPP,ljson$BAM_SEC,ljson$BAM_PRIM,ljson$BAM_ONTARGET,ljson$BAM_OFFTARGET)
    lalignerall<-rep('last',length(lkeyall))

    lall<-data.frame(x=lkeyall,y=lvalueall,z=lalignerall,stringsAsFactors = FALSE)
    listall[['l']]<-lall

    lkeyerr<-c("match", "mismatch", "deletion", "insertion", "soft-clipped")
    lvalueerr<-c(ljson$BAM_CMATCH, ljson$BAM_CDIFF, ljson$BAM_CDEL, ljson$BAM_CINS, ljson$BAM_CSOFT_CLIP)
    lalignererr<-rep('ngmlr',length(lkeyerr))

    lerror<-data.frame(x=lalignererr,y=lvalueerr,z=lkeyerr,stringsAsFactors = FALSE)
    listerr[['l']]<-lerror

    lcovkeys<-grep(":", names(ljson), value=TRUE)
    ltmplist<-list()

    for (lk in lcovkeys) {

      y<-as.numeric(as.character(unlist(ljson[lk])))
      x<-rep(lk, length(y))
      z<-rep('last',length(y))
      ltmplist[[lk]]<-data.frame(x=x,y=y,z=z,stringsAsFactors = FALSE)
    }

    listcov[['l']]<-do.call(rbind,ltmplist)

  }
}

if (length(listall) == 0) {
  stop('[Error] [Usage] Rscript plotstats.R -m <minimap2.json> -n <ngmlr.json> -l <last.json>. At least one argument is required')
}


#overall

dfall<-do.call(rbind,listall)
dfall$x<-factor(dfall$x,levels=unique(dfall$x))

pall<-ggplot(dfall, aes(x=as.numeric(x), y=y, fill=z))+
  geom_bar(stat="identity",width = .3,position=position_dodge())+
  theme_bw()+
  scale_fill_brewer(palette="Dark2")+
  ylab('# reads') + 
  xlab('read type')+
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  ggtitle('Overall statistics')+
  facet_zoom(xlim=c(4,6),ylim=c(0,dfall$y[which(dfall$x=="on-target")]+100), horizontal=FALSE )+
  scale_x_continuous(
    breaks = 1:length(unique(dfall$x)),
    label = levels(dfall$x)
  )

#on-target

dferr<-do.call(rbind,listerr)
dferr$x<-factor(dferr$x,levels=unique(dferr$x))

dferr <- ddply(dferr, .(x), transform, percent = y/sum(y) * 100)
dferr <- ddply(dferr, .(x), transform, pos = (cumsum(y) - 0.5 * y))
dferr$label = paste0(sprintf("%.0f", dferr$percent), "%")

perr<-ggplot(dferr, aes(x=x, y=y, fill=z)) +
  geom_bar(stat="identity",width = .1)+
  theme_bw()+
  scale_fill_brewer(palette="Dark2")+
  ylab('# bases') +
  xlab('aligner') +
  theme(legend.title=element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2)+
  ggtitle('Base composition of on-target reads')


#coverage per region

dfcov<-do.call(rbind,listcov)
dfcov$x<-factor(dfcov$x,levels=unique(dfcov$x))

pcov<-ggplot(dfcov, aes(x=x, y=y, fill=z)) +
  geom_boxplot(position=position_dodge(), width=0.3/length(unique(dfcov$x)))+
  scale_fill_brewer(palette="Dark2") + theme_bw()+
  xlab('regions') + 
  ylab('coverage') +
  theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5)) +
  ggtitle('Coverage on targets')

pdf(file.path(opt$output), height=10, width=15, onefile=TRUE)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(pall, vp = vplayout(1, 1))
print(perr, vp = vplayout(1, 2))
print(pcov, vp = vplayout(2, 1:2))
dev.off()
