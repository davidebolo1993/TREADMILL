#!/usr/bin/env Rscript

if (!requireNamespace('optparse', quietly = TRUE))
  install.packages('optparse',repos='http://cran.us.r-project.org')
if (!requireNamespace('rjson', quietly = TRUE))
  install.packages('rjson',repos='http://cran.us.r-project.org')
if (!requireNamespace('ggplot2', quietly = TRUE))
  install.packages('ggplot2',repos='http://cran.us.r-project.org')
if (!requireNamespace('plyr', quietly = TRUE))
  install.packages('plyr',repos='http://cran.us.r-project.org')
if (!requireNamespace('ggrepel', quietly = TRUE))
  install.packages('ggrepel',repos='http://cran.us.r-project.org')
if (!requireNamespace('ggforce', quietly = TRUE))
  install.packages('ggforce',repos='http://cran.us.r-project.org')


suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(grid))

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

defaultW <- getOption("warn")
options(warn = -1)

option_list = list(
  make_option(c('-m', '--minimap2'), action='store', type='character', help='.json file with informations from a minimap2-generated alignment'),
  make_option(c('-n', '--ngmlr'), action='store', type='character', help='.json file with informations from a ngmlr-generated alignment'),
  #make_option(c('-l', '--last'), action='store', type='character', help='.json file with informations from a last-generated alignment'),
  make_option(c('-o', '--output'), action='store', type='character', help='output .pdf file with calculated stats [out.pdf]', default='out.pdf')
)

opt = parse_args(OptionParser(option_list=option_list))


listall<-list()
listerr<-list()
listcov<-list()
liststats<-list()

if (! is.null(opt$minimap2)) {

  now<-Sys.time()
  message('[',now,'][Message] Processing .json file from minimap2')

  if (file.exists(file.path(opt$minimap2))) {
    mjson<-fromJSON(file = file.path(opt$minimap2))
    mkeyall<-c('unmapped', 'supplementary', 'secondary', 'primary', 'on-target', 'off-target')
    mvalueall<-c(mjson$BAM_UNMAP,mjson$BAM_SUPP,mjson$BAM_SEC,mjson$BAM_PRIM,mjson$BAM_ONTARGET,mjson$BAM_OFFTARGET)
    malignerall<-rep('minimap2',length(mkeyall))
    
    mall<-data.frame(x=mkeyall,y=mvalueall,z=malignerall,stringsAsFactors = FALSE)
    listall[['m']]<-mall
      
    mkeyerr<-c('match', 'mismatch', 'deletion', 'insertion', 'soft-clipped')
    mvalueerr<-c(mjson$BAM_CMATCH, mjson$BAM_CDIFF, mjson$BAM_CDEL, mjson$BAM_CINS, mjson$BAM_CSOFT_CLIP)
    malignererr<-rep('minimap2',length(mkeyerr))
    
    merror<-data.frame(x=malignererr,y=mvalueerr,z=mkeyerr,stringsAsFactors = FALSE)
    listerr[['m']]<-merror

    mstats<-data.frame(x=mjson$BAM_LEN, y=mjson$BAM_QUAL, w=mjson$BAM_PID, z=c('minimap2'),stringsAsFactors = FALSE)
    liststats[['m']]<-mstats
    
    mcovkeys<-grep(':', names(mjson), value=TRUE)
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
  

  now<-Sys.time()
  message('[',now,'][Message] Processing .json file from ngmlr')
  
  if (file.exists(file.path(opt$ngmlr))) {
    njson<-fromJSON(file = file.path(opt$ngmlr))
    nkeyall<-c('unmapped', 'supplementary', 'secondary', 'primary', 'on-target', 'off-target')
    nvalueall<-c(njson$BAM_UNMAP,njson$BAM_SUPP,njson$BAM_SEC,njson$BAM_PRIM,njson$BAM_ONTARGET,njson$BAM_OFFTARGET)
    nalignerall<-rep('ngmlr',length(nkeyall))
    
    nall<-data.frame(x=nkeyall,y=nvalueall,z=nalignerall,stringsAsFactors = FALSE)
    listall[['n']]<-nall
    
    nkeyerr<-c('match', 'mismatch', 'deletion', 'insertion', 'soft-clipped')
    nvalueerr<-c(njson$BAM_CMATCH, njson$BAM_CDIFF, njson$BAM_CDEL, njson$BAM_CINS, njson$BAM_CSOFT_CLIP)
    nalignererr<-rep('ngmlr',length(nkeyerr))
    
    nerror<-data.frame(x=nalignererr,y=nvalueerr,z=nkeyerr,stringsAsFactors = FALSE)
    listerr[['n']]<-nerror

    nstats<-data.frame(x=njson$BAM_LEN, y=njson$BAM_QUAL, w=njson$BAM_PID, z=c('ngmlr'),stringsAsFactors = FALSE)
    liststats[['n']]<-nstats
   
    ncovkeys<-grep(':', names(njson), value=TRUE)
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

#if (! is.null(opt$last)) {

  #now<-Sys.time()
  #message('[',now,'][Message] Processing .json file from last')

  #if (file.exists(file.path(opt$last))) {
    #ljson<-fromJSON(file = file.path(opt$last))
    #lkeyall<-c('unmapped', 'supplementary', 'secondary', 'primary', 'on-target', 'off-target')
    #lvalueall<-c(ljson$BAM_UNMAP,ljson$BAM_SUPP,ljson$BAM_SEC,ljson$BAM_PRIM,ljson$BAM_ONTARGET,ljson$BAM_OFFTARGET)
    #lalignerall<-rep('last',length(lkeyall))

    #lall<-data.frame(x=lkeyall,y=lvalueall,z=lalignerall,stringsAsFactors = FALSE)
    #listall[['l']]<-lall

    #lkeyerr<-c('match', 'mismatch', 'deletion', 'insertion', 'soft-clipped')
    #lvalueerr<-c(ljson$BAM_CMATCH, ljson$BAM_CDIFF, ljson$BAM_CDEL, ljson$BAM_CINS, ljson$BAM_CSOFT_CLIP)
    #lalignererr<-rep('ngmlr',length(lkeyerr))

    #lerror<-data.frame(x=lalignererr,y=lvalueerr,z=lkeyerr,stringsAsFactors = FALSE)
    #listerr[['l']]<-lerror

    #lstats<-data.frame(x=ljson$BAM_LEN, y=ljson$BAM_QUAL, w=ljson$BAM_PID, z=c('last'),stringsAsFactors = FALSE)
    #liststats[['l']]<-lstats

    #lcovkeys<-grep(':', names(ljson), value=TRUE)
    #ltmplist<-list()

    #for (lk in lcovkeys) {

      #y<-as.numeric(as.character(unlist(ljson[lk])))
      #x<-rep(lk, length(y))
      #z<-rep('last',length(y))
      #ltmplist[[lk]]<-data.frame(x=x,y=y,z=z,stringsAsFactors = FALSE)
    #}

    #listcov[['l']]<-do.call(rbind,ltmplist)

  #}
#}

if (length(listall) == 0) {
  now<-Sys.time()
  stop('[',now,'][Error] At least one .json from TREADMILL BASIC is required')
}


now<-Sys.time()
message('[',now,'][Message] Plotting')

#overall

dfall<-do.call(rbind,listall)
dfall$x<-factor(dfall$x,levels=unique(dfall$x))

pall<-ggplot(dfall, aes(x=as.numeric(x), y=y, fill=z))+
  geom_bar(stat='identity',width = .3,position=position_dodge())+
  theme_bw()+
  scale_fill_brewer(palette='Dark2')+
  ylab('# reads') + 
  xlab('read type')+
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  ggtitle('Overall statistics')+
  facet_zoom(xlim=c(4,6),ylim=c(0,dfall$y[which(dfall$x=='on-target')]+100), horizontal=FALSE )+
  scale_x_continuous(
    breaks = 1:length(unique(dfall$x)),
    label = levels(dfall$x)
  )

#on-target

dferr<-do.call(rbind,listerr)
dferr$x<-factor(dferr$x,levels=unique(dferr$x))

dferr <- ddply(dferr, .(x), transform, percent = y/sum(y) * 100)
dferr <- ddply(dferr, .(x), transform, pos = (cumsum(y) - 0.5 * y))
dferr$label = paste0(sprintf('%.0f', dferr$percent), '%')

perr<-ggplot(dferr, aes(x=x, y=y, fill=z)) +
  geom_bar(stat='identity',width = .1)+
  theme_bw()+
  scale_fill_brewer(palette='Dark2')+
  ylab('# bp') +
  xlab('aligner') +
  theme(legend.title=element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2)+
  ggtitle('Base composition of on-target reads')


dfstats<-do.call(rbind, liststats)

#pre-compute means per group

qualmean<- ddply(dfstats, .(z), summarise, xmean=mean(y))
lenmean<- ddply(dfstats, .(z), summarise, xmean=mean(x))
idmean<- ddply(dfstats, .(z), summarise, xmean=mean(w))



p_qual_vs_pid<-ggplot(data=dfstats, aes(x=w, y=y)) + 
  geom_point(size = 0.6)+ stat_density2d(aes(col=..level.., alpha=..level..))+
  geom_rug(sides="t", size=0.05, col=rgb(.8,0,0,alpha=.3)) + 
  geom_rug(sides="r", size=0.05, col=rgb(0,0,.8,alpha=.3)) + 
  scale_x_continuous("% identity") + 
  scale_y_continuous("quality (phred)") + 
  theme_bw()+
  geom_vline(data=idmean, aes(xintercept = xmean),linetype='dashed', col='darkgreen')+
  geom_hline(data=qualmean,aes(yintercept = xmean), linetype='dashed', col='darkgreen')+
  scale_color_continuous(low="darkblue",high="darkred")+
  guides(alpha="none",col=guide_legend(title="Density"))+
  theme(legend.position='bottom', legend.background=element_blank(),legend.direction="horizontal", legend.title=element_text(face="bold.italic"),strip.background =element_rect(fill="white"),plot.title = element_text(hjust = 0.5))+
  facet_wrap(~z, scales = "free")+
  ggtitle('% identity vs quality, on-target reads')


p_qual_vs_len<-ggplot(data=dfstats, aes(x=x, y=y)) + 
  geom_point(size = 0.6)+ stat_density2d(aes(col=..level.., alpha=..level..))+
  geom_rug(sides="t", size=0.05, col=rgb(.8,0,0,alpha=.3)) + 
  geom_rug(sides="r", size=0.05, col=rgb(0,0,.8,alpha=.3)) + 
  scale_x_continuous("length (#bp)") + 
  scale_y_continuous("quality (phred)") + 
  theme_bw()+
  geom_vline(data=lenmean, aes(xintercept = xmean),linetype='dashed', col='darkgreen')+
  geom_hline(data=qualmean,aes(yintercept = xmean), linetype='dashed', col='darkgreen')+
  scale_color_continuous(low="darkblue",high="darkred")+
  guides(alpha="none",col=guide_legend(title="Density"))+
  theme(legend.position='bottom', legend.background=element_blank(),legend.direction="horizontal", legend.title=element_text(face="bold.italic"),strip.background =element_rect(fill="white"),plot.title = element_text(hjust = 0.5))+
  facet_wrap(~z, scales = "free")+
  ggtitle('Length vs quality, on-target reads')
  

p_len_vs_pid<-ggplot(data=dfstats, aes(x=w, y=x)) + 
  geom_point(size = 0.6)+ stat_density2d(aes(col=..level.., alpha=..level..))+
  geom_rug(sides="t", size=0.05, col=rgb(.8,0,0,alpha=.3)) + 
  geom_rug(sides="r", size=0.05, col=rgb(0,0,.8,alpha=.3)) + 
  scale_x_continuous("% identity") + 
  scale_y_continuous("length (#bp)") + 
  theme_bw()+
  geom_vline(data=idmean, aes(xintercept = xmean),linetype='dashed', col='darkgreen')+
  geom_hline(data=lenmean,aes(yintercept = xmean), linetype='dashed', col='darkgreen')+  
  scale_color_continuous(low="darkblue",high="darkred")+
  guides(alpha="none",col=guide_legend(title="Density"))+
  theme(legend.position='bottom', legend.background=element_blank(),legend.direction="horizontal", legend.title=element_text(face="bold.italic"),strip.background =element_rect(fill="white"),plot.title = element_text(hjust = 0.5))+
  facet_wrap(~z, scales = "free")+
  ggtitle('% identity vs length, on-target reads')

#coverage per region

dfcov<-do.call(rbind,listcov)
dfcov$x<-factor(dfcov$x,levels=unique(dfcov$x))

pcov<-ggplot(dfcov, aes(x=x, y=y, fill=z)) +
  geom_boxplot(position=position_dodge(0.3), width=0.2/length(unique(dfcov$x)))+
  scale_fill_brewer(palette='Dark2') + theme_bw()+
  xlab('regions') + 
  ylab('coverage') +
  theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5)) +
  ggtitle('Coverage on targets')


pdf(file.path(opt$output), height=20, width=13, onefile=TRUE)
grid.newpage()
pushViewport(viewport(layout = grid.layout(10, 2)))
print(pall, vp = vplayout(1:2, 1))
print(perr, vp = vplayout(1:2, 2))
print(p_qual_vs_len, vp = vplayout(3:4, 1:2))
print(p_qual_vs_pid, vp = vplayout(5:6, 1:2))
print(p_len_vs_pid, vp = vplayout(7:8, 1:2))
print(pcov, vp = vplayout(9:10, 1:2))
dev.off()

options(warn = defaultW)
