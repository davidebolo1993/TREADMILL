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
  make_option(c('-j', '--json'), action='store', type='character', help='on ore more comma-separated .json file from TREADMILL BASIC [required]'),
  make_option(c('-l', '--labels'), action='store', type='character', help='one ore more comma-separated label to identify each experiment [required]'),
  make_option(c('-o', '--output'), action='store', type='character', help='output .pdf file with calculated stats [out.pdf]', default='out.pdf')
)

opt = parse_args(OptionParser(option_list=option_list))


listall<-list()
listerr<-list()
listcov<-list()
liststats<-list()

if (is.null(opt$json)) {

  now<-Sys.time()
  stop('[',now,'][Error] At least one .json from TREADMILL BASIC is required')

}

jsons<-unlist(strsplit(opt$json, ','))
labels<-unlist(strsplit(opt$labels, ','))

if (length(jsons) != length(labels)) {

  now<-Sys.time()
  stop('[',now,'][Error] Number of input .json files does not match number of labels provided')

}

for (i in 1:length(jsons)) {

  now<-Sys.time()
  injson<-jsons[i]
  inlabel<-labels[i]
  message('[',now,'][Message] Processing .json file ', file.path(injson), ' with label ', inlabel)

  if (file.exists(file.path(injson))) {

    mjson<-fromJSON(file = file.path(injson))
    mkeyall<-c('unmapped', 'supplementary', 'secondary', 'primary', 'on-target', 'off-target')
    mvalueall<-c(mjson$BAM_UNMAP,mjson$BAM_SUPP,mjson$BAM_SEC,mjson$BAM_PRIM,mjson$BAM_ONTARGET,mjson$BAM_OFFTARGET)
    mlabelall<-rep(inlabel,length(mkeyall))

    mall<-data.frame(x=mkeyall,y=mvalueall,z=mlabelall,stringsAsFactors = FALSE)
    listall[[inlabel]]<-mall

    mkeyerr<-c('match', 'mismatch', 'deletion', 'insertion', 'soft-clipped')
    mvalueerr<-c(mjson$BAM_CMATCH, mjson$BAM_CDIFF, mjson$BAM_CDEL, mjson$BAM_CINS, mjson$BAM_CSOFT_CLIP)
    mlabelerr<-rep(inlabel,length(mkeyerr))

    merror<-data.frame(x=mlabelerr,y=mvalueerr,z=mkeyerr,stringsAsFactors = FALSE)
    listerr[[inlabel]]<-merror

    mstats<-data.frame(x=mjson$BAM_LEN, y=mjson$BAM_QUAL, w=mjson$BAM_PID, z=inlabel,stringsAsFactors = FALSE)
    liststats[[inlabel]]<-mstats
    
    mcovkeys<-grep(':', names(mjson), value=TRUE)
    mtmplist<-list()
    
    for (mk in mcovkeys) {
      
      y<-as.numeric(as.character(unlist(mjson[mk])))
      x<-rep(mk, length(y))
      z<-rep(inlabel,length(y))       
      mtmplist[[mk]]<-data.frame(x=x,y=y,z=z,stringsAsFactors = FALSE)
    }

    listcov[[inlabel]]<-do.call(rbind,mtmplist)
  
  } else {

    now<-Sys.time()
    stop('[',now,'][Error] .json file ', file.path(injson), ' does not exist')

  }

}



now<-Sys.time()
message('[',now,'][Message] Plotting')

#overall

dfall<-do.call(rbind,listall)
dfall$x<-factor(dfall$x,levels=unique(dfall$x))
dfall$z<-factor(dfall$z,levels=labels)
dfall_2<-subset(dfall,x=='on-target' | x =='off-target')
dfall_mod1 <- ddply(dfall_2, .(z), transform, percent = y/sum(y) * 100)
dfall_3<-subset(dfall,x!='on-target' & x !='off-target')
dfall_3$percent<-0.0
dfall_mod<-rbind(dfall_3,dfall_mod1)
dfall_mod<- ddply(dfall_mod, .(z), transform, pos = (cumsum(y) - 0.5 * y))
dfall_mod$label <- paste0(sprintf('%.3f', dfall_mod$percent), '%')
dfall_mod$zoom <- TRUE

pall<-ggplot(dfall, aes(x=as.numeric(x), y=y, fill=z))+
  geom_bar(stat='identity',width = .3,position="dodge")+
  theme_bw()+
  scale_fill_brewer(palette='Dark2')+
  ylab('# reads') + 
  xlab('read type')+
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  ggtitle('Overall statistics')+
  facet_zoom(xlim=c(4.5,6.5),ylim=c(0,dfall$y[which(dfall$x=='on-target')]+100), horizontal=FALSE, zoom.data=zoom)+
  geom_text(data=dfall_mod, aes(label = label), position=position_dodge(width=.3), size=2,vjust = -0.5)+
  scale_x_continuous(
    breaks = 1:length(unique(dfall$x)),
    label = levels(dfall$x)
  )

#on-target

dferr<-do.call(rbind,listerr)
dferr$x<-factor(dferr$x,levels=unique(dferr$x))
#dferr$z<-factor(dferr$z,levels=labels)

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
dfstats$z<-factor(dfstats$z,levels=labels)
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
  theme(legend.position='bottom', legend.background=element_blank(),legend.direction="horizontal", legend.title=element_text(face="bold.italic"), strip.background =element_rect(fill="grey20"), strip.text = element_text(colour = 'white'), plot.title = element_text(hjust = 0.5))+
  facet_wrap(~z, scales = "fixed")+
  ggtitle('% identity vs quality, on-target reads')

p_length<-ggplot(data=dfstats, aes(x=x, fill=z)) + 
  geom_histogram(position="dodge", bins=50)+
  geom_vline(data=lenmean, aes(xintercept = xmean, color=z),linetype='dashed', show.legend=FALSE)+
  scale_fill_brewer(palette='Dark2') + 
  theme_bw()+
  scale_x_continuous("read length (#bp)") + 
  theme(legend.position='bottom', legend.background=element_blank(),legend.direction="horizontal", legend.title=element_blank(),plot.title = element_text(hjust = 0.5)) +
  ggtitle('Read length distribution, on-target reads')

#p_len_vs_pid<-ggplot(data=dfstats, aes(x=w, y=x)) + 
#  geom_point(size = 0.6)+ stat_density2d(aes(col=..level.., alpha=..level..))+
#  geom_rug(sides="t", size=0.05, col=rgb(.8,0,0,alpha=.3)) + 
#  geom_rug(sides="r", size=0.05, col=rgb(0,0,.8,alpha=.3)) + 
#  scale_x_continuous("% identity") + 
#  scale_y_continuous("length (#bp)") + 
#  theme_bw()+
#  geom_vline(data=idmean, aes(xintercept = xmean),linetype='dashed', col='darkgreen')+
#  geom_hline(data=lenmean,aes(yintercept = xmean), linetype='dashed', col='darkgreen')+  
#  scale_color_continuous(low="darkblue",high="darkred")+
#  guides(alpha="none",col=guide_legend(title="Density"))+
#  theme(legend.position='bottom', legend.background=element_blank(),legend.direction="horizontal", legend.title=element_text(face="bold.italic"),strip.background =element_rect(fill="grey20"), strip.text = element_text(colour = 'white'), plot.title = element_text(hjust = 0.5))+
#  facet_wrap(~z, scales = "free")+
#  ggtitle('% identity vs length, on-target reads')

#p_qual_vs_len<-ggplot(data=dfstats, aes(x=x, y=y)) + 
#  geom_point(size = 0.6)+ stat_density2d(aes(col=..level.., alpha=..level..))+
#  geom_rug(sides="t", size=0.05, col=rgb(.8,0,0,alpha=.3)) + 
#  geom_rug(sides="r", size=0.05, col=rgb(0,0,.8,alpha=.3)) + 
#  scale_x_continuous("length (#bp)") + 
#  scale_y_continuous("quality (phred)") + 
#  theme_bw()+
#  geom_vline(data=lenmean, aes(xintercept = xmean),linetype='dashed', col='darkgreen')+
#  geom_hline(data=qualmean,aes(yintercept = xmean), linetype='dashed', col='darkgreen')+
#  scale_color_continuous(low="darkblue",high="darkred")+
#  guides(alpha="none",col=guide_legend(title="Density"))+
#  theme(legend.position='bottom', legend.background=element_blank(),legend.direction="horizontal", legend.title=element_text(face="bold.italic"),strip.background =element_rect(fill="grey20"), strip.text = element_text(colour = 'white'), plot.title = element_text(hjust = 0.5))+
#  facet_wrap(~z, scales = "free")+
#  ggtitle('Length vs quality, on-target reads')
  


#coverage per region

dfcov<-do.call(rbind,listcov)
dfcov$x<-factor(dfcov$x,levels=unique(dfcov$x))
dfcov$z<-factor(dfcov$z,levels=labels)
pcov<-ggplot(dfcov, aes(x=x, y=y, fill=z)) +
  geom_boxplot(position=position_dodge(0.3), width=0.2/length(unique(dfcov$x)))+
  scale_fill_brewer(palette='Dark2') + theme_bw()+
  xlab('regions') + 
  ylab('coverage') +
  theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5)) +
  ggtitle('Coverage on targets')


pdf(file.path(opt$output), height=20, width=13, onefile=TRUE)
grid.newpage()
pushViewport(viewport(layout = grid.layout(8, 10)))
print(pall, vp = vplayout(1:2, 1:5))
print(perr, vp = vplayout(1:2, 6:10))
print(p_qual_vs_pid, vp = vplayout(3:4, 2:9))
print(p_length, vp = vplayout(5:6, 2:9))
#print(p_len_vs_pid, vp = vplayout(7:8, 2:9))
print(pcov, vp = vplayout(7:8, 1:10))
dev.off()

options(warn = defaultW)
