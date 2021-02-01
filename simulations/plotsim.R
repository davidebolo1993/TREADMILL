#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(ggplot2)
library(tidyr)

rect_df <- data.frame(xmin = c(-Inf,-Inf,-Inf,-Inf), xmax = c(Inf,Inf,Inf,Inf), ymin = c(0,45,55,200), ymax = c(44,54,199,Inf),grp = factor(c('Normal','Greyzone','Premutated','Affected'),levels = c('Normal','Greyzone','Premutated','Affected')))

#overallplot

res<-fread(file.path(args[1]), sep='\t', header=TRUE)
res2<-separate(data=res, col=predicted, into=c("left", "right"))
res2.a<-data.frame('allele' = res2$left, 'coverage' = res2$coverage, 'identity' = res2$identity, 'algorithm' = res2$algorithm)
res2.b<-data.frame('allele' = res2$right, 'coverage' = res2$coverage, 'identity' = res2$identity, 'algorithm' = res2$algorithm)
res2.tot<-rbind(res2.a,res2.b)
res2.tot$coverage<-factor(res2.tot$coverage, levels=c('20X', '50X', '100X', '200X'))
res2.tot$identity<-factor(res2.tot$identity, levels=c('85%', '90%', '95%'))
res2.tot$algorithm<-factor(res2.tot$algorithm, levels=c('DBSCAN', 'AGGLOMERATIVE'))

overview<-ggplot(res2.tot) + geom_jitter(aes(x=coverage, y=as.numeric(allele), group=coverage),size=1, na.rm = TRUE) + scale_y_continuous(limits=c(0,300), breaks = c(0, 20, 50, 100, 150, 200, 250, 300)) + theme_classic() +
  geom_hline(yintercept = 20, linetype='dashed', size=.2) + geom_hline(yintercept= 50,linetype='dashed', size=.2) + geom_hline(yintercept =150,linetype='dashed', size=.2) + geom_hline(yintercept= 250,linetype='dashed', size=.2) +
  labs(x='coverage', y='predicted #CGG repeats')+ geom_rect(data = rect_df,aes(x = NULL,y = NULL,xmin = xmin,xmax = xmax,ymin = ymin,ymax = ymax,fill = grp),alpha = 0.3)+
  scale_fill_manual(values = c('darkgreen','grey50','darkblue', 'darkred')) + theme(legend.title=element_blank(), legend.position = 'bottom', legend.direction = 'horizontal', strip.background =element_rect(fill='grey20'), strip.text = element_text(colour = 'white')) + facet_grid(identity ~ algorithm)

ggsave(file.path(args[1],'overview.pdf'), overview, width=20, height=15)

#truepositiverate

tres<-fread(file.path(args[2]), sep='\t', header=TRUE)

tres$errors<-factor(tres$errors, levels=c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10'))
tres$coverage<-factor(tres$coverage, levels=c('20X', '50X', '100X', '200X'))
tres$identity<-factor(tres$identity, levels=c('85%', '90%', '95%'))
tres$algorithm<-factor(tres$algorithm, levels=c('DBSCAN', 'AGGLOMERATIVE'))
tres$status<-factor(tres$status, levels=c('greyzone', 'premutated', 'affected'))

colors<-c('greyzone' = 'grey50','premutated'='darkblue', 'affected' = 'darkred')
shapes<-c('DBSCAN' = 8, 'AGGLOMERATIVE' = 1)
tpr<-ggplot(tres) + geom_line(aes(x=errors, y=as.numeric(TPR), group=interaction(status,algorithm), col=status)) +geom_point(aes(x=errors, y=as.numeric(TPR), group=interaction(status,algorithm), shape=algorithm))  + scale_color_manual(values=colors) + theme_bw() + theme(legend.title=element_blank(), legend.position = 'bottom', legend.direction = 'horizontal', strip.background =element_rect(fill="grey20"), strip.text = element_text(colour = 'white')) + labs(y='TPR', x='# errors in CGG prediction')+
  scale_shape_manual(values=shapes) + facet_grid(identity~coverage) + scale_y_continuous(breaks = seq(0.00,1.00,by=0.10))

ggsave(file.path(args[1],'tpr.pdf'), tpr, width=20, height=15)


