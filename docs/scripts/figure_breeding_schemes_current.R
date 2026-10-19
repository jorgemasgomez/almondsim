# Reproduce the complete, current scheme figure entirely with base R graphics.
# No SVG template, Python, JavaScript or external packages required.
# Run: Rscript figure_breeding_schemes_current.R [output_directory]
a<-commandArgs(TRUE)
out<-if(length(a))a[1] else 'outputs/ocs_figure_R_20261019'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
row<-function(label,parent,early=FALSE,renew=8,ocs=FALSE)list(label=label,parent=parent,early=early,renew=renew,ocs=ocs)
blocks<-list(
 list(title='Selection accuracy',color='#2677AD',rows=list(row('Pheno\u20139',9),row('Pedigree\u20139',9),row('GS\u20139',9))),
 list(title='Generation interval',color='#5265A1',rows=list(row('Pheno / Pedigree / GS\u20139',9),row('Pheno / Pedigree / GS\u20136',6),row('Pedigree / GS\u20134',4),row('Speed breeding GS\u20132',2),row('GS\u20136 | 50% renewal',6,renew=16),row('GS\u20136 | 75% renewal',6,renew=24),row('GS\u20136 | 100% renewal',6,renew=32))),
 list(title='Selection intensity',color='#B26F21',rows=list(row('Pedigree\u20136 | early1500',6,TRUE),row('GS\u20136 | early1500',6,TRUE),row('Speed breeding GS\u20132 | early1500',2,TRUE))),
 list(title='Optimal Contribution Selection',color='#337B65',rows=lapply(c('Diversity','Balanced','Gain'),function(p)row(paste0('GS\u20136 / speed breeding GS\u20132\nearly1500 | OCS ',p),c(6,2),TRUE,ocs=TRUE))))
draw<-function(){
 par(mar=rep(0,4),xaxs='i',yaxs='i',family='sans',bg='white');plot.new();plot.window(c(0,2000),c(2020,0),asp=NA)
 ink<-'#243C3B';muted<-'#778985';blue<-'#2366A3';orange<-'#DF7922';teal<-'#137B76';purple<-'#8060A8'
 box<-function(x,y,w,h,col)rect(x,y,x+w,y+h,col=col,border=NA)
 tx<-function(x,y,s,size=16,col=ink,bold=FALSE,anchor=.5)text(x,y,s,adj=c(anchor,.2),cex=size*.432/12,font=if(bold)2 else 1,col=col)
 ln<-function(x,y,x2,y2,col,lwd=1)segments(x,y,x2,y2,col=col,lwd=lwd*.6)
 ar<-function(x,y,x2,y2,col,lwd=1){ln(x,y,x2,y2,col,lwd);polygon(c(x2,x2-9,x2-9),c(y2,y2-3.5,y2+3.5),col=col,border=NA)}
 circle<-function(x,y,r,fill,border=NA){t<-seq(0,2*pi,length.out=40);polygon(x+r*cos(t),y+r*sin(t),col=fill,border=border,lwd=.6)}
 cx<-470+(0:11)*101
 for(g in list(list(0,3,'Establishment'),list(4,5,'Preliminary trials'),list(6,8,'Advanced trials'),list(9,11,'Elite trials'))){tx(mean(cx[c(g[[1]],g[[2]])+1]),35,g[[3]],18,bold=TRUE);ln(cx[g[[1]]+1]-40,49,cx[g[[2]]+1]+40,49,'#CCD7E0',1.5)}
 tx(34,83,'Scheme',22,bold=TRUE,anchor=0);tx(1812,83,'Parents',22,bold=TRUE)
 tx(419,122,'Year',15,muted,anchor=1);tx(419,156,'Phenotyping year',15,muted,anchor=1);tx(419,190,'h\u00b2',15,muted,anchor=1)
 codes<-c('F1','Seedlings','J2','J3','PT1','PT2','AT1','AT2','AT3','ET1','ET2','ET3')
 for(i in 0:11){tx(cx[i+1],83,codes[i+1],17,bold=TRUE);tx(cx[i+1],122,i,18);tx(cx[i+1],156,if(i>=4)i-3 else '\u2013',16,muted);tx(cx[i+1],190,if(i<4)'\u2013' else if(i<6)'.20' else if(i<9)'.40' else '.50',16,muted)}
 y<-225
 for(b in blocks){box(25,y,1950,42,b$color);tx(43,y+29,b$title,23,'white',TRUE,0);y<-y+68
 for(ri in seq_along(b$rows)){r<-b$rows[[ri]]
 if(b$title=='Generation interval'&&ri%in%c(1,5)){tx(40,y+12,if(ri==1)'Timing of parental entry' else 'Parental replacement rate',19,b$color,TRUE,0);y<-y+40}
 rh<-if(r$ocs)100 else 76;if(ri%%2==1)box(25,y-9,1950,rh-11,'#F3F6F8')
 labels<-strsplit(r$label,'\n',fixed=TRUE)[[1]]
 if(length(labels)==2){tx(40,y+17,labels[1],17,bold=TRUE,anchor=0);tx(40,y+43,labels[2],17,bold=TRUE,anchor=0)}else tx(40,y+33,r$label,18,bold=TRUE,anchor=0)
 for(age in 0:11){x<-cx[age+1];fill<-if(r$early&&age==1)orange else if(age==6)teal else if(age==9)purple else '#E3ECF1';box(x-39,y+11,78,27,fill);n<-if(r$early&&age==0)1500 else if(age<=5)500 else if(age<=8)25 else 8;tx(x,y+31,n,16,if(fill%in%c(orange,teal,purple))'white' else ink,TRUE);if(age<11)ar(x+42,y+25,x+54,y+25,'#CBD5DC',.8)}
 ar(cx[min(r$parent)+1]+12,y-1,1676,y-1,blue,1.3)
 for(p in r$parent){circle(cx[p+1],y-1,9,'white',blue);tx(cx[p+1],y+4,'P',12,blue,TRUE)}
 if(r$ocs){tx(1803,y+9,'GS6: \u226410% per parent',16,blue,TRUE);tx(1803,y+27,'SB2: \u226516 adult parents',16,blue,TRUE);for(j in 0:4){h<-c(13,26,19,31,15)[j+1];box(1745+j*22,y+63-h,14,h,blue)}}else{
 tx(1803,y+2,sprintf('%s new + %s retained',r$renew,32-r$renew),14,blue,TRUE)
 for(k in 0:31)circle(1755+(k%%8)*13,y+14+(k%/%8)*12,3,if(k<r$renew)blue else '#AAB6C1')}
 y<-y+rh
 };y<-y+20}
 for(k in 1:3){x<-c(40,590,1100)[k];box(x,y+4,20,20,c(orange,teal,purple)[k]);tx(x+31,y+20,c('Early selection: 1500 to 500','PT-to-AT cut: 500 to 25','AT-to-ET cut: 25 to 8')[k],16,muted,anchor=0)}
 circle(1690,y+14,9,'white',blue);tx(1690,y+19,'P',12,blue,TRUE);tx(1710,y+20,'First parental use',16,muted,anchor=0)
}
stem<-file.path(out,'almond_all_schemes_corrected')
png(paste0(stem,'.png'),width=7200,height=7272,res=600,type='cairo',pointsize=12);draw();dev.off()
svg(paste0(stem,'.svg'),width=12,height=12.12,pointsize=12);draw();dev.off()
png(file.path(out,'preview.png'),width=1600,height=1616,res=133.333,type='cairo',pointsize=12);draw();dev.off()
cat('Generated PNG 7200x7272 at 600 dpi and SVG from the same R drawing function.\n')

