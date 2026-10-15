# Manuscript figures: horizontal trajectory triptychs and a vertical endpoint grid.
# Plotting only; no simulation, model fit or OCS allocation is called here.
manuscript_catalogue <- function() {
 ids<-c('Pheno_Y9','Pedigree_Y9','GS_Y9','Pheno_Y6','Pedigree_Y6','GS_Y6',
  'Pedigree_Y4','GS_Y4','GS_SpeedBreeding_Y2',
  'Pedigree_Y6_Early1500','GS_Y6_Early1500','GS_SpeedBreeding_Y2_Early1500',
  paste0('GS_Y6_Early1500_OCS_',c('Diversity','Balanced','Gain')),
  paste0('GS_SpeedBreeding_Y2_Early1500_OCS_',c('Diversity','Balanced','Gain')),
  paste0('GS_Y6_R',c(50,75,100)))
 labels<-c('PH9','PED9','GS9','PH6','PED6','GS6','PED4','GS4','SB2',
  'PED6-E','GS6-E','SB2-E','GS6-OD','GS6-OB','GS6-OG','SB2-OD','SB2-OB','SB2-OG',
  'GS6-R50','GS6-R75','GS6-R100')
 block<-c(rep('Selection accuracy',3),rep('Generation interval',6),
  rep('Selection intensity',3),rep('Optimal contribution selection',6),rep('Parental renewal (supplement)',3))
 family<-ifelse(grepl('^Pheno',ids),'Pheno',ifelse(grepl('^Pedigree',ids),'Pedigree','GS'))
 # Method families are stable across all trajectory figures: orange, purple, blue.
 colors<-c('#B45309','#581C87','#102A56','#F59E0B','#9333EA','#2459A6',
  '#C084FC','#4A87D7','#86B7EA','#7E22CE','#123F72','#3294C6',
  '#111D56','#244DC7','#538DD5','#145A80','#478CAD','#9ACBE2',
  '#123E85','#5086C3','#A3C7EA')
 data.frame(scenario=ids,label=labels,block=block,family=family,color=colors,
  paper=seq_along(ids)<=18L,stringsAsFactors=FALSE)
}

manuscript_blocks <- function() {
 cat<-manuscript_catalogue()
 list(
  list(name='01_selection_accuracy',title='Selection method',ids=cat$scenario[1:3]),
  list(name='02_generation_interval',title='Generation interval',ids=cat$scenario[1:9]),
  list(name='03_selection_intensity',title='Selection intensity',
   ids=c('GS_Y6','GS_Y6_Early1500','Pedigree_Y6','Pedigree_Y6_Early1500',
     'GS_SpeedBreeding_Y2','GS_SpeedBreeding_Y2_Early1500')),
  list(name='04_optimal_contribution_selection',title='Optimal contribution selection',
   ids=c('GS_Y6_Early1500',cat$scenario[13:15],'GS_SpeedBreeding_Y2_Early1500',cat$scenario[16:18])),
  list(name='S01_parental_renewal',title='Parental renewal | supplementary',
   ids=c('GS_Y6',paste0('GS_Y6_R',c(50,75,100)))))
}

# Order each legend top-to-bottom within meaningful, aligned columns.
manuscript_legend <- function(block) {
 ids<-block$ids
 if(block$name=='S01_parental_renewal')return(list(order=ids,nrow=1L,padding=character()))
 if(block$name=='04_optimal_contribution_selection') {
  groups<-list(GS6=ids[grepl('^GS_Y6',ids)],SB2=ids[grepl('^GS_SpeedBreeding',ids)])
 }else if(block$name=='03_selection_intensity') {
  groups<-list(Pedigree=ids[grepl('^Pedigree',ids)],GS6=ids[grepl('^GS_Y6',ids)],
   SB2=ids[grepl('^GS_SpeedBreeding',ids)])
 }else {
  groups<-list(Pheno=ids[grepl('^Pheno',ids)],Pedigree=ids[grepl('^Pedigree',ids)],GS=ids[grepl('^GS',ids)])
 }
 groups<-groups[lengths(groups)>0]
 rows<-max(lengths(groups));order<-character();padding<-character()
 for(i in seq_along(groups)) {
  empty<-if(length(groups[[i]])<rows)paste0('.legend_blank_',i,'_',seq_len(rows-length(groups[[i]])))else character()
  order<-c(order,groups[[i]],empty);padding<-c(padding,empty)
 }
 stopifnot(setequal(setdiff(order,padding),ids))
 list(order=order,nrow=rows,padding=padding)
}

manuscript_read_results <- function(run) {
 cfg<-readRDS(file.path(run,'configuration.rds'))
 frames<-list();coverage<-list();catalogue<-manuscript_catalogue()
 ids<-intersect(catalogue$scenario,cfg$scenarios$scenario_id)
 for(rep in seq_len(cfg$config$reps)) {
  dir<-file.path(run,sprintf('rep_%03d',rep))
  baseline<-tryCatch(readRDS(file.path(dir,'burnin_results.rds')),error=function(e)NULL)
  if(!is.null(baseline)&&nrow(baseline))baseline<-baseline[which.max(baseline$year),,drop=FALSE]
  for(id in ids) {
   file<-file.path(dir,paste0('results_',id,'.rds'))
   x<-if(file.exists(file))tryCatch(readRDS(file),error=function(e)NULL)else NULL
   complete<-!is.null(baseline)&&!is.null(x)&&nrow(x)==cfg$config$future&&
    all(x$rep==rep)&&all(x$scenario==id)&&!anyDuplicated(x$year)&&
    setequal(x$year,cfg$config$burnin+seq_len(cfg$config$future))
   reason<-if(complete)'complete'else if(!file.exists(file))'not available'else 'incomplete or unreadable'
   if(complete&&'designVersion'%in%names(x)&&any(x$designVersion!=cfg$config$design_version)) {
    complete<-FALSE;reason<-'design version mismatch'
   }
   coverage[[length(coverage)+1L]]<-data.frame(rep=rep,scenario=id,complete=complete,reason=reason)
   if(!complete)next
   x$gain<-x$meanG_Seedlings-baseline$meanG_Seedlings
   x$variance<-x$varG_Seedlings;x$future_year<-x$year-cfg$config$burnin
   frames[[length(frames)+1L]]<-x
  }
 }
 list(config=cfg,data=if(length(frames))do.call(rbind,frames)else data.frame(),
  coverage=do.call(rbind,coverage),catalogue=catalogue)
}

manuscript_metric_labels <- function()c('Genetic gain','Genetic variance','Selection accuracy')

manuscript_summary <- function(data) {
 keys<-unique(data[c('scenario','future_year')]);rows<-list()
 for(i in seq_len(nrow(keys))) {
  z<-data[data$scenario==keys$scenario[i]&data$future_year==keys$future_year[i],,drop=FALSE]
  for(j in seq_along(c('gain','variance','accSel'))) {
   v<-z[[c('gain','variance','accSel')[j]]];v<-v[is.finite(v)];n<-length(v)
   value<-if(n)mean(v)else NA_real_
   half<-if(n>1L)qt(.975,n-1L)*sd(v)/sqrt(n)else NA_real_
   rows[[length(rows)+1L]]<-data.frame(scenario=keys$scenario[i],year=keys$future_year[i],
    metric=manuscript_metric_labels()[j],value=value,lower=value-half,upper=value+half,n=n)
  }
 }
 do.call(rbind,rows)
}

manuscript_save <- function(plot,path,width,height,dpi) {
 ggplot2::ggsave(paste0(path,'.png'),plot,width=width,height=height,dpi=dpi,bg='white',limitsize=FALSE)
 ggplot2::ggsave(paste0(path,'.pdf'),plot,width=width,height=height,device=grDevices::cairo_pdf,bg='white',limitsize=FALSE)
}

plot_manuscript <- function(run,with_intervals=TRUE,out=NULL,preview=FALSE,paired=TRUE,dpi=400) {
 if(!requireNamespace('ggplot2',quietly=TRUE))stop('Plotting requires ggplot2.')
 run<-normalizePath(run,winslash='/',mustWork=TRUE)
 input<-manuscript_read_results(run);cat<-input$catalogue
 if(is.null(out))out<-file.path(run,'manuscript_figures',if(with_intervals)'with_CI'else'lines')
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 write.csv(input$coverage,file.path(out,'coverage.csv'),row.names=FALSE)
 write.csv(cat,file.path(out,'scheme_key.csv'),row.names=FALSE)
 if(!nrow(input$data)) {writeLines('No complete trajectories available.',file.path(out,'README.txt'));return(invisible(out))}
 stopifnot(!anyDuplicated(input$data[c('rep','scenario','year')]))
 all_complete<-all(input$coverage$complete)
 version<-input$config$config$design_version
 historical_ocs<-!identical(version,'PT_AT_ET_optiSel_native_offspring_v9')
 subtitle<-if(preview||!all_complete)'PROVISIONAL | complete trajectories only; unavailable schemes remain explicit'else 'Complete trajectories'
 summaries<-list();panel_coverage<-list()
 for(b in manuscript_blocks()) {
  selected<-input$data[input$data$scenario%in%b$ids,,drop=FALSE]
  available<-b$ids[b$ids%in%unique(selected$scenario)]
  if(paired&&length(available)) {
   required<-if(preview)available else b$ids
   sets<-lapply(required,function(id)unique(selected$rep[selected$scenario==id]))
   common<-Reduce(intersect,sets)
   selected<-selected[selected$rep%in%common,,drop=FALSE]
  }
  counts<-setNames(vapply(b$ids,function(id)length(unique(selected$rep[selected$scenario==id])),integer(1)),b$ids)
  panel_coverage[[b$name]]<-data.frame(panel=b$name,scenario=b$ids,n=as.integer(counts),
   paired=paired,preview=preview)
  labs<-setNames(paste0(cat$label[match(b$ids,cat$scenario)],' (n=',counts,')'),b$ids)
  if(b$name=='S01_parental_renewal')labs['GS_Y6']<-paste0('GS6-R25 (n=',counts['GS_Y6'],')')
  colors<-setNames(cat$color[match(b$ids,cat$scenario)],b$ids)
  legend<-manuscript_legend(b)
  colors<-c(colors,setNames(rep('transparent',length(legend$padding)),legend$padding))
  labs<-c(labs,setNames(rep('',length(legend$padding)),legend$padding))
  if(nrow(selected))summary<-manuscript_summary(selected)else summary<-data.frame(
   scenario=character(),year=numeric(),metric=character(),value=numeric(),lower=numeric(),upper=numeric(),n=integer())
  summary$scenario<-factor(summary$scenario,levels=b$ids)
  summary$metric<-factor(summary$metric,levels=manuscript_metric_labels())
  # Keep all three panels even when an entire requested method has no results yet.
  skeleton<-expand.grid(scenario=b$ids,year=1,metric=manuscript_metric_labels(),stringsAsFactors=FALSE)
  skeleton$metric<-factor(skeleton$metric,levels=manuscript_metric_labels());skeleton$value<-NA_real_
  p<-ggplot2::ggplot(summary,ggplot2::aes(year,value,color=scenario,fill=scenario,group=scenario))+
   ggplot2::geom_blank(data=skeleton,na.rm=TRUE)
  if(with_intervals)p<-p+ggplot2::geom_ribbon(ggplot2::aes(ymin=lower,ymax=upper),alpha=.10,color=NA,na.rm=TRUE)
  p<-p+ggplot2::geom_line(linewidth=.9,na.rm=TRUE,show.legend=TRUE)+
   ggplot2::facet_wrap(~metric,nrow=1,scales='free_y',drop=FALSE)+
   ggplot2::scale_color_manual(values=colors,limits=legend$order,breaks=legend$order,labels=labs,drop=FALSE)+
   ggplot2::scale_fill_manual(values=colors,limits=legend$order,breaks=legend$order,labels=labs,drop=FALSE)+
   ggplot2::scale_x_continuous(limits=c(1,input$config$config$future),breaks=c(1,10,20,30,40))+
   ggplot2::labs(title=NULL,subtitle=NULL,caption=NULL,x='Future year',y=NULL,color=NULL,fill=NULL)+
   ggplot2::theme_bw(base_size=13)+ggplot2::theme(
    plot.title=ggplot2::element_blank(),plot.subtitle=ggplot2::element_blank(),
    strip.background=ggplot2::element_rect(fill='#F1F4F6',color='#CED5DB'),
    strip.text=ggplot2::element_text(face='bold',size=12),
    legend.position='bottom',legend.text=ggplot2::element_text(size=10),
    panel.grid.minor=ggplot2::element_blank(),panel.spacing=grid::unit(1.1,'lines'),
    plot.caption=ggplot2::element_blank(),plot.margin=ggplot2::margin(10,16,10,12))+
   ggplot2::guides(color=ggplot2::guide_legend(nrow=legend$nrow,byrow=FALSE),
     fill=ggplot2::guide_legend(nrow=legend$nrow,byrow=FALSE))
  manuscript_save(p,file.path(out,b$name),18,6,dpi)
  summaries[[b$name]]<-summary
 }
 write.csv(do.call(rbind,panel_coverage),file.path(out,'panel_coverage.csv'),row.names=FALSE)
 saveRDS(summaries,file.path(out,'trajectory_summaries.rds'))
 saveRDS(input$data,file.path(out,'complete_trajectory_data.rds'))

 # Last-year boxplots: all 18 main-text schemes, ordered by Methods block.
 # Renewal and exploratory GS4 OCS are excluded. Each point is one replicate.
 paper<-cat[cat$paper,,drop=FALSE]
 final<-input$data[input$data$future_year==input$config$config$future&input$data$scenario%in%paper$scenario,,drop=FALSE]
 rows<-lapply(seq_along(c('gain','variance','accSel')),function(j) {
  z<-final[c('rep','scenario')];z$value<-final[[c('gain','variance','accSel')[j]]]
  z$metric<-manuscript_metric_labels()[j];z
 })
 final_long<-do.call(rbind,rows);final_long$scenario<-factor(final_long$scenario,levels=paper$scenario)
 final_long$metric<-factor(final_long$metric,levels=manuscript_metric_labels())
 final_long$block<-paper$block[match(final_long$scenario,paper$scenario)]
 final_long$x<-match(final_long$scenario,paper$scenario)
 blockcolors<-c('Selection accuracy'='#287DAD','Generation interval'='#556BA0',
   'Selection intensity'='#B16E20','Optimal contribution selection'='#2F8069')
 bands<-data.frame(block=names(blockcolors),xmin=c(.5,3.5,9.5,12.5),xmax=c(3.5,9.5,12.5,18.5))
 counts<-expand.grid(scenario=paper$scenario,metric=manuscript_metric_labels(),stringsAsFactors=FALSE)
 counts$n<-mapply(function(id,m)sum(final_long$scenario==id&as.character(final_long$metric)==m&is.finite(final_long$value)),counts$scenario,counts$metric)
 counts$label<-ifelse(counts$n==0,'n=0',paste0('n=',counts$n))
 counts$label[counts$scenario%in%c('Pheno_Y9','Pheno_Y6')&counts$metric==manuscript_metric_labels()[3]]<-'NA'
 counts$x<-match(counts$scenario,paper$scenario)
 counts$metric<-factor(counts$metric,levels=manuscript_metric_labels())
 finite<-final_long[is.finite(final_long$value),,drop=FALSE]
 headers<-bands;headers$x<-(headers$xmin+headers$xmax)/2
 headers$metric<-factor(manuscript_metric_labels()[1],levels=manuscript_metric_labels())
 q<-ggplot2::ggplot(finite,ggplot2::aes(x,value,group=scenario,fill=block))+
  ggplot2::geom_rect(data=bands,ggplot2::aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf,fill=block),inherit.aes=FALSE,alpha=.055,color=NA)+
  ggplot2::geom_vline(xintercept=c(3.5,9.5,12.5),color='#C9D0D5',linewidth=.4)+
  ggplot2::geom_boxplot(width=.58,outlier.shape=NA,linewidth=.45,alpha=.85,na.rm=TRUE)+
  ggplot2::geom_point(position=ggplot2::position_jitter(width=.10,height=0,seed=20261015),size=1.4,shape=21,color='#263238',stroke=.25,alpha=.8)+
  ggplot2::geom_text(data=counts,ggplot2::aes(x=x,y=Inf,label=label),inherit.aes=FALSE,vjust=1.5,size=2.8,color='#626A70')+
  ggplot2::geom_text(data=headers,ggplot2::aes(x=x,y=Inf,label=block,color=block),inherit.aes=FALSE,vjust=-.8,size=3.3,fontface='bold')+
  ggplot2::facet_wrap(~metric,ncol=1,scales='free_y',strip.position='left',drop=FALSE)+
  ggplot2::scale_x_continuous(breaks=seq_len(nrow(paper)),labels=paper$label,limits=c(.5,nrow(paper)+.5),expand=c(0,0))+
  ggplot2::scale_fill_manual(values=blockcolors,drop=FALSE)+ggplot2::scale_color_manual(values=blockcolors,guide='none')+
  ggplot2::coord_cartesian(clip='off')+
  ggplot2::labs(title=NULL,subtitle=NULL,caption=NULL,x=NULL,y=NULL,fill=NULL)+
  ggplot2::theme_bw(base_size=12)+ggplot2::theme(
   plot.title=ggplot2::element_blank(),plot.subtitle=ggplot2::element_blank(),
   axis.text.x=ggplot2::element_text(angle=45,hjust=1,size=10,color='#263238'),
   strip.placement='outside',strip.background=ggplot2::element_blank(),
   strip.text.y.left=ggplot2::element_text(angle=90,face='bold',size=11),
   panel.grid.major.x=ggplot2::element_blank(),panel.grid.minor=ggplot2::element_blank(),
   panel.spacing=grid::unit(1.8,'lines'),legend.position='none',
   plot.caption=ggplot2::element_blank(),plot.margin=ggplot2::margin(22,18,10,12))
 manuscript_save(q,file.path(out,'05_final_year_boxplots'),18,12.5,dpi)
 write.csv(final_long,file.path(out,'final_year_values.csv'),row.names=FALSE)
 write.csv(counts,file.path(out,'final_year_counts.csv'),row.names=FALSE)
 saveRDS(q,file.path(out,'final_year_boxplot.rds'))
 writeLines(c(paste('Data status:',subtitle),
  paste('Historical OCS:',historical_ocs),
  'No overall titles, subtitles or explanatory captions are printed inside figures.',
  'Four main horizontal triptychs, one supplementary renewal triptych, one vertical endpoint grid.',
  'Trajectory metric order: gain, genetic variance, prediction accuracy. Solid lines only.',
  'Trajectory colours identify method families (Pheno orange, pedigree purple, GS blue).',
  'Endpoint colours identify Methods blocks, in manuscript order. No renewal or GS4 OCS in endpoint grid.',
  'Paired trajectory means require shared complete replicate IDs; previews omit unavailable arms from pairing.',
  'Endpoint boxes show available complete replicates per scheme; they are not paired differences.',
  'Gain is relative to the replicate burn-in; variance uses retained seedlings; accuracy is measured before early selection.',
  'Pheno accuracy is NA, not zero. n=0 is missing; n=1 is a single replicate, not an estimate of uncertainty.',
  paste('Source design:',version),'Historical OCS and native-offspring OCS must not be pooled.'),file.path(out,'README.txt'))
 invisible(out)
}
