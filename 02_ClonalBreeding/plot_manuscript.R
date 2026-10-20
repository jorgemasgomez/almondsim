# Manuscript figures: horizontal trajectory triptychs and a vertical endpoint grid.
# Postprocessing only; endpoint statistical models do not fit breeding predictions or run simulations.
manuscript_catalogue <- function() {
 ids<-c('Pheno_Y9','Pedigree_Y9','GS_Y9','Pheno_Y6','Pedigree_Y6','GS_Y6',
  'Pedigree_Y4','GS_Y4','GS_SpeedBreeding_Y2',
  'Pedigree_Y6_Early1500','GS_Y6_Early1500','GS_SpeedBreeding_Y2_Early1500',
  paste0('GS_Y6_Early1500_OCS_',c('Diversity','Balanced','Gain')),
  paste0('GS_SpeedBreeding_Y2_Early1500_OCS_',c('Diversity','Balanced','Gain')))
 labels<-c('PH9','PED9','GS9','PH6','PED6','GS6','PED4','GS4','SB2',
  'PED6-E','GS6-E','SB2-E','GS6-OD','GS6-OB','GS6-OG','SB2-OD','SB2-OB','SB2-OG')
 block<-c(rep('Selection accuracy',3),rep('Generation interval',6),
  rep('Selection intensity',3),rep('Optimal contribution selection',6))
 family<-ifelse(grepl('^Pheno',ids),'Pheno',ifelse(grepl('^Pedigree',ids),'Pedigree','GS'))
 colors<-c('#B45309','#581C87','#102A56','#F59E0B','#9333EA','#2459A6',
  '#C084FC','#4A87D7','#86B7EA','#7E22CE','#123F72','#3294C6',
  '#111D56','#244DC7','#538DD5','#145A80','#478CAD','#9ACBE2')
 data.frame(scenario=ids,label=labels,block=block,family=family,color=colors,
  paper=TRUE,stringsAsFactors=FALSE)
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
   ids=c('GS_Y6_Early1500',cat$scenario[13:15],'GS_SpeedBreeding_Y2_Early1500',cat$scenario[16:18])))
}

# Order each legend top-to-bottom within meaningful, aligned columns.
manuscript_legend <- function(block) {
 ids<-block$ids
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

manuscript_metric_labels <- function()c('Selection accuracy','Genetic gain','Genetic variance')

manuscript_summary <- function(data) {
 keys<-unique(data[c('scenario','future_year')]);rows<-list()
 for(i in seq_len(nrow(keys))) {
  z<-data[data$scenario==keys$scenario[i]&data$future_year==keys$future_year[i],,drop=FALSE]
  for(j in seq_along(c('accSel','gain','variance'))) {
   v<-z[[c('accSel','gain','variance')[j]]];v<-v[is.finite(v)];n<-length(v)
   value<-if(n)mean(v)else NA_real_
   half<-if(n>1L)qt(.975,n-1L)*sd(v)/sqrt(n)else NA_real_
   rows[[length(rows)+1L]]<-data.frame(scenario=keys$scenario[i],year=keys$future_year[i],
    metric=manuscript_metric_labels()[j],value=value,lower=value-half,upper=value+half,n=n)
  }
 }
 do.call(rbind,rows)
}

# Endpoint inference: one observation per paired replicate and scheme.
# Never infer from partial coverage or silently reduce the comparison family.
manuscript_endpoint_tests <- function(data, ids, expected_reps, out, alpha=.05) {
 empty<-data.frame()
 results<-list(letters=empty,pairs=empty,means=empty,diagnostics=empty)
 status<-list();fits<-list()
 for(metric in unique(as.character(data$metric))) {
  z<-data[as.character(data$metric)==metric&is.finite(data$value),c('rep','scenario','value')]
  z$scenario<-as.character(z$scenario)
  complete<-!anyDuplicated(z[c('rep','scenario')])&&length(expected_reps)>=3L&&
   nrow(z)==length(ids)*length(expected_reps)&&
   all(vapply(ids,function(id)setequal(z$rep[z$scenario==id],expected_reps),logical(1)))
  if(!complete) {
   status[[metric]]<-data.frame(metric=metric,status='pending_complete_paired_coverage',
    detail='Letters require every requested scheme in every planned replicate (minimum three).')
   next
  }
  packages<-c('nlme','emmeans','multcomp','multcompView')
  missing<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) {
   status[[metric]]<-data.frame(metric=metric,status='missing_package',detail=paste(missing,collapse=', '));next
  }
  answer<-tryCatch({
   z$scenario<-factor(z$scenario,levels=ids);z$rep<-factor(z$rep)
   # Relative scaling improves numerical conditioning without changing inference.
   scale<-sd(z$value)
   if(!is.finite(scale)||scale<=0)stop('Endpoint has no variation.')
   z$response<-z$value/scale
   fit<-nlme::lme(response~scenario,random=~1|rep,data=z,method='REML',na.action=na.fail)
   emm<-emmeans::emmeans(fit,~scenario)
   comparisons<-as.data.frame(emmeans::contrast(emm,method='pairwise',adjust='tukey'))
   means<-as.data.frame(emmeans::emmeans(fit,~scenario))
   letters<-as.data.frame(multcomp::cld(emm,alpha=alpha,Letters=base::letters,
     adjust='tukey',sort=TRUE,reversed=TRUE))
   # CLD uses Tukey-adjusted pairwise tests; mean confidence intervals are separate.
   letters$letter<-gsub(' ', '',letters$.group,fixed=TRUE)
   for(nm in intersect(c('emmean','SE','lower.CL','upper.CL'),names(means)))means[[nm]]<-means[[nm]]*scale
   for(nm in intersect(c('estimate','SE'),names(comparisons)))comparisons[[nm]]<-comparisons[[nm]]*scale
   diagnostic<-data.frame(rep=z$rep,scenario=z$scenario,observed=z$value,
     fitted=as.numeric(fitted(fit))*scale,normalized_residual=as.numeric(residuals(fit,type='normalized')))
   list(fit=fit,means=means,pairs=comparisons,
     letters=letters[c('scenario','letter')],diagnostics=diagnostic)
  },error=function(e)e)
  if(inherits(answer,'error')) {
   status[[metric]]<-data.frame(metric=metric,status='analysis_failed',detail=conditionMessage(answer));next
  }
  for(key in names(results)) {
   answer[[key]]$metric<-metric
   results[[key]]<-rbind(results[[key]],answer[[key]])
  }
  fits[[metric]]<-answer$fit
  status[[metric]]<-data.frame(metric=metric,status='complete',
   detail='REML random replicate intercept; fixed scheme; global Tukey pairwise tests, alpha=0.05; common residual variance assumed.')
 }
 results$status<-do.call(rbind,status)
 for(key in names(results))write.csv(results[[key]],file.path(out,paste0('endpoint_tests_',key,'.csv')),row.names=FALSE)
 saveRDS(fits,file.path(out,'endpoint_tests_models.rds'))
 writeLines(c('Each metric is a separate comparison family including every main-text scheme.',
  'Model: endpoint ~ scheme + (1 | replicate), REML; Tukey multiplicity correction, alpha=0.05.',
  'Raw outcome scale; common residual variance and approximately normal conditional errors assumed.',
  'Review exported residuals and model objects before publication; no automatic assumption test selects the model.',
  'Shared letters mean no significant difference was detected, not equivalence.',
  'Boxplots retain all available observations. Letters appear only with full planned paired coverage.'),
  file.path(out,'endpoint_tests_method.txt'))
 results
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
 summaries<-list();panel_coverage<-list();palette_key<-list()
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
  labs<-setNames(cat$label[match(b$ids,cat$scenario)],b$ids)
  colors<-setNames(cat$color[match(b$ids,cat$scenario)],b$ids)
  if(b$name=='04_optimal_contribution_selection') {
   # Separate GS6 and speed breeding by hue in this OCS-only comparison.
   colors[c('GS_Y6_Early1500',paste0('GS_Y6_Early1500_OCS_',c('Diversity','Balanced','Gain')))]<-
    c('#123F72','#172E6B','#356CC0','#7BA5D8')
   colors[c('GS_SpeedBreeding_Y2_Early1500',paste0('GS_SpeedBreeding_Y2_Early1500_OCS_',c('Diversity','Balanced','Gain')))]<-
    c('#00574F','#168979','#40B6A5','#86D5C6')
  }
  palette_key[[b$name]]<-data.frame(panel=b$name,scenario=b$ids,color=unname(colors[b$ids]))
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
 write.csv(do.call(rbind,palette_key),file.path(out,'panel_color_key.csv'),row.names=FALSE)
 write.csv(do.call(rbind,panel_coverage),file.path(out,'panel_coverage.csv'),row.names=FALSE)
 saveRDS(summaries,file.path(out,'trajectory_summaries.rds'))
 saveRDS(input$data,file.path(out,'complete_trajectory_data.rds'))

 # Last-year boxplots: all 18 main-text schemes, ordered by Methods block.
 # Each point is one replicate from one of the 18 retained schemes.
 paper<-cat[cat$paper,,drop=FALSE]
 final<-input$data[input$data$future_year==input$config$config$future&input$data$scenario%in%paper$scenario,,drop=FALSE]
 endpoint_labels<-c('Genetic gain','Genetic variance')
 rows<-lapply(seq_along(c('gain','variance')),function(j) {
  z<-final[c('rep','scenario')];z$value<-final[[c('gain','variance')[j]]]
  z$metric<-endpoint_labels[j];z
 })
 final_long<-do.call(rbind,rows);final_long$scenario<-factor(final_long$scenario,levels=paper$scenario)
 final_long$metric<-factor(final_long$metric,levels=endpoint_labels)
 final_long$block<-paper$block[match(final_long$scenario,paper$scenario)]
 final_long$x<-match(final_long$scenario,paper$scenario)
 inference<-manuscript_endpoint_tests(final_long,paper$scenario,
   seq_len(input$config$config$reps),out)
 annotations<-inference$letters
 if(nrow(annotations)) {
  annotations$x<-match(annotations$scenario,paper$scenario)
  annotations$metric<-factor(annotations$metric,levels=endpoint_labels)
  annotations$y<-vapply(seq_len(nrow(annotations)),function(i) {
   values<-final_long$value[final_long$metric==annotations$metric[i]]
   values<-values[is.finite(values)];span<-max(diff(range(values)),1)
   max(final_long$value[final_long$scenario==annotations$scenario[i]&
     final_long$metric==annotations$metric[i]],na.rm=TRUE)+.05*span
  },numeric(1))
 }
 blockcolors<-c('Selection accuracy'='#287DAD','Generation interval'='#556BA0',
   'Selection intensity'='#B16E20','Optimal contribution selection'='#2F8069')
 bands<-data.frame(block=names(blockcolors),xmin=c(.5,3.5,9.5,12.5),xmax=c(3.5,9.5,12.5,18.5))
 counts<-expand.grid(scenario=paper$scenario,metric=endpoint_labels,stringsAsFactors=FALSE)
 counts$n<-mapply(function(id,m)sum(final_long$scenario==id&as.character(final_long$metric)==m&is.finite(final_long$value)),counts$scenario,counts$metric)
 counts$label<-ifelse(counts$n==0,'n=0',paste0('n=',counts$n))
 counts$x<-match(counts$scenario,paper$scenario)
 counts$metric<-factor(counts$metric,levels=endpoint_labels)
 finite<-final_long[is.finite(final_long$value),,drop=FALSE]
 headers<-bands;headers$x<-(headers$xmin+headers$xmax)/2
 headers$metric<-factor(endpoint_labels[1],levels=endpoint_labels)
 q<-ggplot2::ggplot(finite,ggplot2::aes(x,value,group=scenario,fill=block))+
  ggplot2::geom_rect(data=bands,ggplot2::aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf,fill=block),inherit.aes=FALSE,alpha=.055,color=NA)+
  ggplot2::geom_vline(xintercept=c(3.5,9.5,12.5),color='#C9D0D5',linewidth=.4)+
  ggplot2::geom_boxplot(width=.58,outlier.shape=NA,linewidth=.45,alpha=.85,na.rm=TRUE)+
  ggplot2::geom_point(position=ggplot2::position_jitter(width=.10,height=0,seed=20261015),size=1.4,shape=21,color='#263238',stroke=.25,alpha=.8)+
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
 if(nrow(annotations))q<-q+
  ggplot2::geom_text(data=annotations,ggplot2::aes(x=x,y=y,label=letter),
    inherit.aes=FALSE,size=3.5,fontface='bold',vjust=0)+
  ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(.05,.18)))
 manuscript_save(q,file.path(out,'05_final_year_boxplots'),18,9,dpi)
 write.csv(final_long,file.path(out,'final_year_values.csv'),row.names=FALSE)
 write.csv(counts,file.path(out,'final_year_counts.csv'),row.names=FALSE)
 saveRDS(q,file.path(out,'final_year_boxplot.rds'))
 writeLines(c(paste('Data status:',subtitle),
  paste('Historical OCS:',historical_ocs),
  'No overall titles, subtitles or explanatory captions are printed inside figures.',
  'Four horizontal triptychs and one vertical endpoint grid.',
  'Trajectory metric order: selection accuracy, genetic gain, genetic variance. Solid lines only.',
  'Trajectory colours identify method families (Pheno orange, pedigree purple, GS blue); OCS panel separates GS6 blues from speed-breeding turquoises.',
  'Endpoint colours identify Methods blocks, in manuscript order. All 18 schemes are included in the endpoint grid.',
  'Paired trajectory means require shared complete replicate IDs; previews omit unavailable arms from pairing.',
  'Endpoint grid contains genetic gain and genetic variance only. Boxes show complete replicates per scheme, not paired differences.',
  'Gain is relative to the replicate burn-in; variance uses retained seedlings; accuracy is measured before early selection.',
  'Pheno accuracy is NA, not zero. n=0 is missing; n=1 is a single replicate, not an estimate of uncertainty.',
  paste('Source design:',version),'Historical OCS and native-offspring OCS must not be pooled.'),file.path(out,'README.txt'))
 invisible(out)
}

