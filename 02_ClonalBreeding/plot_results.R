# Read-only result plotting. Only complete scenario/replicate trajectories enter means.
# Output lives under the external run, never beside source code.
plot_experiment <- function(run,with_intervals=TRUE) {
 if(!requireNamespace('ggplot2',quietly=TRUE)) stop('Plotting requires ggplot2.')
 config <- readRDS(file.path(run,'configuration.rds'))
 if(!config$config$design_version %in% c('PT_AT_ET_first_use_v2','PT_AT_ET_matched_cohorts_rrblup_v3','PT_AT_ET_annual_matched_PT_v4','PT_AT_ET_decision_specific_PT_v5','PT_AT_ET_Y9_repeated_calibrated_v6','PT_AT_ET_optiSel_profiles_v7','PT_AT_ET_optiSel_cap10_v8','PT_AT_ET_optiSel_native_offspring_v9')) stop('Use the historical plotter for old simulation runs.')
 out <- file.path(run,if(with_intervals) 'plots' else 'plots_lines')
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 frames <- list();coverage <- list()
 for(rep in seq_len(config$config$reps)) {
  directory <- file.path(run,sprintf('rep_%03d',rep))
  burn <- file.path(directory,'burnin_results.rds')
  baseline <- if(file.exists(burn)) readRDS(burn) else NULL
  if(!is.null(baseline)) baseline <- baseline[which.max(baseline$year),,drop=FALSE]
  for(id in config$scenarios$scenario_id) {
   file <- file.path(directory,paste0('results_',id,'.rds'))
   complete <- !is.null(baseline)&&file.exists(file)
   if(complete) {
    x <- readRDS(file)
    complete <- nrow(x)==config$config$future&&all(x$rep==rep)&&all(x$scenario==id)&&
      setequal(x$year,config$config$burnin+seq_len(config$config$future))&&!anyDuplicated(x$year)
   }
   coverage[[length(coverage)+1L]] <- data.frame(rep=rep,scenario=id,complete=complete)
   if(!complete) next
   x$gain <- x$meanG_Seedlings-baseline$meanG_Seedlings
   x$variance <- x$varG_Seedlings
   x$parent_gain <- x$meanG_Parents-baseline$meanG_Parents
   x$parent_variance <- x$varG_Parents
   x$future_year <- x$year-config$config$burnin
   frames[[length(frames)+1L]] <- x
  }
 }
 coverage <- redesign_bind(coverage)
 write.csv(coverage,file.path(out,'replicate_coverage.csv'),row.names=FALSE)
 data <- redesign_bind(frames)
 if(is.null(data)) {
  writeLines('No complete scenario/replicate trajectories available; see replicate_coverage.csv.',file.path(out,'README.txt'))
  return(invisible(out))
 }
 stopifnot(!anyDuplicated(data[c('rep','scenario','year')]))
 profile <- c('Diversity','Balanced','Gain')
 blocks <- list(
  list(name='01_selection_accuracy_seedlings',title='Selection accuracy | seedlings',
    ids=paste0(c('Pheno','Pedigree','GS'),'_Y9'),labels=c('Pheno-9','Pedigree-9','GS-9')),
  list(name='02_generation_interval_seedlings',title='First parental use and turnover | seedlings',
    ids=c(paste0(c('Pheno','Pedigree','GS'),'_Y9'),paste0(c('Pheno','Pedigree','GS'),'_Y6'),
          'Pedigree_Y4','GS_Y4','GS_SpeedBreeding_Y2',paste0('GS_Y6_R',c(50,75,100))),
    labels=c('Pheno-9','Pedigree-9','GS-9','Pheno-6','Pedigree-6','GS-6',
          'Pedigree-4','GS-4','Speed breeding GS-2','GS-6 renewal 50%','GS-6 renewal 75%','GS-6 renewal 100%')),
  list(name='03_selection_intensity_seedlings',title='Selection intensity | seedlings',
    ids=c('Pedigree_Y6','GS_Y6','GS_SpeedBreeding_Y2','Pedigree_Y6_Early1500','GS_Y6_Early1500','GS_SpeedBreeding_Y2_Early1500'),
    labels=c('Pedigree-6 | 500','GS-6 | 500','Speed GS-2 | 500','Pedigree-6 | 1500 to 500','GS-6 | 1500 to 500','Speed GS-2 | 1500 to 500')),
  list(name='04_ocs_seedlings',title='OCS | seedlings',
    ids=c('GS_Y4',paste0('GS_Y4_Early1500_OCS_',profile),'GS_Y6_Early1500',paste0('GS_Y6_Early1500_OCS_',profile),
          'GS_SpeedBreeding_Y2_Early1500',paste0('GS_SpeedBreeding_Y2_Early1500_OCS_',profile)),
    labels=c('GS-4 | 500, no OCS',paste0('GS-4 early1500 | ',profile),'GS-6 | no OCS',paste0('GS-6 | ',profile),'Speed GS-2 | no OCS',paste0('Speed GS-2 | ',profile))),
  list(name='05_speed_breeding_seedlings',title='Speed breeding | seedlings',
    ids=c('GS_Y6','GS_Y4','GS_SpeedBreeding_Y2','GS_Y6_Early1500','GS_SpeedBreeding_Y2_Early1500'),
    labels=c('GS-6 | 500','GS-4 | 500','Speed GS-2 | 500','GS-6 | 1500 to 500','Speed GS-2 | 1500 to 500')))
 parents <- lapply(blocks[1:2],function(b) {
  b$name <- sub('seedlings','parents',b$name);b$title <- sub('seedlings','parents',b$title);b$parents <- TRUE;b
 })
 blocks <- c(blocks,parents);summaries <- list()
 for(b in blocks) {
  selected <- data[data$scenario%in%b$ids,,drop=FALSE]
  if(!nrow(selected)) next
  metrics <- if(isTRUE(b$parents)) c('accParents','parent_gain','parent_variance') else c('accSel','gain','variance')
  keys <- unique(selected[c('scenario','future_year')])
  summary <- do.call(rbind,lapply(seq_len(nrow(keys)),function(i) {
   z <- selected[selected$scenario==keys$scenario[i]&selected$future_year==keys$future_year[i],,drop=FALSE]
   do.call(rbind,lapply(seq_along(metrics),function(j) {
    values <- z[[metrics[j]]];values <- values[is.finite(values)]
    n <- length(values);value <- if(n) mean(values) else NA_real_
    interval <- if(n>1L) qt(.975,n-1L)*sd(values)/sqrt(n) else NA_real_
    data.frame(scenario=keys$scenario[i],year=keys$future_year[i],metric=c('Selection accuracy','Genetic gain','Genetic variance')[j],value=value,lower=value-interval,upper=value+interval,n=n)
   }))
  }))
  counts <- vapply(b$ids,function(id) length(unique(selected$rep[selected$scenario==id])),integer(1))
  labels <- paste0(b$labels,' (n=',counts,')')
  summary$Scheme <- factor(labels[match(summary$scenario,b$ids)],levels=labels[counts>0])
  summary$metric <- factor(summary$metric,levels=c('Selection accuracy','Genetic gain','Genetic variance'))
  palette <- c(GS_Y4_Early1500_OCS_Diversity='#86B66B',GS_Y4_Early1500_OCS_Balanced='#D9B85A',GS_Y4_Early1500_OCS_Gain='#F28B82',Pheno_Y9='#C2410C',Pheno_Y6='#FB923C',Pedigree_Y9='#581C87',Pedigree_Y6='#9333EA',Pedigree_Y4='#C084FC',Pedigree_Y6_Early1500='#7E22CE',GS_Y9='#172554',GS_Y6='#1D4ED8',GS_Y4='#60A5FA',GS_Y6_Early1500='#0284C7',GS_Y6_R50='#2563EB',GS_Y6_R75='#0891B2',GS_Y6_R100='#67C1D4',GS_SpeedBreeding_Y2='#115E59',GS_SpeedBreeding_Y2_Early1500='#14B8A6',GS_Y6_Early1500_OCS_Diversity='#166534',GS_Y6_Early1500_OCS_Balanced='#CA8A04',GS_Y6_Early1500_OCS_Gain='#DC2626',GS_SpeedBreeding_Y2_Early1500_OCS_Diversity='#4D7C0F',GS_SpeedBreeding_Y2_Early1500_OCS_Balanced='#F59E0B',GS_SpeedBreeding_Y2_Early1500_OCS_Gain='#FB7185')
  colors <- setNames(unname(palette[b$ids]),labels)
  p <- ggplot2::ggplot(summary,ggplot2::aes(year,value,color=Scheme,fill=Scheme,group=Scheme))
  if(with_intervals) p <- p+ggplot2::geom_ribbon(ggplot2::aes(ymin=lower,ymax=upper),alpha=.10,color=NA,na.rm=TRUE)
  p <- p+ggplot2::geom_line(linewidth=.8,na.rm=TRUE)+ggplot2::facet_wrap(~metric,ncol=1,scales='free_y')+
   ggplot2::scale_color_manual(values=colors)+ggplot2::scale_fill_manual(values=colors)+
   ggplot2::labs(title=b$title,x='Future year',y=NULL,color=NULL,fill=NULL)+ggplot2::theme_bw(base_size=12)+
   ggplot2::theme(legend.position='bottom',panel.grid.minor=ggplot2::element_blank())
  ggplot2::ggsave(file.path(out,paste0(b$name,'.png')),p,width=12,height=10,dpi=400)
  ggplot2::ggsave(file.path(out,paste0(b$name,'.pdf')),p,width=12,height=10)
  summaries[[b$name]] <- summary
 }
 saveRDS(data,file.path(out,'plot_data.rds'));saveRDS(summaries,file.path(out,'plot_summaries.rds'))
 final <- data[data$future_year==config$config$future,,drop=FALSE]
 write.csv(final[c('rep','scenario','gain','variance','accSel','parent_gain','parent_variance')],file.path(out,'final_values_by_replicate.csv'),row.names=FALSE)
 writeLines(c('Only complete scenario/replicate trajectories enter summaries. Coverage is explicit.',
  'Seedling accuracy is measured before early selection; gain and variance use retained seedlings.',
  'Pheno seedling accuracy is NA; no zero line is invented. Gains are relative to each replicate burn-in.',
  'Parent accuracy refers to candidates before renewal; parent gain and variance refer to the active crossing pool.',
  'Pointwise 95% Monte Carlo intervals quantify replicate variability, not model-assumption uncertainty.'),file.path(out,'README.txt'))
 invisible(out)
}
